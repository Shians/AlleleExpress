#include "bam_reader.h"
#include <stdexcept>
#include <string>
#include <utility>
#include <sstream>
#include <regex>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <unordered_set>

// Include the htslib headers directly in the implementation file
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

// BAM flags from htslib/sam.h
#define BAM_FPAIRED        1   // paired-end / multiple segment sequencing
#define BAM_FPROPER_PAIR   2   // each segment properly aligned
#define BAM_FUNMAP         4   // segment unmapped
#define BAM_FMUNMAP        8   // next segment unmapped
#define BAM_FREVERSE      16   // SEQ is reverse complemented
#define BAM_FMREVERSE     32   // SEQ of next segment reversed
#define BAM_FREAD1        64   // first segment in template
#define BAM_FREAD2       128   // last segment in template
#define BAM_FSECONDARY   256   // not primary alignment
#define BAM_FQCFAIL      512   // quality control failure
#define BAM_FDUP        1024   // PCR or optical duplicate
#define BAM_FSUPPLEMENTARY 2048 // supplementary alignment

BamReader::BamReader(const std::string& filename)
    : filename_(filename)
    , file_(nullptr, [](htsFile* f) { if (f) hts_close(f); })
    , header_(nullptr, [](sam_hdr_t* h) { if (h) sam_hdr_destroy(h); })
    , index_(nullptr, [](hts_idx_t* i) { if (i) hts_idx_destroy(i); })
    , iterator_(nullptr, [](hts_itr_t* itr) { if (itr) hts_itr_destroy(itr); })
    , read_(bam_init1(), [](bam1_t* r) { if (r) bam_destroy1(r); }) {

    // Open the BAM file
    file_.reset(hts_open(filename.c_str(), "r"));
    if (!file_) {
        throw std::runtime_error("Failed to open BAM file: " + filename);
    }

    // Check if the file is BAM
    if (file_->format.format != htsExactFormat::bam) {
        throw std::runtime_error("Input file must be a BAM or CRAM file");
    }

    // Read the header
    header_.reset(sam_hdr_read(file_.get()));
    if (!header_) {
        throw std::runtime_error("Failed to read header from BAM file");
    }

    // Load the index
    index_.reset(sam_index_load(file_.get(), filename.c_str()));
    if (!index_) {
        throw std::runtime_error("Failed to load index for BAM file. Make sure the BAM file is indexed (e.g., with samtools index)");
    }

    // Initialize filter params with samtools mpileup default values
    filter_params_.min_base_quality = DEFAULT_MIN_BASE_QUALITY;
    filter_params_.min_mapping_quality = DEFAULT_MIN_MAPPING_QUALITY;
    filter_params_.max_depth = DEFAULT_MAX_DEPTH;
    filter_params_.exclude_flags = DEFAULT_EXCLUDE_FLAGS;
    filter_params_.required_flags = 0; // No flags required by default

    // By default, don't count orphans and handle overlapping pairs
    // These match samtools mpileup defaults
    filter_params_.set_flag(ReadFilterFlag::COUNT_ORPHANS, false);
    filter_params_.set_flag(ReadFilterFlag::DISABLE_OVERLAP, false);
}

BamReader::~BamReader() = default;

BamReader::BamReader(BamReader&& other) noexcept
    : filename_(std::move(other.filename_))
    , file_(std::move(other.file_))
    , header_(std::move(other.header_))
    , index_(std::move(other.index_))
    , iterator_(std::move(other.iterator_))
    , read_(std::move(other.read_))
    , has_more_reads_(other.has_more_reads_)
    , filter_params_(std::move(other.filter_params_)) {
}

BamReader& BamReader::operator=(BamReader&& other) noexcept {
    if (this != &other) {
        filename_ = std::move(other.filename_);
        file_ = std::move(other.file_);
        header_ = std::move(other.header_);
        index_ = std::move(other.index_);
        iterator_ = std::move(other.iterator_);
        read_ = std::move(other.read_);
        has_more_reads_ = other.has_more_reads_;
        filter_params_ = std::move(other.filter_params_);
    }
    return *this;
}

void BamReader::set_filter_params(const BamFilterParams& params) {
    filter_params_ = params;
}

const BamFilterParams& BamReader::get_filter_params() const {
    return filter_params_;
}

std::string BamReader::get_sequence_name(int id) const {
    if (!header_ || id < 0 || id >= header_->n_targets) {
        return "";
    }
    return header_->target_name[id];
}

int BamReader::get_sequence_id(const std::string& name) const {
    if (!header_) {
        return -1;
    }
    return bam_name2id(header_.get(), name.c_str());
}

bool BamReader::set_region(const std::string& chrom, int32_t start, int32_t end) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!header_ || !index_) {
        return false;
    }

    // Clean up previous iterator if it exists
    iterator_.reset();

    // Get chromosome ID from name
    int tid = bam_name2id(header_.get(), chrom.c_str());
    if (tid < 0) {
        return false;
    }

    // Create new iterator for the region
    iterator_.reset(sam_itr_queryi(index_.get(), tid, start - 1, end)); // 0-based coordinates
    if (!iterator_) {
        return false;
    }

    // Get the first read
    has_more_reads_ = next_read();
    return has_more_reads_;
}

bool BamReader::set_region(const std::string& region_str) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!header_ || !index_) {
        return false;
    }

    // Clean up previous iterator if it exists
    iterator_.reset();

    // Parse the region string using regex
    std::regex region_regex("([^:]+):(\\d+)-(\\d+)");
    std::smatch matches;
    if (!std::regex_match(region_str, matches, region_regex)) {
        return false;
    }

    std::string chrom = matches[1].str();
    int32_t start = std::stoi(matches[2].str());
    int32_t end = std::stoi(matches[3].str());

    // Call the other set_region without lock since we already have the lock
    int tid = bam_name2id(header_.get(), chrom.c_str());
    if (tid < 0) {
        return false;
    }

    // Create new iterator for the region
    iterator_.reset(sam_itr_queryi(index_.get(), tid, start - 1, end)); // 0-based coordinates
    if (!iterator_) {
        return false;
    }

    // Get the first read
    has_more_reads_ = next_read();
    return has_more_reads_;
}

bool BamReader::has_next_read() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return has_more_reads_;
}

bool BamReader::next_read() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!iterator_ || !read_) {
        has_more_reads_ = false;
        return false;
    }

    int ret;
    do {
        ret = sam_itr_next(file_.get(), iterator_.get(), read_.get());
        if (ret < 0) {
            has_more_reads_ = false;
            return false;
        }
    } while (!passes_filters(read_.get()));

    has_more_reads_ = true;
    return true;
}

bool BamReader::passes_filters(const bam1_t* b) const {
    // Always filter out unmapped reads
    if ((b->core.flag & BAM_FUNMAP) || b->core.tid < 0) {
        return false;
    }

    // Check mapping quality
    if (b->core.qual < filter_params_.min_mapping_quality) {
        return false;
    }

    // Check exclude flags
    if ((b->core.flag & filter_params_.exclude_flags) != 0) {
        return false;
    }

    // Check required flags if any
    if (filter_params_.required_flags != 0 &&
        (b->core.flag & filter_params_.required_flags) != filter_params_.required_flags) {
        return false;
    }

    // Check if we should count orphans - by default, samtools doesn't count orphans (-A to count them)
    if (!filter_params_.has_flag(ReadFilterFlag::COUNT_ORPHANS) &&
        (b->core.flag & BAM_FPAIRED) &&
        !(b->core.flag & BAM_FPROPER_PAIR)) {
        return false;
    }

    // Check if we should filter reads with large clips
    if (filter_params_.has_flag(ReadFilterFlag::FILTER_CLIPPED)) {
        uint32_t *cigar = bam_get_cigar(b);
        uint32_t qlen = b->core.l_qseq;
        uint32_t clip_len = 0;

        // Check first and last CIGAR operations for clips
        if (b->core.n_cigar > 0) {
            // First op
            uint32_t first_op = cigar[0] & BAM_CIGAR_MASK;
            if (first_op == BAM_CSOFT_CLIP || first_op == BAM_CHARD_CLIP) {
                clip_len += cigar[0] >> BAM_CIGAR_SHIFT;
            }

            // Last op
            uint32_t last_op = cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK;
            if (last_op == BAM_CSOFT_CLIP || last_op == BAM_CHARD_CLIP) {
                clip_len += cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT;
            }
        }

        // If more than 20% of read is clipped, filter it out
        if (qlen > 0 && (double)clip_len / qlen > 0.2) {
            return false;
        }
    }

    return true;
}

std::string BamReader::get_tag_value(const bam1_t* b, const char* tag_name) const {
    uint8_t* aux;
    if ((aux = bam_aux_get(b, tag_name)) == nullptr) {
        return "";
    }

    switch (aux[0]) {
        case 'Z':
            return bam_aux2Z(aux);
        case 'A':
            return std::string(1, bam_aux2A(aux));
        case 'c':
            return std::to_string(bam_aux2i(aux));
        case 'C':
            return std::to_string(bam_aux2i(aux));
        case 's':
            return std::to_string(bam_aux2i(aux));
        case 'S':
            return std::to_string(bam_aux2i(aux));
        case 'i':
            return std::to_string(bam_aux2i(aux));
        case 'I':
            return std::to_string(bam_aux2i(aux));
        case 'f':
            return std::to_string(bam_aux2f(aux));
        default:
            return "";
    }
}

char BamReader::base_as_char(uint8_t base) const {
    static const char bases[] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
    return bases[base & 0xf];
}

// Helper method to set up an iterator for a genomic position
bool BamReader::setup_position_iterator(const std::string& chrom, int32_t position) {
    // Requirements check
    if (!header_ || !index_) {
        return false;
    }

    // Clean up previous iterator if it exists
    iterator_.reset();

    // Get chromosome ID from name
    int tid = bam_name2id(header_.get(), chrom.c_str());
    if (tid < 0) {
        return false;
    }

    // Create new iterator for the region (position-1 because BAM uses 0-based coordinates)
    iterator_.reset(sam_itr_queryi(index_.get(), tid, position - 1, position));
    if (!iterator_) {
        return false;
    }

    // Need to initialize read_ before using it
    if (!read_) {
        return false;
    }

    // Get the first read
    int ret;
    do {
        ret = sam_itr_next(file_.get(), iterator_.get(), read_.get());
        if (ret < 0) {
            has_more_reads_ = false;
            return false;
        }
    } while (!passes_filters(read_.get()));

    has_more_reads_ = true;
    return true;
}

// Helper method to fetch the next read that passes filters
bool BamReader::fetch_next_filtered_read(int mapq_threshold) {
    if (!iterator_ || !read_) {
        has_more_reads_ = false;
        return false;
    }

    int ret;
    do {
        ret = sam_itr_next(file_.get(), iterator_.get(), read_.get());
        if (ret < 0) {
            has_more_reads_ = false;
            return false;
        }

        // Skip reads that don't pass the mapping quality threshold
        if (read_->core.qual < mapq_threshold) {
            continue;
        }
    } while (!passes_filters(read_.get()));

    has_more_reads_ = (ret >= 0);
    return has_more_reads_;
}

// Helper method to find the base in a read that corresponds to a reference position
std::optional<std::tuple<char, uint8_t>> BamReader::find_base_at_ref_position(
    const bam1_t* read,
    int pos0
) {
    // Skip if position is not covered by this read
    if (pos0 < read->core.pos || pos0 >= bam_endpos(read)) {
        return std::nullopt;
    }

    // Get base and quality at the position
    uint8_t* seq = bam_get_seq(read);
    uint8_t* qual = bam_get_qual(read);

    // Find the position in the read that aligns to our position of interest
    int read_pos = 0;  // Position within the read
    int ref_pos = read->core.pos;  // Position on the reference
    uint32_t* cigar = bam_get_cigar(read);

    // Traverse the CIGAR string
    for (uint32_t i = 0; i < read->core.n_cigar; ++i) {
        uint32_t op = cigar[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;

        switch(op) {
            case BAM_CMATCH:  // M: match or mismatch
            case BAM_CEQUAL:  // =: match
            case BAM_CDIFF:   // X: mismatch
                // If our position lies within this CIGAR operation
                if (pos0 >= ref_pos && pos0 < ref_pos + static_cast<int>(len)) {
                    int offset = pos0 - ref_pos;
                    read_pos += offset;

                    // Get base at this position
                    int base_idx = read_pos;
                    if (base_idx >= 0 && base_idx < read->core.l_qseq) {
                        uint8_t base = bam_seqi(seq, base_idx);
                        uint8_t base_qual = qual[base_idx];

                        // For Illumina 1.3+ encoding, convert quality
                        if (filter_params_.has_flag(ReadFilterFlag::ILLUMINA13_QUALITY)) {
                            base_qual = (base_qual <= 93) ? base_qual - 31 : 0;
                        }

                        // Return the base and quality
                        return std::make_tuple(base_as_char(base), base_qual);
                    }
                    return std::nullopt; // Position is in the read but base can't be retrieved
                }
                read_pos += len;
                ref_pos += len;
                break;

            case BAM_CINS:    // I: insertion
                read_pos += len;
                break;

            case BAM_CDEL:    // D: deletion
            case BAM_CREF_SKIP:  // N: skipped region
                ref_pos += len;
                break;

            case BAM_CSOFT_CLIP:  // S: soft clipping
                read_pos += len;
                break;

            case BAM_CHARD_CLIP:  // H: hard clipping
                // Hard clips are not in the read sequence
                break;

            case BAM_CPAD:    // P: padding
                break;

            default:
                break;
        }
    }

    // If we get here, the position wasn't found in any CIGAR operation
    return std::nullopt;
}

// Helper method to process a read at the specified position
bool BamReader::process_read_at_position(
    const bam1_t* read,
    int pos0,
    std::vector<const bam1_t*>& reads,
    std::vector<BamAlignment>& alignments,
    std::vector<char>& bases,
    std::vector<uint8_t>& qualities,
    std::unordered_map<const bam1_t*, size_t>& read_indices,
    int& read_count
) {
    // Skip if position is not covered by this read
    if (pos0 < read->core.pos || pos0 >= bam_endpos(read)) {
        return false;
    }

    // Create the alignment object for this read
    BamAlignment alignment;
    alignment.query_name = bam_get_qname(read);
    alignment.chrom = get_sequence_name(read->core.tid);
    alignment.position = read->core.pos + 1;  // Convert to 1-based
    alignment.mapping_quality = read->core.qual;
    alignment.flag = read->core.flag;
    alignment.is_reverse_strand = (read->core.flag & BAM_FREVERSE) != 0;
    alignment.cigar_string = get_cigar_string(read);
    alignment.seq = get_sequence_string(read);
    alignment.qual = get_quality_string(read);

    alignments.push_back(alignment);

    // Find the base at the position of interest
    auto base_info = find_base_at_ref_position(read, pos0);
    if (!base_info) {
        return false;
    }

    // Extract the base and quality
    auto [base_char, base_qual] = *base_info;

    // Store the base and quality for potential overlap detection
    bases.push_back(base_char);
    qualities.push_back(base_qual);
    read_indices[read] = reads.size();
    reads.push_back(read);
    read_count++;

    return true;
}

std::string BamReader::get_cigar_string(const bam1_t* b) const {
    std::ostringstream cigar_str;
    uint32_t* cigar = bam_get_cigar(b);

    for (uint32_t i = 0; i < b->core.n_cigar; ++i) {
        uint32_t op = cigar[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;

        cigar_str << len;

        switch(op) {
            case BAM_CMATCH:  cigar_str << 'M'; break;
            case BAM_CINS:    cigar_str << 'I'; break;
            case BAM_CDEL:    cigar_str << 'D'; break;
            case BAM_CREF_SKIP: cigar_str << 'N'; break;
            case BAM_CSOFT_CLIP: cigar_str << 'S'; break;
            case BAM_CHARD_CLIP: cigar_str << 'H'; break;
            case BAM_CPAD:    cigar_str << 'P'; break;
            case BAM_CEQUAL:  cigar_str << '='; break;
            case BAM_CDIFF:   cigar_str << 'X'; break;
            default:          cigar_str << '?'; break;
        }
    }

    return cigar_str.str();
}

std::string BamReader::get_sequence_string(const bam1_t* b) const {
    std::string seq;
    seq.reserve(b->core.l_qseq);

    uint8_t* bam_seq = bam_get_seq(b);
    for (int i = 0; i < b->core.l_qseq; ++i) {
        seq.push_back(base_as_char(bam_seqi(bam_seq, i)));
    }

    return seq;
}

std::string BamReader::get_quality_string(const bam1_t* b) const {
    std::string qual;
    qual.reserve(b->core.l_qseq);

    uint8_t* bam_qual = bam_get_qual(b);

    // Check if we need to convert from Illumina 1.3+ encoding
    if (filter_params_.has_flag(ReadFilterFlag::ILLUMINA13_QUALITY)) {
        for (int i = 0; i < b->core.l_qseq; ++i) {
            uint8_t q = bam_qual[i];
            q = (q <= 93) ? q - 31 : 0; // Convert Illumina 1.3+ to Sanger quality
            qual.push_back(static_cast<char>(33 + q)); // Convert to Phred+33 ASCII
        }
    } else {
        // Standard Sanger quality
        for (int i = 0; i < b->core.l_qseq; ++i) {
            qual.push_back(static_cast<char>(33 + bam_qual[i])); // Convert to Phred+33 ASCII
        }
    }

    return qual;
}

void BamReader::for_each_base_at_position(
    const std::string& chrom,
    int32_t position,
    std::function<void(char base_char, uint8_t base_quality, const BamAlignment* alignment)> processor,
    int min_base_quality,
    int min_mapping_quality
) {
    // Use filter params if thresholds not specified
    int base_qual_threshold = (min_base_quality >= 0) ? min_base_quality : filter_params_.min_base_quality;
    int mapq_threshold = (min_mapping_quality >= 0) ? min_mapping_quality : filter_params_.min_mapping_quality;

    // Temporary storage for potentially overlapping reads
    std::vector<const bam1_t*> reads;
    std::vector<BamAlignment> alignments;
    std::vector<char> bases;
    std::vector<uint8_t> qualities;
    std::unordered_map<const bam1_t*, size_t> read_indices;
    int read_count = 0;

    // We need to lock for the entire duration of this operation
    // to prevent iterator_ and read_ from being modified by other threads
    std::lock_guard<std::mutex> lock(mutex_);

    // Set up an iterator for this genomic position
    if (!setup_position_iterator(chrom, position)) {
        return;
    }

    // Convert position to 0-based for internal htslib operations
    int pos0 = position - 1;

    // Process each read overlapping the position
    while (has_more_reads_ && (filter_params_.max_depth <= 0 || read_count < filter_params_.max_depth)) {
        auto* b = read_.get();

        // Skip if mapping quality is too low
        if (b->core.qual < mapq_threshold) {
            if (!fetch_next_filtered_read(mapq_threshold)) {
                break;
            }
            continue;
        }

        // Process the current read at the specified position
        if (process_read_at_position(b, pos0, reads, alignments, bases, qualities, read_indices, read_count)) {
            // If we processed a read successfully, fetch the next one
            if (!fetch_next_filtered_read(mapq_threshold)) {
                break;
            }
        } else {
            // If the read doesn't cover the position, move to the next one
            if (!fetch_next_filtered_read(mapq_threshold)) {
                break;
            }
        }
    }

    // Handle overlapping pairs if needed
    if (!filter_params_.has_flag(ReadFilterFlag::DISABLE_OVERLAP)) {
        handle_overlapping_pairs(bases, qualities, read_indices, reads);
    }

    // Release the lock before calling the processor function as it might be slow
    // and we don't want to keep the lock for longer than necessary
    // We need to make a copy of the data that will be processed
    mutex_.unlock();

    // Process the bases with the processor function
    for (size_t i = 0; i < bases.size(); ++i) {
        if (qualities[i] >= base_qual_threshold) {
            processor(bases[i], qualities[i], &alignments[i]);
        }
    }

    // Reacquire the lock before returning (as required by unlock() in lock_guard's destructor)
    mutex_.lock();
}

void BamReader::handle_overlapping_pairs(
    std::vector<char>& bases,
    std::vector<uint8_t>& quals,
    std::unordered_map<const bam1_t*, size_t>& read_indices,
    const std::vector<const bam1_t*>& alignments
) const {
    // Map of read names to indices
    std::unordered_map<std::string, std::vector<size_t>> read_name_to_indices;

    // Group reads by name
    for (size_t i = 0; i < alignments.size(); ++i) {
        const char* qname = bam_get_qname(alignments[i]);
        read_name_to_indices[qname].push_back(i);
    }

    // For each read name with multiple reads (potential overlaps)
    for (const auto& [name, indices] : read_name_to_indices) {
        if (indices.size() < 2) continue;

        // Check each pair of reads with the same name
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = i + 1; j < indices.size(); ++j) {
                const size_t idx1 = indices[i];
                const size_t idx2 = indices[j];

                // Check if the two reads are paired
                const bam1_t* r1 = alignments[idx1];
                const bam1_t* r2 = alignments[idx2];

                if (!(r1->core.flag & BAM_FPAIRED) || !(r2->core.flag & BAM_FPAIRED)) {
                    continue;
                }

                // Check if these are proper pairs (opposite strands)
                bool r1_reverse = r1->core.flag & BAM_FREVERSE;
                bool r2_reverse = r2->core.flag & BAM_FREVERSE;

                if (r1_reverse == r2_reverse) {
                    continue; // Not opposite strands
                }

                // We have an overlapping pair
                char b1 = bases[idx1];
                char b2 = bases[idx2];
                uint8_t q1 = quals[idx1];
                uint8_t q2 = quals[idx2];

                // Determine which base to keep
                if (b1 == b2) {
                    // If they agree, keep the first base but sum the qualities
                    // Cap the quality at 200 to prevent overflow
                    quals[idx1] = std::min<uint16_t>(200, q1 + q2);
                    // Mark the second base for removal with quality 0
                    quals[idx2] = 0;
                } else {
                    // If they disagree, keep the higher quality base with 80% of its quality
                    if (q1 >= q2) {
                        quals[idx1] = 0.8 * q1;
                        quals[idx2] = 0;
                    } else {
                        quals[idx2] = 0.8 * q2;
                        quals[idx1] = 0;
                    }
                }
            }
        }
    }
}
