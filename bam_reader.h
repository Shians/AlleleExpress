#ifndef BAM_READER_H
#define BAM_READER_H

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <optional>
#include <unordered_map>
#include <mutex>
#include "common.h"

// Forward declarations of htslib types
struct htsFile;
struct sam_hdr_t;
struct hts_idx_t;
struct hts_itr_t;
struct bam1_t;

// Struct to represent alignment information
struct BamAlignment {
    std::string query_name;
    std::string chrom;
    int32_t position;  // 1-based position
    uint8_t mapping_quality;
    uint16_t flag;
    bool is_reverse_strand;
    std::string cigar_string;
    std::string seq;
    std::string qual;
};

// Struct to store pileup data at a position
struct PileupData {
    std::string chrom;
    int32_t position;  // 1-based position
    std::vector<char> bases;
    std::vector<uint8_t> qualities;
    std::vector<BamAlignment> alignments;
};

// Structure for BAM read filtering parameters
struct BamFilterParams {
    int min_base_quality = DEFAULT_MIN_BASE_QUALITY;
    int min_mapping_quality = DEFAULT_MIN_MAPPING_QUALITY;
    int max_depth = DEFAULT_MAX_DEPTH;
    uint16_t exclude_flags = DEFAULT_EXCLUDE_FLAGS;
    uint16_t required_flags = 0; // No flags required by default
    int adjust_mq_value = 0; // Default: no mapping quality adjustment
    std::unordered_map<std::string, bool> excluded_read_groups;
    uint32_t flags = 0; // Bit flags for different filtering options

    // Methods to check/set flags
    bool has_flag(ReadFilterFlag flag) const {
        return (flags & static_cast<uint32_t>(flag)) != 0;
    }

    void set_flag(ReadFilterFlag flag, bool value = true) {
        if (value) {
            flags |= static_cast<uint32_t>(flag);
        } else {
            flags &= ~static_cast<uint32_t>(flag);
        }
    }
};

/**
 * @class BamReader
 * @brief Wrapper for htslib's BAM file reading capabilities
 */
class BamReader {
public:
    /**
     * @brief Construct a BAM reader
     * @param filename Path to the BAM file
     * @throws std::runtime_error if the file cannot be opened or is not a BAM file
     */
    BamReader(const std::string& filename);

    /**
     * @brief Destructor
     */
    ~BamReader();

    // Disable copy operations
    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;

    // Move operations
    BamReader(BamReader&& other) noexcept;
    BamReader& operator=(BamReader&& other) noexcept;

    /**
     * @brief Set filter parameters
     * @param params The filter parameters to use
     */
    void set_filter_params(const BamFilterParams& params);

    /**
     * @brief Get current filter parameters
     * @return The currently set filter parameters
     */
    const BamFilterParams& get_filter_params() const;

    /**
     * @brief Set minimum base quality threshold
     * @param quality Minimum base quality (Phred scale)
     */
    void set_min_base_quality(int quality) { filter_params_.min_base_quality = quality; }

    /**
     * @brief Set minimum mapping quality threshold
     * @param quality Minimum mapping quality
     */
    void set_min_mapping_quality(int quality) { filter_params_.min_mapping_quality = quality; }

    /**
     * @brief Set maximum read depth
     * @param depth Maximum number of reads to consider at any position
     */
    void set_max_depth(int depth) { filter_params_.max_depth = depth; }

    /**
     * @brief Set flag to count orphan reads (unpaired reads or improper pairs)
     * @param count_orphans true to count orphans, false to exclude them
     */
    void set_count_orphans(bool count_orphans) {
        filter_params_.set_flag(ReadFilterFlag::COUNT_ORPHANS, count_orphans);
    }

    /**
     * @brief Set flag to disable overlapping read detection
     * @param disable_overlap true to disable overlap handling, false to handle overlaps
     */
    void set_disable_overlap(bool disable_overlap) {
        filter_params_.set_flag(ReadFilterFlag::DISABLE_OVERLAP, disable_overlap);
    }

    /**
     * @brief Set flag for Illumina 1.3+ quality encoding
     * @param use_illumina13 true if qualities are in Illumina 1.3+ format
     */
    void set_illumina13_quality(bool use_illumina13) {
        filter_params_.set_flag(ReadFilterFlag::ILLUMINA13_QUALITY, use_illumina13);
    }

    /**
     * @brief Set flag to disable BAQ (Base Alignment Quality) computation
     * @param disable_baq true to disable BAQ, false to use it
     */
    void set_disable_baq(bool disable_baq) {
        filter_params_.set_flag(ReadFilterFlag::DISABLE_BAQ, disable_baq);
    }

    /**
     * @brief Set flag to recalculate BAQ on the fly (used with -E in samtools)
     * @param recalculate true to recalculate BAQ, false otherwise
     */
    void set_recalculate_baq(bool recalculate) {
        filter_params_.set_flag(ReadFilterFlag::RECALCULATE_BAQ, recalculate);
    }

    /**
     * @brief Set coefficient for adjusting mapping quality (samtools -C option)
     * @param coefficient Value for adjusting mapping quality (0 disables adjustment)
     */
    void set_adjust_mq_coefficient(int coefficient) {
        filter_params_.adjust_mq_value = coefficient;
        filter_params_.set_flag(ReadFilterFlag::ADJUST_MQ, coefficient > 0);
    }

    /**
     * @brief Set flag to filter out reads with large clips (>20% of read length)
     * @param filter_clipped true to filter clipped reads, false otherwise
     */
    void set_filter_clipped(bool filter_clipped) {
        filter_params_.set_flag(ReadFilterFlag::FILTER_CLIPPED, filter_clipped);
    }

    /**
     * @brief Set flag to ignore indels
     * @param no_indels true to ignore indels, false to include them
     */
    void set_no_indels(bool no_indels) {
        filter_params_.set_flag(ReadFilterFlag::NO_INDELS, no_indels);
    }

    /**
     * @brief Set flags to exclude from the analysis
     * @param flags Bit mask of flags that should cause a read to be excluded
     */
    void set_exclude_flags(uint16_t flags) { filter_params_.exclude_flags = flags; }

    /**
     * @brief Set flags that must be present for a read to be included
     * @param flags Bit mask of flags that must be present
     */
    void set_required_flags(uint16_t flags) { filter_params_.required_flags = flags; }

    /**
     * @brief Exclude specified read groups
     * @param read_group_ids List of read group IDs to exclude
     */
    void exclude_read_groups(const std::vector<std::string>& read_group_ids);

    /**
     * @brief Load read group exclusions from a file
     * @param filename File containing read group IDs to exclude (one per line)
     * @return true if successful, false otherwise
     */
    bool load_read_group_exclusions(const std::string& filename);

    /**
     * @brief Get the name of a sequence/chromosome by ID
     * @param id Sequence ID
     * @return Name of the sequence
     */
    std::string get_sequence_name(int id) const;

    /**
     * @brief Get the ID of a sequence/chromosome by name
     * @param name Name of the sequence
     * @return ID of the sequence, or -1 if not found
     */
    int get_sequence_id(const std::string& name) const;

    /**
     * @brief Get the length of a sequence/chromosome by ID
     * @param id Sequence ID
     * @return Length of the sequence, or -1 if not found
     */
    int32_t get_sequence_length(int id) const;

    /**
     * @brief Get all sequence/chromosome names in the BAM file
     * @return Vector of sequence names
     */
    std::vector<std::string> get_sequence_names() const;

    /**
     * @brief Set the region to query
     * @param chrom Chromosome/contig name
     * @param start Start position (1-based)
     * @param end End position (1-based, inclusive)
     * @return true if the region was set successfully, false otherwise
     */
    bool set_region(const std::string& chrom, int32_t start, int32_t end);

    /**
     * @brief Set the region to query using a string
     * @param region_str Region string in the format "chrom:start-end" (1-based)
     * @return true if the region was set successfully, false otherwise
     */
    bool set_region(const std::string& region_str);

    /**
     * @brief Check if there are more reads to retrieve
     * @return true if there are more reads, false otherwise
     */
    bool has_next_read() const;

    /**
     * @brief Move to the next read
     * @return true if a read was found, false if no more reads
     */
    bool next_read();

    /**
     * @brief Get the current alignment
     * @return Optional containing the alignment, or empty if no current alignment
     */
    std::optional<BamAlignment> get_current_alignment() const;

    /**
     * @brief Get all reads that overlap a position
     * @param chrom Chromosome/contig name
     * @param position Position (1-based)
     * @param min_mapping_quality Minimum mapping quality to include (-1 to use filter_params)
     * @return Vector of alignments that overlap the position
     */
    std::vector<BamAlignment> get_overlapping_reads(const std::string& chrom, int32_t position, int min_mapping_quality = -1);

    /**
     * @brief Process each base at a specific position
     * @param chrom Chromosome/contig name
     * @param position Position (1-based)
     * @param processor Function to process each base
     * @param min_base_quality Minimum base quality to consider (-1 to use filter_params)
     * @param min_mapping_quality Minimum mapping quality to consider (-1 to use filter_params)
     */
    void for_each_base_at_position(
        const std::string& chrom,
        int32_t position,
        std::function<void(char base_char, uint8_t base_quality, const BamAlignment* alignment)> processor,
        int min_base_quality = -1,
        int min_mapping_quality = -1
    );

    /**
     * @brief Get pileup data at a specific position
     * @param chrom Chromosome/contig name
     * @param position Position (1-based)
     * @param min_base_quality Minimum base quality to consider
     * @param min_mapping_quality Minimum mapping quality to consider
     * @return PileupData structure with bases and qualities at the position
     */
    PileupData get_pileup_at_position(
        const std::string& chrom,
        int32_t position,
        int min_base_quality = DEFAULT_MIN_BASE_QUALITY,
        int min_mapping_quality = DEFAULT_MIN_MAPPING_QUALITY
    );

    /**
     * @brief Check if the BAM file has an index
     * @return true if the BAM file has an index, false otherwise
     */
    bool has_index() const;

    /**
     * @brief Get the filename of the BAM file
     * @return Filename of the BAM file
     */
    std::string get_filename() const;

    /**
     * @brief Configure samtools mpileup-compatible default settings
     */
    void set_mpileup_defaults() {
        filter_params_.min_base_quality = DEFAULT_MIN_BASE_QUALITY;
        filter_params_.min_mapping_quality = DEFAULT_MIN_MAPPING_QUALITY;
        filter_params_.max_depth = DEFAULT_MAX_DEPTH;
        filter_params_.exclude_flags = DEFAULT_EXCLUDE_FLAGS;
        filter_params_.required_flags = 0;
        filter_params_.adjust_mq_value = 0;
        filter_params_.set_flag(ReadFilterFlag::COUNT_ORPHANS, false);
        filter_params_.set_flag(ReadFilterFlag::DISABLE_OVERLAP, false);
        filter_params_.set_flag(ReadFilterFlag::DISABLE_BAQ, true);
        filter_params_.set_flag(ReadFilterFlag::RECALCULATE_BAQ, false);
        filter_params_.set_flag(ReadFilterFlag::ILLUMINA13_QUALITY, false);
    }

private:
    // Utility methods
    char base_as_char(uint8_t base) const;
    std::string get_cigar_string(const bam1_t* b) const;
    std::string get_sequence_string(const bam1_t* b) const;
    std::string get_quality_string(const bam1_t* b) const;
    std::string get_tag_value(const bam1_t* b, const char* tag_name) const;
    bool passes_filters(const bam1_t* b) const;
    uint8_t adjust_mapping_quality(const bam1_t* b) const;

    // Handle overlapping paired reads
    void handle_overlapping_pairs(
        std::vector<char>& bases,
        std::vector<uint8_t>& quals,
        std::unordered_map<const bam1_t*, size_t>& read_indices,
        const std::vector<const bam1_t*>& alignments
    ) const;

    // Member variables
    std::string filename_;
    std::unique_ptr<htsFile, std::function<void(htsFile*)>> file_;
    std::unique_ptr<sam_hdr_t, std::function<void(sam_hdr_t*)>> header_;
    std::unique_ptr<hts_idx_t, std::function<void(hts_idx_t*)>> index_;
    std::unique_ptr<hts_itr_t, std::function<void(hts_itr_t*)>> iterator_;
    std::unique_ptr<bam1_t, std::function<void(bam1_t*)>> read_;
    bool has_more_reads_ = false;
    BamFilterParams filter_params_;
    mutable std::mutex mutex_; // Mutex for thread safety
};

#endif // BAM_READER_H
