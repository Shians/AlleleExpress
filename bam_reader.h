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
    uint32_t flags = 0; // Bit flags for different filtering options

    // Methods to check/set flags
    bool has_flag(ReadFilterFlag flag) const {
        return (flags & std::to_underlying(flag)) != 0;
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

private:
    // Helper method to set up an iterator for a genomic position
    bool setup_position_iterator(const std::string& chrom, int32_t position);

    // Helper method to fetch the next read that passes filters
    bool fetch_next_filtered_read(int mapq_threshold);

    // Helper method to process a read at the specified position
    bool process_read_at_position(
        const bam1_t* read,
        int pos0,
        std::vector<const bam1_t*>& reads,
        std::vector<BamAlignment>& alignments,
        std::vector<char>& bases,
        std::vector<uint8_t>& qualities,
        std::unordered_map<const bam1_t*, size_t>& read_indices,
        int& read_count
    );

    // Helper method to find the base in a read that corresponds to a reference position
    std::optional<std::tuple<char, uint8_t>> find_base_at_ref_position(
        const bam1_t* read,
        int pos0
    );

    // Utility methods
    char base_as_char(uint8_t base) const;
    std::string get_cigar_string(const bam1_t* b) const;
    std::string get_sequence_string(const bam1_t* b) const;
    std::string get_quality_string(const bam1_t* b) const;
    std::string get_tag_value(const bam1_t* b, const char* tag_name) const;
    bool passes_filters(const bam1_t* b) const;

    // Handle overlapping paired reads
    void handle_overlapping_pairs(
        std::vector<char>& bases,
        std::vector<uint8_t>& quals,
        std::unordered_map<const bam1_t*, size_t>& read_indices,
        const std::vector<const bam1_t*>& alignments
    ) const;

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
