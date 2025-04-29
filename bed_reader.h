#pragma once

#include <string>
#include <memory>
#include <functional>
#include <optional>
#include <fstream>

#include "common.h"
#include "variant.h"

// Forward declaration for htslib FASTA file handling
extern "C" {
    struct faidx_t;
}

/**
 * @class BedReader
 * @brief Reads BED files and provides positions for variant analysis with reference alleles from FASTA
 */
class BedReader {
public:
    explicit BedReader(const std::string& bed_filename, const std::string& ref_filename);
    ~BedReader();

    // Disallow copy operations
    BedReader(const BedReader&) = delete;
    BedReader& operator=(const BedReader&) = delete;

    // Move operations
    BedReader(BedReader&&) noexcept;
    BedReader& operator=(BedReader&&) noexcept;

    // Get the next position from the BED file
    std::optional<Variant> next_position();

private:
    std::ifstream bed_file_;
    std::string current_line_;
    faidx_t* fai_ = nullptr;  // FASTA index pointer

    bool open_bed_file(const std::string& filename);
    bool open_reference(const std::string& filename);
    bool parse_next_line();

    // Get reference base at a specific position
    std::string get_ref_base(const std::string& chrom, int32_t position);
};
