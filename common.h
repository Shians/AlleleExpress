#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <array>
#include <utility>

// Default filter parameters (matching samtools mpileup defaults)
constexpr int DEFAULT_MIN_BASE_QUALITY = 13;    // -Q option in samtools mpileup
constexpr int DEFAULT_MIN_MAPPING_QUALITY = 0;  // -q option in samtools mpileup
constexpr int DEFAULT_MAX_DEPTH = 8000;         // -d option in samtools mpileup

// BAM flags as an enum class for type safety
enum class BamFlag : uint16_t {
    PAIRED           = 0x0001,  // Paired-end / multiple segment sequencing
    PROPER_PAIR      = 0x0002,  // Each segment properly aligned
    UNMAPPED         = 0x0004,  // Segment unmapped
    MATE_UNMAPPED    = 0x0008,  // Next segment unmapped
    REVERSE_STRAND   = 0x0010,  // SEQ is reverse complemented
    MATE_REVERSE     = 0x0020,  // SEQ of next segment reversed
    READ1            = 0x0040,  // First segment in template
    READ2            = 0x0080,  // Last segment in template
    SECONDARY        = 0x0100,  // Not primary alignment
    QC_FAIL          = 0x0200,  // Quality control failure
    DUPLICATE        = 0x0400,  // PCR or optical duplicate
    SUPPLEMENTARY    = 0x0800   // Supplementary alignment
};

// Helper functions for BamFlag
constexpr uint16_t flag_bits(BamFlag flag) {
    return std::to_underlying(flag);
}

inline bool has_flag(uint16_t flags, BamFlag flag) {
    return (flags & flag_bits(flag)) != 0;
}

// Default exclude flags - now excluding secondary, qcfail, duplicates, and supplementary alignments
constexpr uint16_t DEFAULT_EXCLUDE_FLAGS =
    flag_bits(BamFlag::SECONDARY) |
    flag_bits(BamFlag::QC_FAIL) |
    flag_bits(BamFlag::DUPLICATE) |
    flag_bits(BamFlag::SUPPLEMENTARY);

// Read filtering flags enum for additional filtering options beyond flags
enum class ReadFilterFlag : uint16_t {
    COUNT_ORPHANS      = 1 << 1,       // -A option: count anomalous read pairs
    DISABLE_OVERLAP    = 1 << 2,       // -x option: disable read-pair overlap detection
    FILTER_CLIPPED     = 1 << 3,       // Filter out reads with large clips
    ILLUMINA13_QUALITY = 1 << 4        // -6 option: assume quality is in Illumina-1.3+ format
};

// Structure to represent a variant allele
struct AlleleCount {
    char base;
    int count;
    int total_quality;

    AlleleCount(char b, int c, int q) : base(b), count(c), total_quality(q) {}

    // Default compare by count (descending)
    bool operator<(const AlleleCount& other) const {
        return count > other.count;  // Note: reversed for descending sort
    }
};

// Structure to represent a base at a genomic position
struct PositionResult {
    std::string chrom;
    int32_t position;  // 1-based
    char ref_base;
    int depth;
    std::vector<AlleleCount> alleles;

    // Utility methods
    int get_base_count(char base) const {
        for (const auto& allele : alleles) {
            if (allele.base == base) {
                return allele.count;
            }
        }
        return 0;
    }

    double get_allele_frequency(char base) const {
        if (depth == 0) return 0.0;
        return static_cast<double>(get_base_count(base)) / depth;
    }

    double get_average_quality(char base) const {
        for (const auto& allele : alleles) {
            if (allele.base == base && allele.count > 0) {
                return static_cast<double>(allele.total_quality) / allele.count;
            }
        }
        return 0.0;
    }
};

// Constants for standard nucleotides
const std::array<char, 5> STD_BASES = {'A', 'C', 'G', 'T', 'N'};

#endif // COMMON_H
