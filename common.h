#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <array>

// Default filter parameters (matching samtools mpileup defaults)
constexpr int DEFAULT_MIN_BASE_QUALITY = 13;    // -Q option in samtools mpileup
constexpr int DEFAULT_MIN_MAPPING_QUALITY = 0;  // -q option in samtools mpileup
constexpr int DEFAULT_MAX_DEPTH = 8000;         // -d option in samtools mpileup

// BAM flags
constexpr uint16_t BAM_FPAIRED = 0x1;
constexpr uint16_t BAM_FPROPER_PAIR = 0x2;
constexpr uint16_t BAM_FUNMAP = 0x4;
constexpr uint16_t BAM_FMUNMAP = 0x8;
constexpr uint16_t BAM_FREVERSE = 0x10;
constexpr uint16_t BAM_FMREVERSE = 0x20;
constexpr uint16_t BAM_FREAD1 = 0x40;
constexpr uint16_t BAM_FREAD2 = 0x80;
constexpr uint16_t BAM_FSECONDARY = 0x100;
constexpr uint16_t BAM_FQCFAIL = 0x200;
constexpr uint16_t BAM_FDUP = 0x400;
constexpr uint16_t BAM_FSUPPLEMENTARY = 0x800;

// Default exclude flags (matching samtools mpileup defaults: SECONDARY,QCFAIL,DUP)
constexpr uint16_t DEFAULT_EXCLUDE_FLAGS = BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;

// Read filtering flags enum for additional filtering options beyond flags
enum class ReadFilterFlag : uint32_t {
    DISABLE_BAQ = 1 << 0,              // -B option: disable BAQ (Base Alignment Quality)
    RECALCULATE_BAQ = 1 << 1,          // -E option: recalculate BAQ on the fly
    ILLUMINA13_QUALITY = 1 << 2,       // -6 option: assume quality is in Illumina-1.3+ format
    COUNT_ORPHANS = 1 << 3,            // -A option: count anomalous read pairs
    DISABLE_OVERLAP = 1 << 4,          // -x option: disable read-pair overlap detection
    IGNORE_RG = 1 << 5,                // -R option: ignore read groups specified with -r
    ADJUST_MQ = 1 << 6,                // -C option: adjust mapping quality
    FILTER_UNMAPPED = 1 << 7,          // Always exclude unmapped reads
    NO_INDELS = 1 << 8,                // -I option: don't include indels
    FILTER_CLIPPED = 1 << 9,           // Filter out reads with large clips
    MIN_DEPTH_FILTER = 1 << 10,        // Apply minimum depth filter
    MAX_DEPTH_FILTER = 1 << 11,        // Apply maximum depth filter
    REPORT_DELETIONS = 1 << 12,        // Report deletions as extended CIGAR ('*')
    REPORT_INSERTIONS = 1 << 13,       // Report insertions as extended CIGAR ('+[ACGT]+')
    REPORT_MQUAL = 1 << 14             // Report mapping quality
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
