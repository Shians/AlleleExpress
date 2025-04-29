#pragma once

#include <string>
#include <cstdint>

/**
 * @class Variant
 * @brief Represents a genomic variant (SNP, indel, etc.)
 */
class Variant {
public:
    // Default constructor
    Variant() = default;

    // Constructor with chromosome and position
    Variant(const std::string& chrom, int32_t pos);

    // Constructor with all fields
    Variant(const std::string& chrom, int32_t pos,
            const std::string& ref, const std::string& alt,
            const std::string& id = ".", double qual = 0.0,
            const std::string& filt = "PASS", const std::string& info = ".");

    // Return variant as a string in VCF-like format
    std::string to_string() const;

    // Public member variables
    std::string chromosome;
    int32_t position = 0;       // 1-based position
    std::string id = ".";       // Variant ID or "." if unknown
    std::string ref_allele;     // Reference allele
    std::string alt_allele;     // Alternative allele
    double quality = 0.0;       // Quality score
    std::string filter = "PASS"; // Filter status
    std::string info = ".";     // Additional information
};

/**
 * @class VariantResult
 * @brief Extends Variant with allele count results
 */
class VariantResult {
public:
    // Default constructor
    VariantResult() = default;

    // Constructor from a Variant
    explicit VariantResult(const Variant& v);

    // Calculate and return variant allele frequency
    double vaf() const;

    // Calculate VAF and store in the object
    void calculate_vaf();

    // Get total read depth
    int total_depth() const { return ref_count + alt_count + other_count; }

    // Return result as a string
    std::string to_string() const;

    // Public member variables
    Variant variant;                // The variant
    int ref_count = 0;             // Count of reference allele
    int alt_count = 0;             // Count of alternate allele
    int other_count = 0;           // Count of other alleles
    double variant_allele_frequency = 0.0;  // Calculated VAF
};
