#include "allele_counter.h"
#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <map>

AlleleCounter::AlleleCounter(BamReader* bam_reader)
    : bam_reader_(bam_reader) {
    if (!bam_reader_) {
        throw std::invalid_argument("BamReader cannot be null");
    }
}

bool AlleleCounter::is_same_base(char base1, char base2) {
    return std::toupper(base1) == std::toupper(base2);
}

VariantResult AlleleCounter::count_alleles(const Variant& variant,
                                           int min_base_quality,
                                           int min_mapping_quality) {
    VariantResult result;
    result.variant = variant;

    // Check that we have a valid reference base
    if (variant.ref_allele.length() != 1) {
        return result; // Return empty result for non-SNPs
    }

    // Extract the reference base
    char ref_base = variant.ref_allele[0];

    // Map to count occurrences of each base (using uppercase for consistency)
    std::map<char, int> base_counts;

    // Process all reads at this position
    bam_reader_->for_each_base_at_position(
        variant.chromosome,
        variant.position,
        [&](char base_char, uint8_t base_quality) {
            // Ignore 'N' bases
            if (base_char == 'N') {
                return;
            }

            // Convert to uppercase for consistent counting
            char base = std::toupper(base_char);

            // Count this base
            base_counts[base]++;

            // Update reference counts immediately
            if (is_same_base(base_char, ref_base)) {
                ++result.ref_count;
            }
        },
        min_base_quality,
        min_mapping_quality
    );

    // Find the most common non-reference base
    char alt_base = '.';
    int max_count = 0;

    for (const auto& entry : base_counts) {
        // Skip the reference base
        if (is_same_base(entry.first, ref_base)) {
            continue;
        }

        // If this base is more common, or it's a tie and alphabetically earlier
        if (entry.second > max_count || (entry.second == max_count && entry.first < alt_base)) {
            alt_base = entry.first;
            max_count = entry.second;
        }
    }

    // Update the variant with the determined alt allele
    if (max_count > 0) {
        result.variant.alt_allele = alt_base;
        result.alt_count = base_counts[alt_base];
    } else {
        result.variant.alt_allele = ".";
        result.alt_count = 0;
    }

    // Calculate other counts (bases that are neither ref nor alt)
    for (const auto& entry : base_counts) {
        if (!is_same_base(entry.first, ref_base) && !is_same_base(entry.first, alt_base)) {
            result.other_count += entry.second;
        }
    }

    // Calculate the variant allele frequency
    result.calculate_vaf();

    return result;
}
