#pragma once

#include <string>
#include <memory>

#include "common.h"
#include "bam_reader.h"
#include "variant.h"

class AlleleCounter {
public:
    explicit AlleleCounter(BamReader* bam_reader);

    // Count REF and ALT alleles for a specific variant
    VariantResult count_alleles(const Variant& variant,
                                int min_base_quality = DEFAULT_MIN_BASE_QUALITY,
                                int min_mapping_quality = DEFAULT_MIN_MAPPING_QUALITY);

private:
    // Case-insensitive character comparison for allele matching
    static bool is_same_base(char base1, char base2);

    // Non-owning pointer to BAM reader
    BamReader* bam_reader_;
};
