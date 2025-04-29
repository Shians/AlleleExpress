#include "variant.h"
#include <sstream>
#include <iomanip>

// Constructor with chromosome and position
Variant::Variant(const std::string& chrom, int32_t pos)
    : chromosome(chrom), position(pos) {
}

// Constructor with all fields
Variant::Variant(const std::string& chrom, int32_t pos,
                 const std::string& ref, const std::string& alt,
                 const std::string& id, double qual,
                 const std::string& filt, const std::string& info)
    : chromosome(chrom), position(pos), id(id), ref_allele(ref),
      alt_allele(alt), quality(qual), filter(filt), info(info) {
}

// Return variant as a string in VCF-like format
std::string Variant::to_string() const {
    std::stringstream ss;
    ss << chromosome << '\t'
       << position << '\t'
       << id << '\t'
       << ref_allele << '\t'
       << alt_allele << '\t';

    // Format quality with precision
    if (quality > 0) {
        ss << std::fixed << std::setprecision(2) << quality << '\t';
    } else {
        ss << "." << '\t';
    }

    ss << filter << '\t'
       << info;

    return ss.str();
}

// Constructor from a Variant
VariantResult::VariantResult(const Variant& v)
    : variant(v) {
}

// Calculate and return variant allele frequency
double VariantResult::vaf() const {
    int depth = total_depth();
    if (depth == 0) {
        return 0.0;
    }
    return static_cast<double>(alt_count) / depth;
}

// Calculate VAF and store in the object
void VariantResult::calculate_vaf() {
    variant_allele_frequency = vaf();
}

// Return result as a string
std::string VariantResult::to_string() const {
    std::stringstream ss;
    ss << variant.to_string() << '\t'
       << ref_count << '\t'
       << alt_count << '\t'
       << other_count << '\t'
       << total_depth() << '\t'
       << std::fixed << std::setprecision(4) << vaf();

    return ss.str();
}
