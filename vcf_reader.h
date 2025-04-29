#pragma once

#include <string>
#include <memory>
#include <functional>
#include <optional>

#include "common.h"

// Forward declarations for htslib types
extern "C" {
    struct htsFile;
    struct bcf_hdr_t;
    struct bcf1_t;
}

class VcfReader {
public:
    explicit VcfReader(const std::string& filename);
    ~VcfReader();

    // Disallow copy operations
    VcfReader(const VcfReader&) = delete;
    VcfReader& operator=(const VcfReader&) = delete;

    // Move operations
    VcfReader(VcfReader&&) noexcept;
    VcfReader& operator=(VcfReader&&) noexcept;

    // Get the next variant from the VCF file
    std::optional<Variant> next_variant();

    // Get header information
    std::string get_sample_name(int idx) const;
    int get_num_samples() const;

private:
    // Utility function to convert htslib's BCF record to our Variant structure
    Variant convert_bcf_to_variant(bcf1_t* bcf_record) const;

    // Private implementation details
    std::unique_ptr<htsFile, std::function<void(htsFile*)>> file_;
    std::unique_ptr<bcf_hdr_t, std::function<void(bcf_hdr_t*)>> header_;
    std::unique_ptr<bcf1_t, std::function<void(bcf1_t*)>> record_;
};
