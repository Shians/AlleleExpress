#include "vcf_reader.h"
#include <stdexcept>
#include <vector>
#include <utility>
#include <string_view>

// Include the htslib headers
#include <htslib/vcf.h>
#include <htslib/hts.h>

VcfReader::VcfReader(const std::string& filename)
    : file_(nullptr, [](htsFile* f) { if (f) hts_close(f); })
    , header_(nullptr, [](bcf_hdr_t* h) { if (h) bcf_hdr_destroy(h); })
    , record_(bcf_init(), [](bcf1_t* r) { if (r) bcf_destroy(r); }) {

    // Open the VCF file
    file_.reset(hts_open(filename.c_str(), "r"));
    if (!file_) {
        throw std::runtime_error("Failed to open VCF file: " + filename);
    }

    // Check if the file is VCF/BCF
    if (file_->format.format != htsExactFormat::vcf && file_->format.format != htsExactFormat::bcf) {
        throw std::runtime_error("Input file must be a VCF or BCF file");
    }

    // Read the header
    header_.reset(bcf_hdr_read(file_.get()));
    if (!header_) {
        throw std::runtime_error("Failed to read header from VCF file");
    }
}

VcfReader::~VcfReader() = default;

VcfReader::VcfReader(VcfReader&& other) noexcept
    : file_(std::move(other.file_))
    , header_(std::move(other.header_))
    , record_(std::move(other.record_)) {
}

VcfReader& VcfReader::operator=(VcfReader&& other) noexcept {
    if (this != &other) {
        file_ = std::move(other.file_);
        header_ = std::move(other.header_);
        record_ = std::move(other.record_);
    }
    return *this;
}

std::string VcfReader::get_sample_name(int idx) const {
    if (!header_ || idx < 0 || idx >= bcf_hdr_nsamples(header_.get())) {
        return "";
    }
    return header_->samples[idx];
}

int VcfReader::get_num_samples() const {
    if (!header_) {
        return 0;
    }
    return bcf_hdr_nsamples(header_.get());
}

Variant VcfReader::convert_bcf_to_variant(bcf1_t* bcf_record) const {
    Variant variant;

    // Get chromosome name
    variant.chromosome = bcf_seqname(header_.get(), bcf_record);

    // Get 1-based position
    variant.position = bcf_record->pos + 1; // Convert from 0-based to 1-based

    // Get variant ID
    if (bcf_record->d.id && bcf_record->d.id[0] != '.' && bcf_record->d.id[0] != '\0') {
        variant.id = bcf_record->d.id;
    } else {
        variant.id = ".";
    }

    // Get quality
    variant.quality = bcf_record->qual;

    // Get ref allele
    variant.ref_allele = bcf_record->d.allele[0];

    // Get alt allele (assuming only one alt allele for simplicity)
    if (bcf_record->n_allele > 1) {
        variant.alt_allele = bcf_record->d.allele[1];
    } else {
        variant.alt_allele = ".";
    }

    // Get filter information
    int32_t* filters = nullptr;
    int n_filters = 0;
    if (bcf_get_info_int32(header_.get(), bcf_record, "FILTER", &filters, &n_filters) > 0) {
        std::string filter_str;
        for (int i = 0; i < n_filters; ++i) {
            if (i > 0) filter_str += ";";
            filter_str += bcf_hdr_int2id(header_.get(), BCF_DT_ID, filters[i]);
        }
        variant.filter = filter_str;
        free(filters);
    } else {
        variant.filter = "PASS";
    }

    // Get INFO field as a string (simplified)
    // In a real-world application, you might want to parse specific INFO fields
    variant.info = ".";

    return variant;
}

std::optional<Variant> VcfReader::next_variant() {
    if (!file_ || !header_ || !record_) {
        return std::nullopt;
    }

    // Read the next record
    int ret = bcf_read(file_.get(), header_.get(), record_.get());
    if (ret < 0) {
        return std::nullopt; // End of file or error
    }

    // Unpack the record
    bcf_unpack(record_.get(), BCF_UN_ALL);

    // Check if this is a SNP (we're only interested in SNPs)
    if (record_->n_allele != 2 ||
        strlen(record_->d.allele[0]) != 1 ||
        strlen(record_->d.allele[1]) != 1) {
        // Skip non-SNPs, try next variant recursively
        return next_variant();
    }

    return convert_bcf_to_variant(record_.get());
}
