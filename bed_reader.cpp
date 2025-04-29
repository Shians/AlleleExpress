#include "bed_reader.h"
#include <stdexcept>
#include <sstream>
#include <utility>

// Include htslib FASTA file handling
#include <htslib/faidx.h>

BedReader::BedReader(const std::string& bed_filename, const std::string& ref_filename) {
    if (!open_bed_file(bed_filename)) {
        throw std::runtime_error("Failed to open BED file: " + bed_filename);
    }

    if (!open_reference(ref_filename)) {
        throw std::runtime_error("Failed to open reference genome: " + ref_filename);
    }
}

BedReader::~BedReader() {
    if (fai_) {
        fai_destroy(fai_);
    }
}

BedReader::BedReader(BedReader&& other) noexcept
    : bed_file_(std::move(other.bed_file_)),
      current_line_(std::move(other.current_line_)),
      fai_(other.fai_) {
    other.fai_ = nullptr;
}

BedReader& BedReader::operator=(BedReader&& other) noexcept {
    if (this != &other) {
        bed_file_ = std::move(other.bed_file_);
        current_line_ = std::move(other.current_line_);

        if (fai_) {
            fai_destroy(fai_);
        }
        fai_ = other.fai_;
        other.fai_ = nullptr;
    }
    return *this;
}

bool BedReader::open_bed_file(const std::string& filename) {
    bed_file_.open(filename);
    return bed_file_.is_open();
}

bool BedReader::open_reference(const std::string& filename) {
    fai_ = fai_load(filename.c_str());
    return fai_ != nullptr;
}

bool BedReader::parse_next_line() {
    while (std::getline(bed_file_, current_line_)) {
        // Skip empty lines and comments
        if (current_line_.empty() || current_line_[0] == '#') {
            continue;
        }
        return true;
    }
    return false;
}

std::string BedReader::get_ref_base(const std::string& chrom, int32_t position) {
    // position is 1-based in our Variant struct
    // faidx_fetch_seq returns sequence starting from 0-based position
    int fetch_pos = position - 1;

    // Fetch a single base at the position
    int len = 0;
    char* base = faidx_fetch_seq(fai_, chrom.c_str(), fetch_pos, fetch_pos, &len);

    if (!base || len <= 0) {
        // If we can't get the reference base, return N
        return "N";
    }

    std::string result(base, len);
    free(base);

    return result;
}

std::optional<Variant> BedReader::next_position() {
    if (!parse_next_line()) {
        return std::nullopt;  // End of file or error
    }

    std::istringstream ss(current_line_);
    Variant variant;

    // BED format has chromosome, start, end (and possibly more fields)
    // We only need chromosome and start position
    std::string start_pos_str;

    if (!(ss >> variant.chromosome >> start_pos_str)) {
        throw std::runtime_error("Invalid BED file format: " + current_line_);
    }

    try {
        // BED positions are 0-based, convert to 1-based for compatibility
        variant.position = std::stoi(start_pos_str) + 1;
    } catch (const std::exception& e) {
        throw std::runtime_error("Invalid position in BED file: " + start_pos_str);
    }

    // Get the reference base from the FASTA file
    variant.ref_allele = get_ref_base(variant.chromosome, variant.position);

    // For BED files, we don't have actual variants, so use placeholders for non-reference fields
    variant.id = ".";
    variant.alt_allele = ".";  // No alternate allele
    variant.quality = 0.0;
    variant.filter = "PASS";
    variant.info = ".";

    return variant;
}
