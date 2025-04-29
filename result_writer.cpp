#include "result_writer.h"
#include <iostream>
#include <iomanip>
#include <sstream>

void ResultWriter::write_header(std::ostream& out) const {
    out << "chr\tpos\tref\talt\tref_count\talt_count\tother_count" << std::endl;
}

std::string ResultWriter::format_result(const VariantResult& result) const {
    std::stringstream ss;

    ss << result.variant.chromosome << "\t"
       << result.variant.position << "\t"
       << result.variant.ref_allele << "\t"
       << result.variant.alt_allele << "\t"
       << result.ref_count << "\t"
       << result.alt_count << "\t"
       << result.other_count;

    return ss.str();
}

void ResultWriter::write_results(const std::vector<VariantResult>& results, const std::string& filename) {
    if (filename.empty()) {
        // Write to stdout if no filename provided
        write_header(std::cout);
        for (const auto& result : results) {
            std::cout << format_result(result) << std::endl;
        }
    } else {
        // Write to the specified file
        std::ofstream out_file(filename);
        if (!out_file) {
            throw std::runtime_error("Failed to open output file: " + filename);
        }

        write_header(out_file);
        for (const auto& result : results) {
            out_file << format_result(result) << std::endl;
        }

        out_file.close();
    }
}
