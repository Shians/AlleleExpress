#include <iostream>
#include <string>
#include <memory>
#include <vector>

#include "argparse/argparse.hpp"
#include "bam_reader.h"
#include "bed_reader.h"
#include "allele_counter.h"
#include "result_writer.h"

int main(int argc, char* argv[]) {
    // Set up argument parser
    argparse::ArgumentParser program("allele-express", "1.0.0");

    program.add_description("AlleleExpress: A C++ utility for calculating allele frequencies from BAM files at specified genomic positions");

    program.add_argument("bam_file")
        .help("Input BAM file (must be indexed)");

    program.add_argument("-b", "--bed-file")
        .help("Input BED file with target positions")
        .required();

    program.add_argument("-r", "--reference")
        .help("Reference genome in FASTA format (must be indexed)")
        .required();

    program.add_argument("-o", "--output")
        .help("Output file path (stdout if not specified)")
        .default_value(std::string(""));

    program.add_argument("-q", "--min-base-quality")
        .help("Minimum base quality to include in counts")
        .scan<'i', int>()
        .default_value(DEFAULT_MIN_BASE_QUALITY);

    program.add_argument("-m", "--min-mapping-quality")
        .help("Minimum read mapping quality to include in counts")
        .scan<'i', int>()
        .default_value(DEFAULT_MIN_MAPPING_QUALITY);

    program.add_argument("-v", "--version")
        .help("Print version information and exit")
        .default_value(false)
        .implicit_value(true);

    try {
        program.parse_args(argc, argv);

        if (program.get<bool>("--version")) {
            std::cout << "AlleleExpress version 1.0.0" << std::endl;
            return 0;
        }

        // Get parsed arguments
        std::string bam_filename = program.get("bam_file");
        std::string bed_filename = program.get("--bed-file");
        std::string ref_filename = program.get("--reference");
        std::string output_filename = program.get("--output");
        int min_base_quality = program.get<int>("--min-base-quality");
        int min_mapping_quality = program.get<int>("--min-mapping-quality");

        // Initialize readers
        auto bam_reader = std::make_unique<BamReader>(bam_filename);

        auto bed_reader = std::make_unique<BedReader>(bed_filename, ref_filename);

        // Initialize counter
        auto allele_counter = std::make_unique<AlleleCounter>(bam_reader.get());

        // Process BED positions
        std::vector<VariantResult> results;
        while (auto position = bed_reader->next_position()) {
            auto result = allele_counter->count_alleles(
                *position,
                min_base_quality,
                min_mapping_quality
            );
            results.push_back(result);
        }

        // Write results
        auto writer = std::make_unique<ResultWriter>();
        writer->write_results(results, output_filename);

        std::cerr << "Successfully processed " << results.size() << " positions." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    return 0;
}
