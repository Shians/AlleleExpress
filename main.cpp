#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <future>

#include "argparse/argparse.hpp"
#include "bam_reader.h"
#include "bed_reader.h"
#include "allele_counter.h"
#include "result_writer.h"
#include "thread_pool.h"

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

    program.add_argument("-t", "--threads")
        .help("Number of threads to use (default: number of available CPU cores)")
        .scan<'i', int>()
        .default_value(DEFAULT_N_THREADS);

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
        int num_threads = program.get<int>("--threads");

        std::cerr << "Using " << num_threads << " threads for variant processing." << std::endl;

        // Initialize readers
        auto bam_reader = std::make_unique<BamReader>(bam_filename);
        auto bed_reader = std::make_unique<BedReader>(bed_filename, ref_filename);

        // Create a threadpool with the specified number of threads
        ThreadPool thread_pool(num_threads);

        // Read all positions into memory (could be optimized for very large BED files)
        std::vector<Variant> positions;
        while (auto position = bed_reader->next_position()) {
            positions.push_back(*position);
        }

        std::cerr << "Processing " << positions.size() << " positions using " << num_threads << " threads..." << std::endl;

        // Submit jobs to the thread pool and collect futures
        std::vector<std::future<VariantResult>> futures;
        futures.reserve(positions.size());

        // Create allele counters for each thread (thread-safe approach)
        std::vector<std::unique_ptr<AlleleCounter>> counters;
        for (int i = 0; i < num_threads; ++i) {
            counters.push_back(std::make_unique<AlleleCounter>(bam_reader.get()));
        }

        // Process positions in parallel
        for (size_t i = 0; i < positions.size(); ++i) {
            // Use counter from the thread's index mod number of threads
            int counter_index = i % num_threads;

            auto future = thread_pool.enqueue(
                [&counters, counter_index, &positions, i, min_base_quality, min_mapping_quality]() {
                    return counters[counter_index]->count_alleles(
                        positions[i],
                        min_base_quality,
                        min_mapping_quality
                    );
                }
            );
            futures.push_back(std::move(future));
        }

        // Collect results in original order
        std::vector<VariantResult> results;
        results.reserve(futures.size());

        for (auto& future : futures) {
            results.push_back(future.get());
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
