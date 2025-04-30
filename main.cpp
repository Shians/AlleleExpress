#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <future>
#include <map>

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
        std::vector<std::future<std::vector<VariantResult>>> futures;

        // Create allele counters for each thread (thread-safe approach)
        std::vector<std::unique_ptr<AlleleCounter>> counters;
        for (int i = 0; i < num_threads; ++i) {
            counters.push_back(std::make_unique<AlleleCounter>(bam_reader.get()));
        }

        // Group positions by chromosome for batch processing
        std::map<std::string, std::vector<Variant>> chromosome_positions;
        for (const auto& pos : positions) {
            chromosome_positions[pos.chromosome].push_back(pos);
        }

        std::cerr << "Grouped positions across " << chromosome_positions.size() << " chromosomes." << std::endl;

        // Counter for tracking assigned batches
        size_t batch_counter = 0;

        // Process positions in chromosome-based batches
        for (auto& [chromosome, chrom_positions] : chromosome_positions) {
            // Process up to 50 positions per batch
            for (size_t batch_start = 0; batch_start < chrom_positions.size(); batch_start += 50) {
                // Calculate the end position of the current batch
                size_t batch_end = std::min(batch_start + 50, chrom_positions.size());

                // Select counter index based on batch counter
                int counter_index = batch_counter % num_threads;
                batch_counter++;

                // Create a batch of positions for processing
                std::vector<Variant> batch(chrom_positions.begin() + batch_start, chrom_positions.begin() + batch_end);

                auto future = thread_pool.enqueue(
                    [&counters, counter_index, batch, min_base_quality, min_mapping_quality]() {
                        std::vector<VariantResult> batch_results;
                        batch_results.reserve(batch.size());

                        for (const auto& pos : batch) {
                            batch_results.push_back(counters[counter_index]->count_alleles(
                                pos,
                                min_base_quality,
                                min_mapping_quality
                            ));
                        }

                        return batch_results;
                    }
                );

                futures.push_back(std::move(future));
            }
        }

        // Collect and flatten results
        std::vector<VariantResult> results;
        for (auto& future : futures) {
            auto batch_results = future.get();
            results.insert(results.end(), batch_results.begin(), batch_results.end());
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
