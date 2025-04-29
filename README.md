# AlleleExpress

AlleleExpress is a fast C++ utility for calculating allele frequencies from BAM files at specified genomic positions. It efficiently processes sequencing data to identify variants and quantify their respective allele frequencies.

## Features

- Fast processing of BAM files at target positions specified in BED files
- Accurate allele frequency calculation with customizable quality thresholds
- Output in a simple tab-delimited format for easy downstream analysis

## Requirements

- C++23 compatible compiler (GCC 7+ or Clang 5+)
- CMake 3.14 or higher
- HTSlib 1.9 or higher

## Installation

### From Source

1. Clone the repository:

```bash
git clone https://github.com/username/allele-express.git
cd allele-express
```

2. Create a build directory and compile:

```bash
cmake -S . -B build
cmake --build build
```

3. Install (optional):

```bash
sudo make install
```

## Usage

Basic usage:

```bash
allele-express <bam_file> -b <bed_file> -r <reference.fa> [options]
```

### Required Arguments

- `<bam_file>`: Input BAM file (must be indexed)
- `-b, --bed-file`: Input BED file with target positions
- `-r, --reference`: Reference genome in FASTA format (must be indexed)

### Optional Arguments

- `-o, --output`: Output file path (stdout if not specified)
- `-q, --min-base-quality`: Minimum base quality to include in counts (default: 13)
- `-m, --min-mapping-quality`: Minimum read mapping quality (default: 0)
- `-h, --help`: Show help message
- `-v, --version`: Print version information and exit

## Output Format

The output is a tab-delimited file with the following columns:

```
chr  pos  ref  alt  ref_count  alt_count  other_count
```

Where:
- `chr`: Chromosome name
- `pos`: 1-based position
- `ref`: Reference allele
- `alt`: Alternative allele
- `ref_count`: Number of reads supporting the reference allele
- `alt_count`: Number of reads supporting the alternative allele
- `other_count`: Number of reads with bases other than ref or alt

## License

This project is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This project uses the [argparse](https://github.com/p-ranav/argparse) library for command-line argument parsing
- HTSlib for BAM file handling
