#pragma once

#include <string>
#include <vector>
#include <fstream>

#include "common.h"

class ResultWriter {
public:
    ResultWriter() = default;

    // Write the results to either a file or stdout
    void write_results(const std::vector<VariantResult>& results, const std::string& filename = "");

private:
    // Helper function to format a single result as a string
    std::string format_result(const VariantResult& result) const;

    // Write header line
    void write_header(std::ostream& out) const;
};
