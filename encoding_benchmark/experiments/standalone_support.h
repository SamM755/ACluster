#pragma once

#include "utils/csv_reader.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

inline void verify_rows(
    std::vector<std::vector<double>> original,
    std::vector<std::vector<double>> decoded
) {
    std::sort(original.begin(), original.end());
    std::sort(decoded.begin(), decoded.end());
    if (original.size() != decoded.size()) throw std::runtime_error("Decoded row count mismatch");
    for (std::size_t row = 0; row < original.size(); ++row) {
        if (original[row].size() != decoded[row].size()) {
            throw std::runtime_error("Decoded column count mismatch");
        }
        for (std::size_t column = 0; column < original[row].size(); ++column) {
            const double left = original[row][column];
            const double right = decoded[row][column];
            if (std::isnan(left) != std::isnan(right)) throw std::runtime_error("Decoded NaN mismatch");
            if (!std::isnan(left) && std::abs(left - right) > 1e-9) {
                throw std::runtime_error("Decoded value mismatch");
            }
        }
    }
}

inline void write_string_rows(
    const std::vector<std::vector<std::string>>& rows,
    const std::string& path
) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream output(path);
    if (!output.is_open()) throw std::runtime_error("Cannot write temporary CSV: " + path);
    for (const auto& row : rows) {
        for (std::size_t column = 0; column < row.size(); ++column) {
            if (column != 0) output << ',';
            output << row[column];
        }
        output << '\n';
    }
}

inline double parse_number(const std::string& value) {
    try {
        const double parsed = std::stod(value);
        return std::isfinite(parsed) ? parsed : 0.0;
    } catch (...) {
        return 0.0;
    }
}
