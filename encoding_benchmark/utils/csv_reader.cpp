#include "csv_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<std::vector<double>> read_csv_columns(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::vector<std::vector<double>> columns;
    std::string line;
    bool first_line = true;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        int current_col = 0;

        if (first_line) {
            while (std::getline(ss, cell, ',')) {
                columns.push_back({});
            }
            first_line = false;
            ss.clear();
            ss.str(line);
        }

        while (std::getline(ss, cell, ',')) {
            if (current_col < columns.size()) {
                try {
                    columns[current_col].push_back(std::stod(cell));
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Warning: Could not convert '" << cell << "' to double. Skipping." << std::endl;
                }
            }
            current_col++;
        }
    }

    return columns;
}

std::vector<std::vector<std::string>> read_csv_rows_as_strings(const std::string& filename) {
    std::vector<std::vector<std::string>> rows;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            row.push_back(value);
        }
        rows.push_back(row);
    }
    return rows;
}

std::vector<std::vector<double>> read_csv_rows(const std::string& file_path) {
    std::ifstream file_stream(file_path);
    if (!file_stream.is_open()) throw std::runtime_error("Cannot open file: " + file_path);
    std::vector<std::vector<double>> all_rows;
    std::string line;
    while (std::getline(file_stream, line)) {
        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) continue;
        std::stringstream line_stream(line);
        std::string cell;
        std::vector<double> current_row;
        while (std::getline(line_stream, cell, ',')) {
            try {
                current_row.push_back(std::stod(cell));
            } catch (const std::exception&) {
                // Handle non-numeric values if necessary, e.g., push back 0.0 or NaN
                current_row.push_back(0.0);
            }
        }
        if (!current_row.empty()) all_rows.push_back(current_row);
    }
    return all_rows;
}