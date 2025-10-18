#ifndef CSV_READER_H
#define CSV_READER_H

#include <string>
#include <vector>

std::vector<std::vector<double>> read_csv_columns(const std::string& filename);

std::vector<std::vector<std::string>> read_csv_rows_as_strings(const std::string& filename);

std::vector<std::vector<double>> read_csv_rows(const std::string& file_path);
#endif // CSV_READER_H