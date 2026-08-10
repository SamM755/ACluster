#include "experiments.h"
void write_binary_file(const std::string& path, const std::vector<uint8_t>& bytes) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot create stream");
    out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    out.close();
    if (!out) throw std::runtime_error("Cannot write stream");
}
std::vector<DatasetEntry> load_selected_dataset_manifest(const std::string& path, const std::string&) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open dataset manifest");
    std::string line;
    std::getline(in, line);
    std::vector<DatasetEntry> result;
    while (std::getline(in, line)) {
        std::stringstream stream(line);
        std::vector<std::string> fields;
        std::string field;
        while (std::getline(stream, field, ',')) fields.push_back(field);
        if (fields.size() != 5) throw std::runtime_error("Invalid dataset manifest");
        result.push_back({fields[0], std::stoi(fields[3]), std::stoi(fields[1])});
    }
    return result;
}
