#ifndef ELF_H
#define ELF_H

#include <vector>
#include <cstdint>
#include <string>

class Elf {
public:
    std::vector<uint8_t> encode(const std::vector<double>& values);
    std::vector<double> decode(const std::vector<uint8_t>& compressed_data);
    std::string get_name() const { return "Elf"; }
};

#endif // ELF_H