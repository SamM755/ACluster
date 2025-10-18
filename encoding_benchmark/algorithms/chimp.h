#ifndef CHIMP_H
#define CHIMP_H

#include <vector>
#include <cstdint>
#include <string>

class Chimp {
public:
    std::vector<uint8_t> encode(const std::vector<double>& values);
    std::vector<double> decode(const std::vector<uint8_t>& compressed_data);
    std::string get_name() const { return "Chimp"; }
};

#endif // CHIMP_H