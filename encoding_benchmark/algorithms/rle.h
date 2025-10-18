// in algorithms/rle.h
#ifndef RLE_H
#define RLE_H

#include <vector>
#include <cstdint>
#include <string>

class Rle {
public:
    
    std::vector<uint8_t> encode(const std::vector<double>& data);
    std::vector<double> decode(const std::vector<uint8_t>& compressed_data);
    
    std::string get_name() const { return "Rle"; }
};

#endif // RLE_H