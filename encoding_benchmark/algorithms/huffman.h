// in algorithms/huffman.h
#ifndef HUFFMAN_H
#define HUFFMAN_H

#include <vector>
#include <cstdint>
#include <string>

class Huffman {
public:
    
    std::vector<uint8_t> encode(const std::vector<double>& data);
    std::vector<double> decode(const std::vector<uint8_t>& compressed_data);

    std::string get_name() const { return "Huffman"; }
};

#endif // HUFFMAN_H