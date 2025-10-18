#ifndef GORILLA_H
#define GORILLA_H

#include <vector>
#include <cstdint>
#include <string>

class Gorilla {
public:
    
    std::vector<uint8_t> encode(const std::vector<double>& values);
    std::vector<double> decode(const std::vector<uint8_t>& compressed_data);
    std::string get_name() const { return "Gorilla"; }
};

#endif // GORILLA_H