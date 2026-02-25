#ifndef ALP_H
#define ALP_H

#include <vector>
#include <string>
#include <cstdint>

class Alp {
public:
    std::string get_name() const { return "ALP"; }

    std::vector<uint8_t> encode(const std::vector<double>& data);

    std::vector<double> decode(const std::vector<uint8_t>& data);
};

#endif // ALP_H