#ifndef REGER_H
#define REGER_H

#include <vector>
#include <cstdint>
#include <string>

class Reger {
public:

    std::vector<uint8_t> encode_multidim(
        const std::vector<std::vector<int64_t>>& rows,
        const std::vector<double>& scaling_factors,
        std::vector<size_t>* output_row_order = nullptr
    );

    std::pair<std::vector<std::vector<int64_t>>, std::vector<double>>
    decode_multidim(const std::vector<uint8_t>& compressed_data);

    std::string get_name() const { return "Reger"; }
};

#endif
