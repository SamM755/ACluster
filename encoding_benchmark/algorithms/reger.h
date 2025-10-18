#ifndef REGER_H
#define REGER_H

#include <vector>
#include <cstdint>
#include <string>

// The class now operates on a matrix of scaled integer values
class Reger {
public:
    // Input: A vector of rows, where each row is a vector of scaled int64_t values.
    // Also takes the scaling factors for metadata storage.
    std::vector<uint8_t> encode_multidim(
        const std::vector<std::vector<int64_t>>& rows,
        const std::vector<double>& scaling_factors
    );

    // Output: Decoded rows and their corresponding scaling factors.
    std::pair<std::vector<std::vector<int64_t>>, std::vector<double>>
    decode_multidim(const std::vector<uint8_t>& compressed_data);

    std::string get_name() const { return "Reger"; }
};

#endif // REGER_H