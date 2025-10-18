// in algorithms/huffman.cpp

#include "huffman.h"
#include "../utils/bit_stream.h"
#include "../utils/type_converter.h"
#include <map>
#include <queue>
#include <vector>
#include <memory>
#include <algorithm>
#include <stdexcept>
#include <cmath> // For std::pow, std::round


namespace { // Anonymous namespace for implementation details

struct ScalingResult {
    std::vector<int64_t> scaled_data;
    double scaling_factor;
};

ScalingResult scale_doubles(const std::vector<double>& data) {
    if (data.empty()) return {{}, 1.0};
    int max_decimals = 0;
    for (double val : data) {
        if (std::isnan(val) || std::isinf(val)) continue;
        std::string s = std::to_string(val);
        auto dot_pos = s.find('.');
        if (dot_pos != std::string::npos) {
            auto last_digit = s.find_last_not_of('0');
            if(last_digit != std::string::npos && last_digit > dot_pos) {
                max_decimals = std::max(max_decimals, static_cast<int>(last_digit - dot_pos));
            }
        }
    }
    double factor = std::pow(10.0, max_decimals);
    std::vector<int64_t> scaled_data;
    scaled_data.reserve(data.size());
    for (double val : data) {
        scaled_data.push_back(static_cast<int64_t>(std::round(val * factor)));
    }
    return {scaled_data, factor};
}

// =================================================================================
//                     Part 1: Classic Huffman Tree Building
// =================================================================================

struct HuffmanNode {
    uint8_t data;
    size_t freq;
    std::unique_ptr<HuffmanNode> left, right;
    HuffmanNode(uint8_t d, size_t f) : data(d), freq(f), left(nullptr), right(nullptr) {}
};

struct CompareNodes {
    bool operator()(const std::unique_ptr<HuffmanNode>& l, const std::unique_ptr<HuffmanNode>& r) {
        return l->freq > r->freq;
    }
};

void getCodeLengths(const std::unique_ptr<HuffmanNode>& root, uint8_t length, std::map<uint8_t, uint8_t>& codeLengths) {
    if (!root) return;
    if (root->left == nullptr && root->right == nullptr) {
        codeLengths[root->data] = (length == 0) ? 1 : length; // Handle single-node tree
        return;
    }
    getCodeLengths(root->left, length + 1, codeLengths);
    getCodeLengths(root->right, length + 1, codeLengths);
}

// =================================================================================
//                  Part 2: Canonical Huffman Code Generation
// =================================================================================

// Helper struct to hold symbol, code length, and the generated canonical code
struct CanonicalCode {
    uint8_t symbol;
    uint8_t length;
    uint32_t code; // Using uint32_t is safe as code length won't exceed 32

    bool operator<(const CanonicalCode& other) const {
        if (length != other.length) return length < other.length;
        return symbol < other.symbol;
    }
};

// Generates a canonical code map from a map of code lengths
std::vector<CanonicalCode> generateCanonicalCodes(const std::map<uint8_t, uint8_t>& codeLengths) {
    std::vector<CanonicalCode> codes;
    codes.reserve(codeLengths.size());
    for (auto const& [symbol, length] : codeLengths) {
        codes.push_back({symbol, length, 0});
    }

    std::sort(codes.begin(), codes.end());

    uint32_t current_code = 0;
    if (!codes.empty()) {
        for (size_t i = 0; i < codes.size(); ++i) {
            codes[i].code = current_code;
            if (i + 1 < codes.size()) {
                current_code = (current_code + 1) << (codes[i+1].length - codes[i].length);
            }
        }
    }
    return codes;
}

// =================================================================================
//                  Part 3: Decoding Helper for Canonical Codes
// =================================================================================
// A specialized tree/table for decoding canonical codes efficiently.
class CanonicalDecoderTable {
public:
    explicit CanonicalDecoderTable(const std::vector<CanonicalCode>& canonicalCodes) {
        if (canonicalCodes.empty()) return;

        max_len = 0;
        for(const auto& c : canonicalCodes) {
            max_len = std::max(max_len, c.length);
        }
        
        // Create a direct lookup table for fast decoding
        // The table maps a code of a certain length to its symbol
        lookup.resize(max_len + 1);
        for(const auto& c : canonicalCodes) {
            lookup[c.length][c.code] = c.symbol;
        }
    }

    // Decodes one symbol from the bitstream using the new, simpler logic
    uint8_t decodeSymbol(InputBitStream& in) {
        uint32_t current_code = 0;
        for (uint8_t current_len = 1; current_len <= max_len; ++current_len) {
            // Read one more bit
            current_code = (current_code << 1) | in.read_bit();
            
            // Check if this code exists for this length
            if (current_len < lookup.size() && lookup[current_len].count(current_code)) {
                return lookup[current_len].at(current_code);
            }
        }
        // If the loop finishes, it means we've read max_len bits and still haven't found a valid code.
        throw std::runtime_error("Huffman decode error: invalid code in bitstream.");
    }

private:
    uint8_t max_len = 0;
    // A vector of maps: lookup[length][code] -> symbol
    std::vector<std::map<uint32_t, uint8_t>> lookup;
};


} // end anonymous namespace

// =================================================================================
//                          PUBLIC ENCODE / DECODE METHODS
// =================================================================================

std::vector<uint8_t> huffman_core_encode(const std::vector<uint8_t>& data) {    if (data.empty()) return {};

    std::map<uint8_t, size_t> freq;
    for (uint8_t byte : data) freq[byte]++;

    // Handle single symbol case (no compression possible, but must be decodable)
    if (freq.size() <= 1) {
        OutputBitStream out;
        out.write_bits(data.size(), 64);
        out.write_bits(freq.size(), 8); // Write 0 or 1
        if (!freq.empty()) {
            out.write_bits(freq.begin()->first, 8); // Write the symbol
        }
        return out.get_bytes();
    }
    
    // --- Step 1: Build classic Huffman tree ---
    std::priority_queue<std::unique_ptr<HuffmanNode>, std::vector<std::unique_ptr<HuffmanNode>>, CompareNodes> minHeap;
    for (auto const& [byte, f] : freq) minHeap.push(std::make_unique<HuffmanNode>(byte, f));
    while (minHeap.size() > 1) {
        // Correctly and safely extract nodes from the priority queue
        std::unique_ptr<HuffmanNode> left = std::move(const_cast<std::unique_ptr<HuffmanNode>&>(minHeap.top()));
        minHeap.pop();
        
        std::unique_ptr<HuffmanNode> right = std::move(const_cast<std::unique_ptr<HuffmanNode>&>(minHeap.top()));
        minHeap.pop();

        auto top = std::make_unique<HuffmanNode>('$', left->freq + right->freq);
        top->left = std::move(left);
        top->right = std::move(right);
        minHeap.push(std::move(top));
    }
    // --- Step 2: Get code lengths from the tree ---
    std::map<uint8_t, uint8_t> codeLengths;
    if (!minHeap.empty()) getCodeLengths(minHeap.top(), 0, codeLengths);
    
    // --- Step 3: Generate canonical codes from lengths ---
    std::vector<CanonicalCode> canonicalCodes = generateCanonicalCodes(codeLengths);
    std::map<uint8_t, CanonicalCode> codeMap;
    for(const auto& c : canonicalCodes) codeMap[c.symbol] = c;

    // --- Step 4: Write new, compact header ---
    OutputBitStream out;
    out.write_bits(data.size(), 64); // Original size
    out.write_bits((canonicalCodes.size() == 256) ? 0 : canonicalCodes.size(), 8); // Number of codes
    
    for (const auto& code_info : canonicalCodes) {
        out.write_bits(code_info.symbol, 8);
        out.write_bits(code_info.length, 5); // Use 5 bits for length (max 32)
    }

    // --- Step 5: Encode data using canonical codes ---
    for (uint8_t byte : data) {
        const auto& code_info = codeMap.at(byte);
        out.write_bits(code_info.code, code_info.length);
    }

    return out.get_bytes();
}

std::vector<uint8_t> huffman_core_decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    InputBitStream in(compressed_data);
    
    size_t original_size = 0;
    try { original_size = in.read_bits(64); } catch (const std::out_of_range&) { return {}; }
    if (original_size == 0) return {};

    size_t num_codes = in.read_bits(8);
    if (num_codes == 0 && original_size > 0) num_codes = 256;

    if (num_codes == 0) return {};
    if (num_codes == 1) {
        uint8_t symbol = in.read_bits(8);
        return std::vector<uint8_t>(original_size, symbol);
    }
    
    // --- Step 1: Read the compact code length header ---
    std::map<uint8_t, uint8_t> codeLengths;
    for (size_t i = 0; i < num_codes; ++i) {
        uint8_t symbol = in.read_bits(8);
        uint8_t length = in.read_bits(5);
        codeLengths[symbol] = length;
    }
    
    // --- Step 2: Rebuild the exact same canonical code table ---
    std::vector<CanonicalCode> canonicalCodes = generateCanonicalCodes(codeLengths);
    
    // --- Step 3: Create a fast decoder table ---
    CanonicalDecoderTable decoder_table(canonicalCodes);

    // --- Step 4: Decode the data ---
    std::vector<uint8_t> decoded_data;
    decoded_data.reserve(original_size);
    for (size_t i = 0; i < original_size; ++i) {
        decoded_data.push_back(decoder_table.decodeSymbol(in));
    }

    return decoded_data;
}

std::vector<uint8_t> Huffman::encode(const std::vector<double>& data) {
    if (data.empty()) return {};

    // --- Step 1: Pre-processing (Scaling + Zigzag) ---
    ScalingResult scaling_res = scale_doubles(data);
    
    // --- Step 2: Convert scaled integers to a byte stream ---
    // We will apply Huffman to the bytes of the zigzag-encoded integers.
    std::vector<uint8_t> byte_stream;
    byte_stream.reserve(scaling_res.scaled_data.size() * sizeof(int64_t));
    for (int64_t val : scaling_res.scaled_data) {
        uint64_t zigzag_val = (static_cast<uint64_t>(val) << 1) ^ (val >> 63);
        // Write bytes in a consistent order (little-endian)
        for (size_t i = 0; i < sizeof(uint64_t); ++i) {
            byte_stream.push_back((zigzag_val >> (i * 8)) & 0xFF);
        }
    }

    // --- Step 3: Apply core Huffman encoding to the byte stream ---
    std::vector<uint8_t> compressed_bytes = huffman_core_encode(byte_stream);

    // --- Step 4: Prepend our own header (scaling factor) ---
    OutputBitStream final_output;
    final_output.write_bits(double_to_long(scaling_res.scaling_factor), 64);
    final_output.write_bytes(compressed_bytes);

    return final_output.get_bytes();
}

std::vector<double> Huffman::decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.size() < sizeof(double)) return {};

    InputBitStream in(compressed_data);

    // --- Step 1: Read our custom header ---
    double scaling_factor = long_to_double(in.read_bits(64));
    
    // The rest of the data is the Huffman-compressed byte stream
    std::vector<uint8_t> huffman_payload = in.read_remaining_bytes(); // Assuming this exists

    // --- Step 2: Apply core Huffman decoding ---
    std::vector<uint8_t> decoded_bytes = huffman_core_decode(huffman_payload);

    // --- Step 3: Convert byte stream back to scaled integers ---
    if (decoded_bytes.size() % sizeof(int64_t) != 0) {
        throw std::runtime_error("Huffman decoded byte stream size is invalid.");
    }
    
    std::vector<double> result;
    result.reserve(decoded_bytes.size() / sizeof(int64_t));
    for (size_t i = 0; i < decoded_bytes.size(); i += sizeof(int64_t)) {
        uint64_t zigzag_val = 0;
        for (size_t j = 0; j < sizeof(int64_t); ++j) {
            zigzag_val |= static_cast<uint64_t>(decoded_bytes[i + j]) << (j * 8);
        }
        
        // --- Step 4: Zigzag decode and Unscale ---
        int64_t scaled_val = (zigzag_val >> 1) ^ (-(zigzag_val & 1));
        result.push_back(static_cast<double>(scaled_val) / scaling_factor);
    }
    
    return result;
}