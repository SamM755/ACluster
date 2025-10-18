// in algorithms/chimp.cpp

#include "chimp.h"
#include "../utils/bit_stream.h"
#include "../utils/debug.h"
#include "../utils/type_converter.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <vector>
#include <stdexcept>

// Platform-specific functions for counting zeros
#ifdef _MSC_VER
#include <intrin.h>
#define count_leading_zeros(x) ( (x) == 0 ? 64 : __lzcnt64(x) )
#define count_trailing_zeros(x) ( (x) == 0 ? 64 : _tzcnt_u64(x) )
#else
#define count_leading_zeros(x) ( (x) == 0 ? 64 : __builtin_clzll(x) )
#define count_trailing_zeros(x) ( (x) == 0 ? 64 : __builtin_ctzll(x) )
#endif


namespace {

constexpr int PREVIOUS_VALUES = 128;
constexpr int PREVIOUS_VALUES_LOG2 = 7;
constexpr int THRESHOLD = 6 + PREVIOUS_VALUES_LOG2;
constexpr int SET_LSB = (1 << (THRESHOLD + 1)) - 1;
constexpr int CASE_ZERO_METADATA_LENGTH = PREVIOUS_VALUES_LOG2 + 2;
constexpr int CASE_ONE_METADATA_LENGTH = PREVIOUS_VALUES_LOG2 + 11;

const uint64_t END_SIGN = double_to_long(std::numeric_limits<double>::quiet_NaN());

constexpr uint8_t leading_representation[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

constexpr uint8_t leading_round[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 12, 12, 12, 12, 16, 16, 18, 18, 20, 20, 22, 22, 24, 24, 24,
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24};

// Decoder's leading representation is different
constexpr uint8_t leading_decode[] = {0, 8, 12, 16, 18, 20, 22, 24};


// ChimpEncoderImpl
class ChimpEncoderImpl {
public:
    ChimpEncoderImpl(OutputBitStream& stream) 
        : out(stream), first(true), stored_leading_zeros(std::numeric_limits<int>::max()),
          index(0), current(0)
    {
        stored_values.resize(PREVIOUS_VALUES);
        indices.resize(1 << (THRESHOLD + 1), 0);
    }

    void add_value(double value) {
        uint64_t long_val = double_to_long(value);
        if (first) {
            write_first(long_val);
            first = false;
        } else {
            compress_value(long_val);
        }
    }
    
    void close() {
        add_value(std::numeric_limits<double>::quiet_NaN());
        // Your original padding logic, which is correct for your framework.
        size_t bits_written = out.get_byte_count() * 8 + (out.get_bit_position() % 8);
        int padding = (8 - (bits_written % 8)) % 8;
        if (padding > 0) {
            out.write_bits(0, padding);
        }
    }

private:
    void write_first(uint64_t value) {
        stored_values[current] = value;
        out.write_bits(value, 64);
        indices[static_cast<int>(value & SET_LSB)] = index;
    }

    void compress_value(uint64_t value) {
        // Direct translation of LongChimpEncoder.compressValue
        int key = static_cast<int>(value & SET_LSB);
        uint64_t xor_val;
        int previous_index;
        int trailing_zeros = 0;
        int curr_index = indices[key];

        if ((index - curr_index) < PREVIOUS_VALUES) {
            uint64_t temp_xor = value ^ stored_values[curr_index % PREVIOUS_VALUES];
            trailing_zeros = count_trailing_zeros(temp_xor);
            if (trailing_zeros > THRESHOLD) {
                previous_index = curr_index % PREVIOUS_VALUES;
                xor_val = temp_xor;
            } else {
                previous_index = index % PREVIOUS_VALUES;
                xor_val = stored_values[previous_index] ^ value;
            }
        } else {
            previous_index = index % PREVIOUS_VALUES;
            xor_val = stored_values[previous_index] ^ value;
        }
        
        if (xor_val == 0) {
            // case 00:
            out.write_bits(previous_index, CASE_ZERO_METADATA_LENGTH);
            stored_leading_zeros = 64 + 1;
        } else {
            int leading_zeros = leading_round[count_leading_zeros(xor_val)];
            if (trailing_zeros > THRESHOLD) {
                // case 01:
                int significant_bits = 64 - leading_zeros - trailing_zeros;
                uint64_t temp = 512LL * (PREVIOUS_VALUES + previous_index) + 64LL * leading_representation[leading_zeros] + significant_bits;
                out.write_bits(temp, CASE_ONE_METADATA_LENGTH);
                out.write_bits(xor_val >> trailing_zeros, significant_bits);
                stored_leading_zeros = 64 + 1;
            } else if (leading_zeros == stored_leading_zeros) {
                // case 10: (Note: TSFile writes `writeBit(out); skipBit(out);` which is `10` in reverse)
                out.write_bits(0b10, 2); 
                int significant_bits = 64 - leading_zeros;
                out.write_bits(xor_val, significant_bits);
            } else {
                // case 11:
                stored_leading_zeros = leading_zeros;
                int significant_bits = 64 - leading_zeros;
                out.write_bits(24LL + leading_representation[leading_zeros], 5);
                out.write_bits(xor_val, significant_bits);
            }
        }
        
        current = (current + 1) % PREVIOUS_VALUES;
        stored_values[current] = value;
        index++;
        indices[key] = index;
    }

    OutputBitStream& out;
    bool first;
    int stored_leading_zeros;
    
    std::vector<uint64_t> stored_values;
    std::vector<int> indices;
    int index;
    int current;
};


// ChimpDecoderImpl 
class ChimpDecoderImpl {
public:
    ChimpDecoderImpl(InputBitStream& stream) 
        : in(stream), first(true), stored_leading_zeros(std::numeric_limits<int>::max()), 
          stored_trailing_zeros(0), current(0)
    {
        stored_values.resize(PREVIOUS_VALUES);
    }

    double read_value() {
        uint64_t value;
        if (first) {
            value = in.read_bits(64);
            stored_values[current] = value;
            first = false;
        } else {
            value = next_value();
        }
        
        if (value == END_SIGN) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return long_to_double(value);
    }

private:
    uint64_t next_value() {
        uint64_t control_bits = in.read_bits(1);
        uint64_t value;
        
        if (control_bits == 0) {
            control_bits = in.read_bits(1);
            if (control_bits == 0) {
                // case 00
                int previous_index = in.read_bits(PREVIOUS_VALUES_LOG2);
                value = stored_values[previous_index];
                current = (current + 1) % PREVIOUS_VALUES;
                stored_values[current] = value;
            } else {
                // case 01
                int fill = PREVIOUS_VALUES_LOG2 + 9;
                uint64_t temp = in.read_bits(fill);
                int index = (temp >> (fill - PREVIOUS_VALUES_LOG2)) & ((1 << PREVIOUS_VALUES_LOG2) - 1);
                fill -= PREVIOUS_VALUES_LOG2;
                stored_leading_zeros = leading_decode[(temp >> (fill - 3)) & ((1 << 3) - 1)];
                fill -= 3;
                int significant_bits = temp & ((1 << 6) - 1);
                
                value = stored_values[index];
                if (significant_bits == 0) {
                    significant_bits = 64;
                }
                stored_trailing_zeros = 64 - significant_bits - stored_leading_zeros;
                
                uint64_t xor_part = in.read_bits(64 - stored_leading_zeros - stored_trailing_zeros);
                xor_part <<= stored_trailing_zeros;
                value ^= xor_part;
                
                current = (current + 1) % PREVIOUS_VALUES;
                stored_values[current] = value;
            }
        } else {
            control_bits = in.read_bits(1);
            if (control_bits == 0) {
                // case 10
                uint64_t xor_part = in.read_bits(64 - stored_leading_zeros);
                value = stored_values[current] ^ xor_part;
                current = (current + 1) % PREVIOUS_VALUES;
                stored_values[current] = value;
            } else {
                // case 11
                stored_leading_zeros = leading_decode[in.read_bits(3)];
                uint64_t xor_part = in.read_bits(64 - stored_leading_zeros);
                value = stored_values[current] ^ xor_part;
                current = (current + 1) % PREVIOUS_VALUES;
                stored_values[current] = value;
            }
        }
        return value;
    }

    InputBitStream& in;
    bool first;
    int stored_leading_zeros;
    int stored_trailing_zeros;
    
    std::vector<uint64_t> stored_values;
    int current;
};

} // end anonymous namespace

// Public Interface
std::vector<uint8_t> Chimp::encode(const std::vector<double>& values) {
    if (values.empty()) return {};
    OutputBitStream out;
    ChimpEncoderImpl encoder(out);
    for (double val : values) {
        encoder.add_value(val);
    }
    encoder.close();
    return out.get_bytes();
}

std::vector<double> Chimp::decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    InputBitStream in(compressed_data);
    ChimpDecoderImpl decoder(in);
    std::vector<double> values;
    
    while(in.has_more()) {
        try {
            double val = decoder.read_value();
            if (std::isnan(val)) {
                break;
            }
            values.push_back(val);
        } catch (const std::out_of_range& e) {
            break;
        }
    }
    return values;
}