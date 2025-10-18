// in algorithms/gorilla.cpp

#include "gorilla.h"
#include "../utils/bit_stream.h"
#include "../utils/debug.h"
#include "../utils/type_converter.h" // Assuming double_to_long is here
#include <cmath>
#include <limits>
#include <cstring>

// Platform-specific functions for counting zeros
#ifdef _MSC_VER
#include <intrin.h>
#define count_leading_zeros(x) ( (x) == 0 ? 64 : __lzcnt64(x) )
#define count_trailing_zeros(x) ( (x) == 0 ? 64 : _tzcnt_u64(x) )
#else
#define count_leading_zeros(x) ( (x) == 0 ? 64 : __builtin_clzll(x) )
#define count_trailing_zeros(x) ( (x) == 0 ? 64 : __builtin_ctzll(x) )
#endif


namespace { // Anonymous namespace

const uint64_t END_SIGN = double_to_long(std::numeric_limits<double>::quiet_NaN());

class GorillaEncoderImpl {
public:
    GorillaEncoderImpl(OutputBitStream& stream) 
        : out(stream), stored_leading_zeros(std::numeric_limits<int>::max()), 
          stored_trailing_zeros(0), stored_val(0), first(true) {}

    void addValue(double value) {
        uint64_t long_val = double_to_long(value);
        if (first) {
            writeFirst(long_val);
            first = false;
        } else {
            compressValue(long_val);
        }
    }

    void close() {
        addValue(std::numeric_limits<double>::quiet_NaN());
        size_t bits_written = out.get_byte_count() * 8 + (out.get_bit_position() % 8);
        int padding = (8 - (bits_written % 8)) % 8;
        if (padding > 0) {
            out.write_bits(0, padding);
        }
    }

private:
    void writeFirst(uint64_t value) {
        stored_val = value;
        out.write_bits(stored_val, 64);
    }

    void compressValue(uint64_t value) {
        uint64_t xor_val = stored_val ^ value;

        if (xor_val == 0) {
            // Case 1: Value is the same, write a single '0' bit.
            out.write_bit(false);
        } else {
            int leading_zeros = count_leading_zeros(xor_val);
            int trailing_zeros = count_trailing_zeros(xor_val);

            // Using stored_leading_zeros != max() to check if it's the second value
            if (stored_leading_zeros != std::numeric_limits<int>::max() &&
                leading_zeros >= stored_leading_zeros && 
                trailing_zeros >= stored_trailing_zeros) {
                // Case 2: Control bits '10'. Use stored LZ/TZ.
                out.write_bits(0b10, 2);
                int significant_bits = 64 - stored_leading_zeros - stored_trailing_zeros;
                out.write_bits(xor_val >> stored_trailing_zeros, significant_bits);
            } else {
                // Case 3: Control bits '11'. Write new LZ/TZ.
                out.write_bits(0b11, 2);
                
                out.write_bits(leading_zeros, 6); 
                
                int significant_bits = 64 - leading_zeros - trailing_zeros;
                if (significant_bits >= 64) { // handle case where sig_bits == 64
                    out.write_bits(0, 6);
                    out.write_bits(xor_val, 64);
                } else {
                    out.write_bits(significant_bits, 6);
                    out.write_bits(xor_val >> trailing_zeros, significant_bits);
                }
                
                // Update stored LZ/TZ
                stored_leading_zeros = leading_zeros;
                stored_trailing_zeros = trailing_zeros;
            }
        }
        stored_val = value;
    }

    OutputBitStream& out;
    int stored_leading_zeros;
    int stored_trailing_zeros;
    uint64_t stored_val;
    bool first;
};


class GorillaDecompressorImpl {
public:
    GorillaDecompressorImpl(InputBitStream& stream)
        : in(stream), 
          stored_leading_zeros(0), 
          stored_trailing_zeros(0), 
          stored_val(0), 
          first(true) {}

    double readValue() {
        uint64_t value;
        if (first) {
            value = in.read_bits(64);
            first = false;
        } else {
            value = nextValue();
        }
        
        if (value == END_SIGN) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        stored_val = value;
        return long_to_double(stored_val);
    }

private:
    uint64_t nextValue() {
        uint64_t value;
        
        if (!in.read_bit()) {
            // Case 1: Value is the same.
            return stored_val;
        }
        
        if (!in.read_bit()) {
            // Case 2: Control bits '10'. Use stored LZ/TZ.
            int significant_bits = 64 - stored_leading_zeros - stored_trailing_zeros;
            uint64_t xor_part = in.read_bits(significant_bits);
            xor_part <<= stored_trailing_zeros;
            value = stored_val ^ xor_part;
        } else {
            // Case 3: Control bits '11'. Read new LZ/TZ.
            // Read 6 bits for leading zeros.
            stored_leading_zeros = in.read_bits(6);
            
            int significant_bits = in.read_bits(6);
            if (significant_bits == 0) {
                significant_bits = 64;
            }
            stored_trailing_zeros = 64 - stored_leading_zeros - significant_bits;
            
            uint64_t xor_part = in.read_bits(significant_bits);
            xor_part <<= stored_trailing_zeros;
            value = stored_val ^ xor_part;
        }
        
        return value;
    }

    InputBitStream& in;
    int stored_leading_zeros;
    int stored_trailing_zeros;
    uint64_t stored_val;
    bool first;
};

} // end anonymous namespace


std::vector<uint8_t> Gorilla::encode(const std::vector<double>& values) {
    if (values.empty()) return {};
    OutputBitStream out;
    GorillaEncoderImpl encoder(out);
    for (double val : values) {
        encoder.addValue(val);
    }
    encoder.close();
    return out.get_bytes();
}

std::vector<double> Gorilla::decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    InputBitStream in(compressed_data);
    GorillaDecompressorImpl decoder(in);
    std::vector<double> values;

    while(in.has_more()) {
        try {
            double val = decoder.readValue();
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