#include "alp.h"
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <iostream>

namespace {

// ==========================================
// 1. ALP Constants
// ==========================================
constexpr size_t VECTOR_SIZE = 1024;
constexpr uint8_t MAX_EXPONENT = 18;

// FRAC_ARR = 10^{-f} (ALP calls this i_F10)
constexpr double FRAC_ARR[] = {
    1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001,
    0.0000000001, 0.00000000001, 0.000000000001, 0.0000000000001, 0.00000000000001,
    0.000000000000001, 0.0000000000000001, 0.00000000000000001, 0.000000000000000001
};

// EXP_ARR = 10^e (ALP calls this F10)
constexpr double EXP_ARR[] = {
    1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0, 100000000.0,
    1000000000.0, 10000000000.0, 100000000000.0, 1000000000000.0, 10000000000000.0,
    100000000000000.0, 1000000000000000.0, 10000000000000000.0, 100000000000000000.0,
    1000000000000000000.0
};

// ==========================================
// 2. Robust Bit Manipulation
// ==========================================

class BitWriter {
public:
    std::vector<uint8_t> data;
    uint64_t buffer = 0;
    int bits_in_buffer = 0;

    void write_bits(uint64_t value, uint8_t width) {
        if (width == 0) return;
        
        value &= (~0ULL >> (64 - width)); // Mask garbage

        if (bits_in_buffer + width <= 64) {
            // Fits in current buffer
            buffer |= (value << bits_in_buffer);
            bits_in_buffer += width;
        } else {
            // Split across buffer boundary
            uint32_t first_chunk_size = 64 - bits_in_buffer;
            uint64_t first_chunk_mask = (~0ULL >> (64 - first_chunk_size));
            
            buffer |= ((value & first_chunk_mask) << bits_in_buffer);
            flush_full_buffer();
            
            // Write remaining bits
            uint32_t second_chunk_size = width - first_chunk_size;
            buffer = (value >> first_chunk_size);
            bits_in_buffer = second_chunk_size;
        }

        if (bits_in_buffer >= 64) flush_full_buffer(); 
    }

    void flush_full_buffer() {
        for(int i=0; i<8; ++i) {
            data.push_back(static_cast<uint8_t>(buffer & 0xFF));
            buffer >>= 8;
        }
        bits_in_buffer = 0;
        buffer = 0; // Clean slate, though shifted out anyway
    }

    void flush_final() {
        while (bits_in_buffer > 0) {
            data.push_back(static_cast<uint8_t>(buffer & 0xFF));
            buffer >>= 8;
            bits_in_buffer -= 8;
        }
        buffer = 0;
        bits_in_buffer = 0;
    }

    // Aligned writes
    void write_u8(uint8_t v) { flush_final(); data.push_back(v); }
    void write_u16(uint16_t v) { flush_final(); data.push_back(v & 0xFF); data.push_back(v >> 8); }
    void write_u64(uint64_t v) { 
        flush_final(); 
        for(int i=0; i<8; ++i) data.push_back((v >> (i*8)) & 0xFF);
    }
    void write_double(double v) {
        uint64_t tmp;
        std::memcpy(&tmp, &v, sizeof(double));
        write_u64(tmp);
    }
};

class BitReader {
public:
    const std::vector<uint8_t>& data;
    size_t byte_pos = 0;
    uint64_t buffer = 0;
    int bits_in_buffer = 0;

    BitReader(const std::vector<uint8_t>& d) : data(d) {}

    uint64_t read_bits(uint8_t width) {
        if (width == 0) return 0;

        // Ensure enough bits
        while (bits_in_buffer < width) {
            if (byte_pos < data.size()) {
                uint64_t next_byte = data[byte_pos++];
                buffer |= (next_byte << bits_in_buffer);
                bits_in_buffer += 8;
            } else {
                bits_in_buffer += 8; // Padding with 0s conceptually
            }
        }

        uint64_t mask = (~0ULL >> (64 - width));
        uint64_t result = buffer & mask;
        
        buffer >>= width;
        bits_in_buffer -= width;
        return result;
    }

    void align() {
        buffer = 0;
        bits_in_buffer = 0;
    }

    uint8_t read_u8() { align(); if(byte_pos >= data.size()) return 0; return data[byte_pos++]; }
    
    uint16_t read_u16() { 
        align(); 
        if(byte_pos+1 >= data.size()) { byte_pos+=2; return 0; }
        uint16_t v = data[byte_pos] | (data[byte_pos+1] << 8); 
        byte_pos += 2; return v; 
    }

    uint64_t read_u64() {
        align();
        if(byte_pos+7 >= data.size()) { byte_pos+=8; return 0; }
        uint64_t v = 0;
        for(int i=0; i<8; ++i) v |= (static_cast<uint64_t>(data[byte_pos+i]) << (i*8));
        byte_pos += 8; return v;
    }

    double read_double() {
        uint64_t v = read_u64();
        double d;
        std::memcpy(&d, &v, sizeof(double));
        return d;
    }
};

// ==========================================
// 3. Core ALP Logic
// ==========================================

// CRITICAL FIX: Use std::llround instead of cast
// Formula: value * 10^exp * 10^-fac
inline int64_t encode_val(double value, uint8_t fac_idx, uint8_t exp_idx) {
    return std::llround(value * EXP_ARR[exp_idx] * FRAC_ARR[fac_idx]);
}

inline double decode_val(int64_t encoded, uint8_t fac_idx, uint8_t exp_idx) {
    // Formula: encoded * 10^fac * 10^-exp
    return static_cast<double>(encoded) * EXP_ARR[fac_idx] * FRAC_ARR[exp_idx];
}

uint8_t get_bit_width(uint64_t val) {
    if (val == 0) return 0;
    // __builtin_clzll is undefined for 0, handled above
#ifdef _MSC_VER
    unsigned long idx;
    if (_BitScanReverse64(&idx, val)) return (uint8_t)(idx + 1);
    return 0;
#else
    return (uint8_t)(64 - __builtin_clzll(val));
#endif
}

struct AlpBlock {
    bool is_raw = false; // Flag to skip compression if ineffective
    uint8_t exp;
    uint8_t fac;
    uint8_t bit_width;
    uint16_t count;
    std::vector<uint16_t> exc_pos;
    std::vector<double> exc_val;
    int64_t base_value; 
    std::vector<int64_t> encoded_ints;
};

// Analyze a block to find best parameters
AlpBlock compress_block(const double* data, size_t n) {
    AlpBlock best_block;
    best_block.is_raw = true; // Default to raw until we find good compression
    double min_bits_total = std::numeric_limits<double>::max();
    
    // Raw size reference (64 bits per value)
    double raw_size_bits = (double)n * 64.0;

    // Search space
    for (uint8_t exp = 0; exp <= MAX_EXPONENT; ++exp) {
        for (uint8_t fac = 0; fac <= exp; ++fac) {
            
            int64_t min_v = std::numeric_limits<int64_t>::max();
            int64_t max_v = std::numeric_limits<int64_t>::min();
            int exceptions = 0;
            
            // Temporary storage for checking range
            // We only need to know min/max of SUCCESSFUL conversions for Frame-Of-Reference
            std::vector<int64_t> successful_values;
            successful_values.reserve(n);

            for (size_t i = 0; i < n; ++i) {
                int64_t enc = encode_val(data[i], fac, exp);
                double dec = decode_val(enc, fac, exp);
                
                // CRITICAL FIX: Compare with tolerance for round-trip safety
                // ALP ensures lossless by storing exception if exact match fails
                if (dec != data[i]) {
                    exceptions++;
                } else {
                    if (enc < min_v) min_v = enc;
                    if (enc > max_v) max_v = enc;
                    successful_values.push_back(enc);
                }
            }

            uint8_t bw = 0;
            int64_t current_base = 0;
            
            if (exceptions == n) {
                // All failed, this combo is useless
                bw = 0; 
                current_base = 0;
            } else {
                current_base = min_v;
                // Calculate max delta
                // Note: We use uint64 for delta calculation to safely handle full int64 range
                // if min_v is negative and max_v positive.
                // However, standard ALP implementation ensures values fit in int64.
                // We strictly cast to u64 for bit width calc.
                uint64_t max_delta = (uint64_t)(max_v - current_base);
                bw = get_bit_width(max_delta);
            }

            // Calculate estimated cost
            // Overhead: 16(pos) + 64(val) = 80 bits per exception
            double current_bits = (double)exceptions * 80.0 + (double)n * bw;
            
            // Header overheads are constant, ignored for comparison
            
            if (current_bits < min_bits_total) {
                min_bits_total = current_bits;
                best_block.exp = exp;
                best_block.fac = fac;
                best_block.bit_width = bw;
                best_block.base_value = current_base;
                best_block.is_raw = false;
            }
        }
    }

    // Heuristic: If compression is barely better than raw, or worse, use raw.
    // 50 bits overhead for headers roughly.
    if (min_bits_total > raw_size_bits || best_block.is_raw) {
        best_block.is_raw = true;
        best_block.count = (uint16_t)n;
        return best_block;
    }

    // Finalize the best block
    best_block.count = (uint16_t)n;
    best_block.encoded_ints.resize(n);
    
    for (size_t i = 0; i < n; ++i) {
        int64_t enc = encode_val(data[i], best_block.fac, best_block.exp);
        double dec = decode_val(enc, best_block.fac, best_block.exp);
        
        if (dec != data[i]) {
            best_block.exc_pos.push_back((uint16_t)i);
            best_block.exc_val.push_back(data[i]);
            // ALP Logic: Store the 'base_value' in the integer slot for exceptions
            // so the delta is 0, minimizing bit width impact.
            best_block.encoded_ints[i] = best_block.base_value; 
        } else {
            best_block.encoded_ints[i] = enc;
        }
    }

    return best_block;
}

} // namespace

// ==========================================
// 4. Interface Implementation
// ==========================================

std::vector<uint8_t> Alp::encode(const std::vector<double>& data) {
    if (data.empty()) return {};

    BitWriter writer;
    
    // Global Header: Total Count (4 bytes)
    uint32_t total_count = (uint32_t)data.size();
    writer.write_u8((total_count >> 0) & 0xFF);
    writer.write_u8((total_count >> 8) & 0xFF);
    writer.write_u8((total_count >> 16) & 0xFF);
    writer.write_u8((total_count >> 24) & 0xFF);

    for (size_t i = 0; i < data.size(); i += VECTOR_SIZE) {
        size_t chunk_size = std::min(VECTOR_SIZE, data.size() - i);
        const double* chunk_ptr = &data[i];

        AlpBlock block = compress_block(chunk_ptr, chunk_size);

        // Serialize Block
        // Mode flag: 1 bit (0 = Raw, 1 = Compressed) - NO, let's just use a special Exp value or Header.
        // Simple Header: [IsRaw(1 byte)]
        
        if (block.is_raw) {
            writer.write_u8(1); // 1 means Raw
            writer.write_u16(block.count);
            for(size_t j=0; j<chunk_size; ++j) {
                writer.write_double(chunk_ptr[j]);
            }
        } else {
            writer.write_u8(0); // 0 means Compressed
            // [Count(16)][Exp(8)][Fac(8)][BW(8)][ExcCount(16)]
            writer.write_u16(block.count);
            writer.write_u8(block.exp);
            writer.write_u8(block.fac);
            writer.write_u8(block.bit_width);
            writer.write_u16((uint16_t)block.exc_pos.size());
            
            // Base Value (64 bits)
            writer.write_u64(static_cast<uint64_t>(block.base_value));

            // Exceptions
            for (uint16_t pos : block.exc_pos) writer.write_u16(pos);
            for (double val : block.exc_val) writer.write_double(val);

            // Packed Integers (FOR)
            for (int64_t val : block.encoded_ints) {
                uint64_t delta = static_cast<uint64_t>(val - block.base_value);
                writer.write_bits(delta, block.bit_width);
            }
        }
    }

    writer.flush_final();
    return writer.data;
}


std::vector<double> Alp::decode(const std::vector<uint8_t>& data) {
    if (data.empty()) return {};

    BitReader reader(data);
    std::vector<double> result;

    uint32_t total_count = 0;
    total_count |= reader.read_u8();
    total_count |= (static_cast<uint32_t>(reader.read_u8()) << 8);
    total_count |= (static_cast<uint32_t>(reader.read_u8()) << 16);
    total_count |= (static_cast<uint32_t>(reader.read_u8()) << 24);
    
    result.resize(total_count);
    size_t current_idx = 0;

    while (current_idx < total_count) {
        uint8_t is_raw = reader.read_u8();
        uint16_t chunk_size = reader.read_u16();

        if (is_raw) {
            for(size_t i=0; i<chunk_size; ++i) {
                result[current_idx + i] = reader.read_double();
            }
        } else {
            uint8_t exp = reader.read_u8();
            uint8_t fac = reader.read_u8();
            uint8_t bw = reader.read_u8();
            uint16_t exc_count = reader.read_u16();

            int64_t base_value = static_cast<int64_t>(reader.read_u64());

            std::vector<uint16_t> exc_pos(exc_count);
            for(int i=0; i<exc_count; ++i) exc_pos[i] = reader.read_u16();

            std::vector<double> exc_val(exc_count);
            for(int i=0; i<exc_count; ++i) exc_val[i] = reader.read_double();

            for (size_t i = 0; i < chunk_size; ++i) {
                uint64_t delta = reader.read_bits(bw);
                int64_t val = base_value + static_cast<int64_t>(delta);
                result[current_idx + i] = decode_val(val, fac, exp);
            }

            for (int i = 0; i < exc_count; ++i) {
                result[current_idx + exc_pos[i]] = exc_val[i];
            }
        }
        current_idx += chunk_size;
    }

    return result;
}