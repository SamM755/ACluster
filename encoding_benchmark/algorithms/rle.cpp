// in algorithms/rle.cpp

#include "rle.h"
#include "../utils/bit_stream.h"
#include "../utils/type_converter.h"
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <numeric>




namespace { // Anonymous namespace

// =====================================================================================
//                 PRE-PROCESSING: Double <-> Scaled Int64
// =====================================================================================

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

inline int get_bit_width(uint64_t n) {
    if (n == 0) return 1;
#ifdef _MSC_VER
    unsigned long index;
    _BitScanReverse64(&index, n);
    return static_cast<int>(index) + 1;
#else
    return 64 - __builtin_clzll(n);
#endif
}

class LongPacker {
private:
    int bit_width;
public:
    LongPacker(int width) : bit_width(width) {}

    void pack8Values(const std::vector<int64_t>& in_buf, int offset, std::vector<uint8_t>& out_buf) {
        uint64_t temp[8];
        for (int i = 0; i < 8; i++) temp[i] = static_cast<uint64_t>(in_buf[offset + i]);
        
        out_buf.assign(bit_width, 0);
        for (int i = 0; i < bit_width; i++) {
            out_buf[i] = (uint8_t)(
                ((temp[0] >> i) & 1) << 7 | ((temp[1] >> i) & 1) << 6 |
                ((temp[2] >> i) & 1) << 5 | ((temp[3] >> i) & 1) << 4 |
                ((temp[4] >> i) & 1) << 3 | ((temp[5] >> i) & 1) << 2 |
                ((temp[6] >> i) & 1) << 1 | ((temp[7] >> i) & 1));
        }
    }

    void unpackAllValues(const std::vector<uint8_t>& in_buf, int length, std::vector<int64_t>& out_buf) {
        int in_offset = 0;
        int out_offset = 0;
        int group_count = length / bit_width;

        for (int i = 0; i < group_count; i++) {
            unpack8Values(in_buf, in_offset, out_buf, out_offset);
            in_offset += bit_width;
            out_offset += 8;
        }
    }
private:
    void unpack8Values(const std::vector<uint8_t>& in_buf, int in_offset, std::vector<int64_t>& out_buf, int out_offset) {
        for (int i = 0; i < 8; i++) out_buf[out_offset + i] = 0;
        for (int i = 0; i < bit_width; i++) {
            uint8_t byte = in_buf[in_offset + i];
            for (int j = 0; j < 8; j++) {
                out_buf[out_offset + j] |= (int64_t)(((byte >> (7 - j)) & 1)) << i;
            }
        }
    }
};

// =====================================================================================
//                 TSFile-style RLE/Bit-Packing Encoder
// =====================================================================================
class RleEncoderImpl {
private:
    static constexpr int RLE_MIN_REPEATED_NUM = 8;
    int bit_width;
    LongPacker packer;

    // State machine variables
    std::vector<int64_t> buffered_values; // For bit-packing
    int num_buffered_values;
    int64_t pre_value;
    int repeat_count;

    void writeRleRun(OutputBitStream& out) {
        out.write_unsigned_varint(repeat_count << 1);
        write_little_endian(out, pre_value, (bit_width + 7) / 8);
        num_buffered_values = 0;
        repeat_count = 0;
    }

    void convertBuffer(OutputBitStream& out) {
        std::vector<uint8_t> bytes(bit_width);
        packer.pack8Values(buffered_values, 0, bytes);
        
        out.write_unsigned_varint((RLE_MIN_REPEATED_NUM << 1) | 1);
        for(uint8_t byte : bytes) {
            out.write_bits(byte, 8);
        }
        num_buffered_values = 0;
    }

    void write_little_endian(OutputBitStream& out, int64_t value, int num_bytes) {
        for(int i = 0; i < num_bytes; ++i) {
            out.write_bits((value >> (i*8)) & 0xFF, 8);
        }
    }

public:
    RleEncoderImpl(int width)
        : bit_width(width), packer(width), num_buffered_values(0), pre_value(0), repeat_count(0)
    {
        buffered_values.resize(RLE_MIN_REPEATED_NUM, 0);
    }
    
    // This is the precise C++ translation of RleEncoder.encodeValue(T value)
    void encode(int64_t value, OutputBitStream& out) {
        if (repeat_count == 0 && num_buffered_values == 0) {
            // First value
            pre_value = value;
            repeat_count = 1;
            return;
        }

        if (value == pre_value) {
            repeat_count++;
            if (repeat_count == RLE_MIN_REPEATED_NUM) {
                // We have a full RLE run starting. If there were any
                // previous values in the bit-pack buffer, write them out first.
                if (num_buffered_values > 0) {
                     out.write_unsigned_varint((num_buffered_values << 1) | 1);
                     for(int i=0; i<num_buffered_values; ++i) write_little_endian(out, buffered_values[i], (bit_width+7)/8);
                }
                num_buffered_values = 0;
            }
        } else {
            // New value arrived. Handle the previous sequence.
            if (repeat_count >= RLE_MIN_REPEATED_NUM) {
                // Previous sequence was a long RLE run
                writeRleRun(out);
            } else {
                // Previous sequence was a short run, treat as bit-packed
                for(int i = 0; i < repeat_count; ++i) {
                    buffered_values[num_buffered_values++] = pre_value;
                    if (num_buffered_values == RLE_MIN_REPEATED_NUM) {
                        convertBuffer(out);
                    }
                }
            }
            // Start new sequence with the current value
            pre_value = value;
            repeat_count = 1;
        }
    }
    
    void flush(OutputBitStream& out) {
        if (repeat_count >= RLE_MIN_REPEATED_NUM) {
            writeRleRun(out);
        } else if (repeat_count > 0) {
            // Not a long enough run, so add the repeated values to the buffer
             for(int i = 0; i < repeat_count; ++i) {
                buffered_values[num_buffered_values++] = pre_value;
                if (num_buffered_values == RLE_MIN_REPEATED_NUM) {
                    convertBuffer(out);
                }
            }
        }
        
        // Write any remaining values in the bit-packing buffer
        if (num_buffered_values > 0) {
            out.write_unsigned_varint((num_buffered_values << 1) | 1);
            for (int i = 0; i < num_buffered_values; ++i) {
                write_little_endian(out, buffered_values[i], (bit_width + 7) / 8);
            }
        }
        num_buffered_values = 0;
        repeat_count = 0;
    }
};

class RleDecoderImpl {
private:
    InputBitStream& in;
    int bit_width;
    enum Mode { RLE, BIT_PACKED };
    Mode mode;
    int current_count;
    int64_t rle_value;
    std::vector<int64_t> bit_packed_values;
    int bit_packed_idx;
    LongPacker packer;

    void readNextGroup() {
        uint32_t header = in.read_unsigned_varint();
        mode = ((header & 1) == 0) ? RLE : BIT_PACKED;
        current_count = header >> 1;
        
        if (mode == RLE) {
            rle_value = read_little_endian((bit_width + 7) / 8);
        } else {
            bit_packed_values.assign(current_count, 0);
            int num_groups = current_count / 8;
            if (num_groups > 0) {
                int bytes_to_read = num_groups * bit_width;
                std::vector<uint8_t> bytes(bytes_to_read);
                for (int i = 0; i < bytes_to_read; ++i) bytes[i] = in.read_bits(8);
                packer.unpackAllValues(bytes, bytes_to_read, bit_packed_values);
            }
            int remainder = current_count % 8;
            for (int i = 0; i < remainder; ++i) {
                bit_packed_values[num_groups * 8 + i] = read_little_endian((bit_width + 7) / 8);
            }
            bit_packed_idx = 0;
        }
    }

    int64_t read_little_endian(int num_bytes) {
        int64_t val = 0;
        for(int i = 0; i < num_bytes; ++i) {
            val |= (int64_t)in.read_bits(8) << (i*8);
        }
        return val;
    }

public:
    RleDecoderImpl(InputBitStream& stream, int width) 
        : in(stream), bit_width(width), current_count(0), bit_packed_idx(0), packer(width) {}

    int64_t read() {
        if (current_count == 0) readNextGroup();
        current_count--;
        if (mode == RLE) {
            return rle_value;
        } else {
            return bit_packed_values[bit_packed_idx++];
        }
    }
};
}

// =====================================================================================
//                                  Public Rle Class
// =====================================================================================
std::vector<uint8_t> Rle::encode(const std::vector<double>& data) {
    if (data.empty()) return {};
    ScalingResult scaling_res = scale_doubles(data);
    const auto& scaled_data = scaling_res.scaled_data;
    uint64_t max_val_unsigned = 0;
    for (int64_t val : scaled_data) {
        uint64_t zigzag_val = (static_cast<uint64_t>(val) << 1) ^ (val >> 63);
        if (zigzag_val > max_val_unsigned) max_val_unsigned = zigzag_val;
    }
    int bit_width = get_bit_width(max_val_unsigned);
    OutputBitStream out;
    out.write_bits(double_to_long(scaling_res.scaling_factor), 64);
    out.write_bits(bit_width, 8);
    out.write_bits(data.size(), 32);
    
    RleEncoderImpl encoder(bit_width);
    if (!scaled_data.empty()) {
        for (int64_t val : scaled_data) {
            uint64_t zigzag_val = (static_cast<uint64_t>(val) << 1) ^ (val >> 63);
            encoder.encode(zigzag_val, out);
        }
        encoder.flush(out);
    }
    
    return out.get_bytes();
}

std::vector<double> Rle::decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    InputBitStream in(compressed_data);
    double scaling_factor = long_to_double(in.read_bits(64));
    int bit_width = in.read_bits(8);
    uint32_t num_values = in.read_bits(32);
    
    std::vector<double> decoded_values;
    decoded_values.reserve(num_values);
    
    RleDecoderImpl decoder(in, bit_width);
    
    for (uint32_t i = 0; i < num_values; ++i) {
        uint64_t zigzag_val = decoder.read();
        int64_t scaled_val = (zigzag_val >> 1) ^ (-(zigzag_val & 1));
        decoded_values.push_back(static_cast<double>(scaled_val) / scaling_factor);
    }
    return decoded_values;
}