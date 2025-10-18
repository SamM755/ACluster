#ifndef BIT_STREAM_H
#define BIT_STREAM_H

#include <vector>
#include <cstdint>
#include <stdexcept>

class OutputBitStream {
public:
    OutputBitStream();

    void write_bit(bool bit);
    void write_bits(uint64_t value, int num_bits);

    const std::vector<uint8_t>& get_buffer() const;
    std::vector<uint8_t> get_bytes(); 
    size_t get_byte_count() const;
    int get_bit_position() const { return bit_position_; } 

    void write_bits_padded(uint64_t value, int bits);
    void write_unsigned_varint(uint32_t value);

    void flush();
    void write_bytes(const std::vector<uint8_t>& bytes);

private:
    std::vector<uint8_t> buffer_;
    uint8_t current_byte_;
    int bit_position_;
};

class InputBitStream {
public:
    InputBitStream(const std::vector<uint8_t>& buffer);

    bool read_bit();
    uint64_t read_bits(int num_bits);
    bool has_more() const;

    uint64_t read_bits_padded(int bits);
    uint32_t read_unsigned_varint();

    std::vector<uint8_t> read_remaining_bytes();

private:
    const std::vector<uint8_t>& buffer_;
    size_t byte_index_;
    int bit_position_;
};

#endif // BIT_STREAM_H