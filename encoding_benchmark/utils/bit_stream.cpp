#include "bit_stream.h"
#include "debug.h"

// --- OutputBitStream Implementation ---

OutputBitStream::OutputBitStream() : current_byte_(0), bit_position_(0) {
    buffer_.reserve(1024);
}

void OutputBitStream::write_bit(bool bit) {
    if (bit) {
        current_byte_ |= (1 << (7 - bit_position_));
    }
    bit_position_++;
    if (bit_position_ == 8) {
        buffer_.push_back(current_byte_);
        current_byte_ = 0;
        bit_position_ = 0;
    }
}


void OutputBitStream::write_bits(uint64_t value, int num_bits) {
    DEBUG_BITS("Write", value, num_bits);
    if (num_bits > 64) {
        throw std::invalid_argument("Cannot write more than 64 bits at a time.");
    }
    for (int i = num_bits - 1; i >= 0; --i) {
        write_bit((value >> i) & 1);
    }
}

const std::vector<uint8_t>& OutputBitStream::get_buffer() const {
    return buffer_;
}


std::vector<uint8_t> OutputBitStream::get_bytes() {
    std::vector<uint8_t> result = buffer_;
    if (bit_position_ > 0) {
        result.push_back(current_byte_);
    }
    return result;
}


size_t OutputBitStream::get_byte_count() const {
    return buffer_.size() + (bit_position_ > 0 ? 1 : 0);
}

void OutputBitStream::write_bits_padded(uint64_t value, int bits) {
    // This function writes 'bits' LSBs of value, byte by byte.
    // It's different from write_bits which is bit-by-bit.
    for (int i = 0; i < bits; i += 8) {
        int bits_to_write = std::min(8, bits - i);
        uint8_t byte = (value >> i) & ((1 << bits_to_write) - 1);
        write_bits(byte, bits_to_write);
    }
}

void OutputBitStream::write_unsigned_varint(uint32_t value) {
    while (value >= 0x80) {
        write_bits((value & 0x7F) | 0x80, 8);
        value >>= 7;
    }
    write_bits(value & 0x7F, 8);
}

void OutputBitStream::write_bytes(const std::vector<uint8_t>& bytes) {
    if (bit_position_ == 0) {
        
        flush(); 
        buffer_.insert(buffer_.end(), bytes.begin(), bytes.end());
    } else {
        
        for (uint8_t byte : bytes) {
            write_bits(byte, 8);
        }
    }
}


void OutputBitStream::flush() {
    if (bit_position_ > 0) {
        buffer_.push_back(current_byte_);
        current_byte_ = 0;
        bit_position_ = 0;
    }
}

std::vector<uint8_t> InputBitStream::read_remaining_bytes() {
    std::vector<uint8_t> remaining;
    
    
    if (bit_position_ != 0) {
        
        while (has_more()) {
            try {
                remaining.push_back(static_cast<uint8_t>(read_bits(8)));
            } catch (const std::out_of_range&) {
                
                break;
            }
        }
    } else {
        
        if (byte_index_ < buffer_.size()) {
            remaining.insert(remaining.end(), buffer_.begin() + byte_index_, buffer_.end());
            byte_index_ = buffer_.size(); 
        }
    }
    return remaining;
}

// --- InputBitStream Implementation ---
InputBitStream::InputBitStream(const std::vector<uint8_t>& buffer)
    : buffer_(buffer), byte_index_(0), bit_position_(0) {}

bool InputBitStream::read_bit() {
    if (byte_index_ >= buffer_.size()) {
        throw std::out_of_range("Reading past end of buffer.");
    }
    bool bit = (buffer_[byte_index_] >> (7 - bit_position_)) & 1;
    bit_position_++;
    if (bit_position_ == 8) {
        byte_index_++;
        bit_position_ = 0;
    }
    return bit;
}

uint64_t InputBitStream::read_bits(int num_bits) {
    if (num_bits > 64) {
        throw std::invalid_argument("Cannot read more than 64 bits at a time.");
    }
    uint64_t value = 0;
    for (int i = 0; i < num_bits; ++i) {
        value <<= 1;
        if (read_bit()) {
            value |= 1;
        }
    }
    DEBUG_BITS("Read ", value, num_bits);
    return value;
}

bool InputBitStream::has_more() const {
    return byte_index_ < buffer_.size();
}

uint64_t InputBitStream::read_bits_padded(int bits) {
    uint64_t value = 0;
    for (int i = 0; i < bits; i += 8) {
        int bits_to_read = std::min(8, bits - i);
        uint64_t byte = read_bits(bits_to_read);
        value |= (byte << i);
    }
    return value;
}

uint32_t InputBitStream::read_unsigned_varint() {
    uint32_t value = 0;
    int shift = 0;
    uint8_t byte;
    do {
        byte = read_bits(8);
        value |= (byte & 0x7F) << shift;
        shift += 7;
    } while ((byte & 0x80) != 0);
    return value;
}