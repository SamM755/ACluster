#ifndef SECONDARY_COMPRESSOR_H
#define SECONDARY_COMPRESSOR_H

#include <vector>
#include <string>
#include <stdexcept>
#include "../utils/miniz.h"

class SecondaryCompressor {
public:
    virtual ~SecondaryCompressor() = default;
    virtual std::vector<uint8_t> compress(const std::vector<uint8_t>& input) = 0;
    virtual std::vector<uint8_t> decompress(const std::vector<uint8_t>& input) = 0;
    virtual std::string get_name() const = 0;
};

class NoneCompressor : public SecondaryCompressor {
public:
    std::vector<uint8_t> compress(const std::vector<uint8_t>& input) override {
        return input; 
    }
    std::vector<uint8_t> decompress(const std::vector<uint8_t>& input) override {
        return input; 
    }
    std::string get_name() const override {
        return ""; 
    }
};

class GzipCompressor : public SecondaryCompressor {
public:
    std::vector<uint8_t> compress(const std::vector<uint8_t>& input) override {
        if (input.empty()) return {};

        uLong sourceLen = (uLong)input.size();
        uLong destLen = compressBound(sourceLen);
        std::vector<uint8_t> output(destLen);

        int status = mz_compress(output.data(), &destLen, input.data(), sourceLen);
        
        if (status != MZ_OK) {
            throw std::runtime_error("GZIP compression failed with error code: " + std::to_string(status));
        }

        output.resize(destLen);
        return output;
    }

    std::vector<uint8_t> decompress(const std::vector<uint8_t>& input) override {
        if (input.empty()) return {};

        uLong sourceLen = (uLong)input.size();
        uLong destLen = sourceLen * 5;
        if (destLen < 1024) destLen = 1024;
        
        std::vector<uint8_t> output(destLen);
        
        while (true) {
            int status = mz_uncompress(output.data(), &destLen, input.data(), sourceLen);
            
            if (status == MZ_OK) {
                output.resize(destLen);
                return output;
            } 
            else if (status == MZ_BUF_ERROR) {
                destLen *= 2;
                output.resize(destLen);
            } 
            else {
                throw std::runtime_error("GZIP decompression failed with error code: " + std::to_string(status));
            }
        }
    }

    std::string get_name() const override {
        return "+GZIP";
    }
};

#endif // SECONDARY_COMPRESSOR_H