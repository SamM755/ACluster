#ifndef TYPE_CONVERTER_H
#define TYPE_CONVERTER_H

#include <cstdint>
#include <cstring> // For std::memcpy


inline uint64_t double_to_long(double val) {
    uint64_t res;
    
    std::memcpy(&res, &val, sizeof(double));
    return res;
}


inline double long_to_double(uint64_t val) {
    double res;
    std::memcpy(&res, &val, sizeof(double));
    return res;
}

#endif // TYPE_CONVERTER_H