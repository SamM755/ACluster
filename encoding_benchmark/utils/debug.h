#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <iomanip>
#include <bitset>

constexpr bool DEBUG_ENABLED = false; 


#define DEBUG_LOG(x) \
    do { \
        if (DEBUG_ENABLED) { \
            std::cout << "[DEBUG] " << x << std::endl; \
        } \
    } while (0)


#define DEBUG_BITS(action, value, num_bits) \
    do { \
        if (DEBUG_ENABLED) { \
            if (num_bits > 0) { \
                std::cout << "[DEBUG] " << action << " " << num_bits << " bits: " \
                          << std::bitset<64>(value).to_string().substr(64 - num_bits) \
                          << " (val=" << value << ")" << std::endl; \
            } \
        } \
    } while (0)

#endif // DEBUG_H