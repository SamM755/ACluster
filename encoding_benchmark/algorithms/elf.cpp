#include "elf.h"
#include "../utils/bit_stream.h"
#include "../utils/debug.h"
#include <cmath>
#include <limits>
#include <cstring>
#include <array>
#include <numeric>

inline uint64_t double_to_long(double val) {
    uint64_t res;
    std::memcpy(&res, &val, sizeof(val));
    return res;
}

inline double long_to_double(uint64_t val) {
    double res;
    std::memcpy(&res, &val, sizeof(val));
    return res;
}

#ifdef _MSC_VER
#include <intrin.h>
#define count_leading_zeros(x) (x == 0 ? 64 : __lzcnt64(x))
#define count_trailing_zeros(x) (x == 0 ? 64 : _tzcnt_u64(x))
#else
#define count_leading_zeros(x) (x == 0 ? 64 : __builtin_clzll(x))
#define count_trailing_zeros(x) (x == 0 ? 64 : __builtin_ctzll(x))
#endif

namespace {

// =====================================================================================
//                                    Elf64Utils
// =====================================================================================
namespace Elf64Utils {
    // Forward declarations
    static double get10iP(int i);

    // αlog_2(10) for look-up
    static constexpr std::array<int, 21> f = {0, 4, 7, 10, 14, 17, 20, 24, 27, 30, 34, 37, 40, 44, 47, 50, 54, 57, 60, 64, 67};
    static constexpr std::array<double, 21> map10iP = {1.0, 1.0E1, 1.0E2, 1.0E3, 1.0E4, 1.0E5, 1.0E6, 1.0E7, 1.0E8, 1.0E9, 1.0E10, 1.0E11, 1.0E12, 1.0E13, 1.0E14, 1.0E15, 1.0E16, 1.0E17, 1.0E18, 1.0E19, 1.0E20};
    static constexpr std::array<double, 21> map10iN = {1.0, 1.0E-1, 1.0E-2, 1.0E-3, 1.0E-4, 1.0E-5, 1.0E-6, 1.0E-7, 1.0E-8, 1.0E-9, 1.0E-10, 1.0E-11, 1.0E-12, 1.0E-13, 1.0E-14, 1.0E-15, 1.0E-16, 1.0E-17, 1.0E-18, 1.0E-19, 1.0E-20};
    static constexpr std::array<long long, 10> mapSPGreater1 = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
    static constexpr std::array<double, 11> mapSPLess1 = {1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.00001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001};

    static const double LOG_2_10 = std::log(10) / std::log(2);

    int getFAlpha(int alpha) {
        if (alpha < 0) {
            throw std::invalid_argument("The argument should be greater than 0");
        }
        if (alpha >= f.size()) {
            return static_cast<int>(std::ceil(alpha * LOG_2_10));
        } else {
            return f[alpha];
        }
    }

    double get10iP(int i) {
        if (i < 0) {
            throw std::invalid_argument("The argument should be greater than 0");
        }
        if (i >= map10iP.size()) {
            return std::stod("1.0E" + std::to_string(i));
        } else {
            return map10iP[i];
        }
    }
    
    double get10iN(int i) {
        if (i < 0) {
            throw std::invalid_argument("The argument should be greater than 0");
        }
        if (i >= map10iN.size()) {
            return std::stod("1.0E-" + std::to_string(i));
        } else {
            return map10iN[i];
        }
    }

    double roundUp(double v, int alpha) {
        double scale = get10iP(alpha);
        if (v < 0) {
            return std::floor(v * scale) / scale;
        } else {
            return std::ceil(v * scale) / scale;
        }
    }
    
    // Corresponds to Java's getSPAnd10iNFlag
    std::array<int, 2> getSPAnd10iNFlag(double v) {
        std::array<int, 2> spAnd10iNFlag = {0, 0};
        if (v >= 1) {
            for (size_t i = 0; i < mapSPGreater1.size() - 1; ++i) {
                if (v < mapSPGreater1[i + 1]) {
                    spAnd10iNFlag[0] = i;
                    return spAnd10iNFlag;
                }
            }
        } else {
            for (size_t i = 1; i < mapSPLess1.size(); ++i) {
                if (v >= mapSPLess1[i]) {
                    spAnd10iNFlag[0] = -i;
                    spAnd10iNFlag[1] = (v == mapSPLess1[i]) ? 1 : 0;
                    return spAnd10iNFlag;
                }
            }
        }
        double log10v = std::log10(v);
        spAnd10iNFlag[0] = static_cast<int>(std::floor(log10v));
        spAnd10iNFlag[1] = (log10v == static_cast<long long>(log10v)) ? 1 : 0;
        return spAnd10iNFlag;
    }
    
    // Corresponds to Java's getSignificantCount
    int getSignificantCount(double v, int sp, int lastBetaStar) {
        int i;
        if (lastBetaStar != std::numeric_limits<int>::max() && lastBetaStar != 0) {
            i = std::max(lastBetaStar - sp - 1, 1);
        } else if (lastBetaStar == std::numeric_limits<int>::max()) {
            i = 17 - sp - 1;
        } else if (sp >= 0) {
            i = 1;
        } else {
            i = -sp;
        }

        double temp = v * get10iP(i);
        long long tempLong = static_cast<long long>(temp);
        while (tempLong != temp) {
            i++;
            temp = v * get10iP(i);
            tempLong = static_cast<long long>(temp);
        }

        constexpr double epsilon = 1e-9; // Tolerance for float comparison
        if (std::abs(temp / get10iP(i) - v) > epsilon) {
            return 17;
        } else {
            while (i > 0 && tempLong % 10 == 0) {
                i--;
                tempLong /= 10;
            }
            return sp + i + 1;
        }
    }
    
    // Corresponds to Java's getAlphaAndBetaStar
    std::array<int, 2> getAlphaAndBetaStar(double v, int lastBetaStar) {
        if (v < 0) {
            v = -v;
        }
        std::array<int, 2> alphaAndBetaStar;
        std::array<int, 2> spAnd10iNFlag = getSPAnd10iNFlag(v);
        int beta = getSignificantCount(v, spAnd10iNFlag[0], lastBetaStar);
        alphaAndBetaStar[0] = beta - spAnd10iNFlag[0] - 1;
        alphaAndBetaStar[1] = spAnd10iNFlag[1] == 1 ? 0 : beta;
        return alphaAndBetaStar;
    }

} // namespace Elf64Utils


// =====================================================================================
//                                ElfXORCompressorImpl
// =====================================================================================
class ElfXORCompressorImpl {
public:
    ElfXORCompressorImpl(OutputBitStream& stream) 
        : out(stream), storedLeadingZeros(std::numeric_limits<int>::max()), 
          storedTrailingZeros(std::numeric_limits<int>::max()), storedVal(0), first(true) {}

    OutputBitStream& getOutputStream() { return out; }

    void addValue(long long value) {
        if (first) {
            writeFirst(value);
        } else {
            compressValue(value);
        }
    }

    void close() {
        addValue(double_to_long(std::numeric_limits<double>::quiet_NaN()));
    }

private:
    void writeFirst(long long value) {
        first = false;
        storedVal = value;
        int trailingZeros = count_trailing_zeros(value);
        out.write_bits(trailingZeros, 7);
        if (trailingZeros < 64) {
            out.write_bits(storedVal >> (trailingZeros + 1), 63 - trailingZeros);
        }
    }

    void compressValue(long long value) {
        long long xor_val = storedVal ^ value;

        if (xor_val == 0) {
            out.write_bits(0b01, 2); // case 01
        } else {
            int leadingZeros = count_leading_zeros(xor_val);
            int trailingZeros = count_trailing_zeros(xor_val);
            int roundedLeadingZeros = ElfXORCompressorImpl::leading_round[leadingZeros];

            if (roundedLeadingZeros == storedLeadingZeros && trailingZeros >= storedTrailingZeros) {
                // case 00
                int centerBits = 64 - storedLeadingZeros - storedTrailingZeros;
                out.write_bits(0b00, 2);
                out.write_bits(xor_val >> storedTrailingZeros, centerBits);
            } else {
                storedLeadingZeros = roundedLeadingZeros;
                storedTrailingZeros = trailingZeros;
                int centerBits = 64 - storedLeadingZeros - storedTrailingZeros;

                if (centerBits > 0 && centerBits <= 16) {
                    // case 10
                    uint64_t header = (0b10 << 7) | (ElfXORCompressorImpl::leading_representation[leadingZeros] << 4) | (centerBits & 0xf);
                    out.write_bits(header, 9);
                    out.write_bits(xor_val >> (storedTrailingZeros + 1), centerBits - 1);
                } else {
                    // case 11
                    uint64_t header = (0b11 << 9) | (ElfXORCompressorImpl::leading_representation[leadingZeros] << 6) | (centerBits & 0x3f);
                    out.write_bits(header, 11);
                    if (centerBits > 0)
                        out.write_bits(xor_val >> (storedTrailingZeros + 1), centerBits - 1);
                }
            }
            storedVal = value;
        }
    }

    static constexpr uint8_t leading_representation[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};
    static constexpr uint8_t leading_round[] = {0, 0, 0, 0, 0, 0, 0, 0, 8, 8, 8, 8, 12, 12, 12, 12, 16, 16, 18, 18, 20, 20, 22, 22, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24};

    OutputBitStream& out;
    int storedLeadingZeros;
    int storedTrailingZeros;
    long long storedVal;
    bool first;
};


// =====================================================================================
//                                ElfCompressorImpl
// =====================================================================================
class ElfCompressorImpl {
public:
    ElfCompressorImpl(OutputBitStream& stream) : xorCompressor(stream), lastBetaStar(std::numeric_limits<int>::max()) {}

    void addValue(double v) {
        long long vLong = double_to_long(v);
        long long vPrimeLong;

        if (v == 0.0 || std::isinf(v)) {
            xorCompressor.getOutputStream().write_bits(0b10, 2);
            vPrimeLong = vLong;
        } else if (std::isnan(v)) {
            xorCompressor.getOutputStream().write_bits(0b10, 2);
            vPrimeLong = 0x7ff8000000000000LL;
        } else {
            std::array<int, 2> alphaAndBetaStar = Elf64Utils::getAlphaAndBetaStar(v, lastBetaStar);
            int e = static_cast<int>((vLong >> 52) & 0x7ff);
            int gAlpha = Elf64Utils::getFAlpha(alphaAndBetaStar[0]) + e - 1023;
            
            if (gAlpha < 0) gAlpha = 0;
            int eraseBits = 52 - gAlpha;

            long long mask = (eraseBits >= 64) ? 0 : (0xffffffffffffffffLL << eraseBits);
            long long delta = (~mask) & vLong;

            if (delta != 0 && eraseBits > 4) {
                if (alphaAndBetaStar[1] == lastBetaStar) {
                    xorCompressor.getOutputStream().write_bit(false); // case 0
                } else {
                    uint64_t header = (0b11 << 4) | (alphaAndBetaStar[1] & 0xf);
                    xorCompressor.getOutputStream().write_bits(header, 6); // case 11
                    lastBetaStar = alphaAndBetaStar[1];
                }
                vPrimeLong = mask & vLong;
            } else {
                xorCompressor.getOutputStream().write_bits(0b10, 2); // case 10
                vPrimeLong = vLong;
            }
        }
        xorCompressor.addValue(vPrimeLong);
    }
    
    void close() {
        xorCompressor.getOutputStream().write_bits(0b10, 2);
        xorCompressor.close();
    }

private:
    ElfXORCompressorImpl xorCompressor;
    int lastBetaStar;
};


// =====================================================================================
//                                ElfXORDecompressorImpl
// =====================================================================================
class ElfXORDecompressorImpl {
public:
    ElfXORDecompressorImpl(InputBitStream& stream) 
        : in(stream), storedVal(0), storedLeadingZeros(std::numeric_limits<int>::max()), 
          storedTrailingZeros(std::numeric_limits<int>::max()), first(true), endOfStream(false) {}

    double readValue() {
        next();
        if (endOfStream) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return long_to_double(storedVal);
    }

    InputBitStream& getInputStream() { return in; }

private:
    void next() {
        if (first) {
            first = false;
            int trailingZeros = in.read_bits(7);
            if (trailingZeros < 64) {
                uint64_t bits = in.read_bits(63 - trailingZeros);
                storedVal = ((bits << 1) | 1) << trailingZeros;
            } else {
                storedVal = 0;
            }
            if (storedVal == double_to_long(std::numeric_limits<double>::quiet_NaN())) {
                endOfStream = true;
            }
        } else {
            nextValue();
        }
    }

    void nextValue() {
        int flag = in.read_bits(2);
        switch (flag) {
            case 0b11: { // case 11
                uint64_t leadAndCenter = in.read_bits(9);
                storedLeadingZeros = leading_decode[leadAndCenter >> 6];
                int centerBits = leadAndCenter & 0x3f;
                if (centerBits == 0) centerBits = 64;
                storedTrailingZeros = 64 - storedLeadingZeros - centerBits;
                long long value = 0;
                if (centerBits > 0) {
                    uint64_t bits = in.read_bits(centerBits - 1);
                    value = ((bits << 1) | 1) << storedTrailingZeros;
                }
                storedVal ^= value;
                break;
            }
            case 0b10: { // case 10
                uint64_t leadAndCenter = in.read_bits(7);
                storedLeadingZeros = leading_decode[leadAndCenter >> 4];
                int centerBits = leadAndCenter & 0xf;
                if (centerBits == 0) centerBits = 16;
                storedTrailingZeros = 64 - storedLeadingZeros - centerBits;
                uint64_t bits = in.read_bits(centerBits - 1);
                long long value = ((bits << 1) | 1) << storedTrailingZeros;
                storedVal ^= value;
                break;
            }
            case 0b01: // case 01, same value
                break;
            default: { // case 00
                int centerBits = 64 - storedLeadingZeros - storedTrailingZeros;
                long long value = in.read_bits(centerBits) << storedTrailingZeros;
                storedVal ^= value;
                break;
            }
        }
        if (storedVal == double_to_long(std::numeric_limits<double>::quiet_NaN())) {
            endOfStream = true;
        }
    }

    static constexpr uint8_t leading_decode[] = {0, 8, 12, 16, 18, 20, 22, 24};
    InputBitStream& in;
    long long storedVal;
    int storedLeadingZeros;
    int storedTrailingZeros;
    bool first;
    bool endOfStream;
};


// =====================================================================================
//                                ElfDecompressorImpl (CORRECTED)
// =====================================================================================
class ElfDecompressorImpl {
public:
    ElfDecompressorImpl(InputBitStream& stream) : xorDecompressor(stream), lastBetaStar(std::numeric_limits<int>::max()) {}

    double nextValue() {
        if (in().read_bits(1) == 0) {
            return recoverVByBetaStar();
        } else {
            if (in().read_bits(1) == 0) {
                double val = xorDecompressor.readValue();
                if (std::isnan(val)) return std::numeric_limits<double>::quiet_NaN();
                return val;
            } else {
                lastBetaStar = in().read_bits(4);
                return recoverVByBetaStar();
            }
        }
    }

private:
    double recoverVByBetaStar() {
        double vPrime = xorDecompressor.readValue();
        if (std::isnan(vPrime)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        std::array<int, 2> spAndFlag = Elf64Utils::getSPAnd10iNFlag(std::abs(vPrime));
        int sp = spAndFlag[0];
        double v;
        if (lastBetaStar == 0) {
            v = Elf64Utils::get10iN(-sp - 1); // Note: Java code calls get10iN here
            if (vPrime < 0) {
                v = -v;
            }
        } else {
            int alpha = lastBetaStar - sp - 1;
            v = Elf64Utils::roundUp(vPrime, alpha);
        }
        return v;
    }

    InputBitStream& in() { return xorDecompressor.getInputStream(); }

    ElfXORDecompressorImpl xorDecompressor;
    int lastBetaStar;
};


} // end anonymous namespace

// =====================================================================================
//                                  Public Elf Class
// =====================================================================================

std::vector<uint8_t> Elf::encode(const std::vector<double>& values) {
    if (values.empty()) {
        return {};
    }
    OutputBitStream out;
    ElfCompressorImpl encoder(out);
    for (double val : values) {
        encoder.addValue(val);
    }
    encoder.close();
    return out.get_bytes();
}

std::vector<double> Elf::decode(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) {
        return {};
    }
    InputBitStream in(compressed_data);
    ElfDecompressorImpl decoder(in);
    std::vector<double> values;
    
    while (true) {
        try {
            double val = decoder.nextValue();
            if (std::isnan(val)) {
                break;
            }
            values.push_back(val);
        } catch (const std::out_of_range& e) {
            DEBUG_LOG("ELF Decoding stopped: " << e.what());
            break;
        }
    }
    return values;
}