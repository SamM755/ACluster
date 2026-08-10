#ifndef VORTEX_H
#define VORTEX_H

#include <cstdint>
#include <string>
#include <vector>

class Vortex {
public:
    std::string get_name() const { return "Vortex"; }

    std::vector<uint32_t> compute_permutation(const std::vector<std::vector<double>>& rows) const;

    std::vector<std::vector<double>> reorder_rows(const std::vector<std::vector<double>>& rows) const;
};

#endif
