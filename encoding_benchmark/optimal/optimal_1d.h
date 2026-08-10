#pragma once
#include "cluster_common.h"

namespace ClusterLogic {
struct Optimal1DResult {
    ClusteringResult result;
    long long residual_cost = 0;
    long long reference_cost = 0;
    long long objective = 0;
    int distinct_values = 0;
    int dp_states = 0;
    uint64_t reference_candidates = 0;
    long long sort_ns = 0, prefix_ns = 0, search_ns = 0, dp_ns = 0;
};

int paperScalarCost(long long x, long long b = 0);

Optimal1DResult optimal1D(const SoAData<long long>& data, int k = 0);
}
