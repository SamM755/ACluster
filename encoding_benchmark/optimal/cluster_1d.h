#pragma once
#include "cluster_encoder_logic.h"

namespace ClusterLogic {

KMedoidResult kMedoid1D(const SoAData<long long>& data, int k, int max_iter,
    double tolerance, uint32_t seed);
AClusterResult adaptiveGreedy1D(const SoAData<long long>& data);
}
