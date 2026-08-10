#pragma once
#include "algorithms/cluster_encoder_logic.h"
struct ExperimentStats {
 long long initial_k=0,k=0,reference_bits=0,residual_bits=0,paper_reference_bits=0,paper_residual_bits=0;
 long long before_bits=0,after_bits=0,reassigned=0,deleted=0,queries=0,deletion_queries=0,distance_evaluations=0;
 long long workspace_bytes=0,checked=0,skipped=0;
};
extern int experimentBudget;
extern bool experimentVerify;
extern ExperimentStats experimentStats;
ClusteringResult experimentSelect(const SoAData<long long>& data);
void experimentSelfTest();

namespace ClusterLogic {
struct ExperimentEncodedPage { std::vector<uint8_t> bytes; size_t valid_bits; };
ExperimentEncodedPage encodeExperimentPage(const std::vector<std::string>& lines, int& dim);
ClusteringResult sortClusteringBySize(std::vector<std::vector<long long>>& refs,
    std::vector<int>& assignments, std::vector<long long>& sizes);
}
