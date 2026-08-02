#include "experiment_api.h"

#include <filesystem>

std::vector<DatasetConfig> public_datasets() {
    return {
        {"SSD", "SSD.csv", 1, 8926},
        {"BC", "BC.csv", 1, 99999},
        {"PI", "PI.csv", 1, 14825},
        {"WS", "WS.csv", 1, 100000},
        {"CT", "CT.csv", 3, 581012},
        {"Crop", "Crop.csv", 20, 24000},
        {"Gas", "Gas.csv", 48, 13910},
        {"BT", "BT.csv", 49, 583250},
        {"Musk", "Musk.csv", 111, 6598},
    };
}

std::string dataset_path(const ExperimentContext& context, const DatasetConfig& dataset) {
    return (std::filesystem::path(context.data_root) / dataset.file_name).string();
}

std::string experiment_output(const ExperimentContext& context, const std::string& name) {
    return (std::filesystem::path(context.output_root) / name).string();
}
