#pragma once
#include "duckdb.hpp"
namespace duckdb {
class ClusterExtension : public Extension {
public:
	void Load(ExtensionLoader &loader) override;
	std::string Name() override { return "cluster"; }
	std::string Version() const override { return "0.1.0"; }
};
}
