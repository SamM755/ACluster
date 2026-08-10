#include "optimal_1d.h"
#include "cluster_encoder_logic.h"
#include <algorithm>
#include <chrono>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace ClusterLogic {
int paperScalarCost(long long x, long long b) {

    const uint64_t distance = x >= b ? uint64_t(x) - uint64_t(b) : uint64_t(b) - uint64_t(x);
    return distance == 0 ? 1 : 65 - __builtin_clzll(distance);
}

Optimal1DResult optimal1D(const SoAData<long long>& data, int k) {
    using Clock = std::chrono::steady_clock;
    auto ns = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::nanoseconds>(b-a).count(); };
    const auto start = Clock::now();
    if (data.dim != 1 || data.columns.size() != 1 || data.num_rows < 0 ||
        data.num_rows > 10000 || data.columns[0].size() != static_cast<size_t>(data.num_rows) || k < 0)
        throw std::invalid_argument("Optimal requires valid 1D input of at most 10000 rows.");
    Optimal1DResult result;
    const int n = data.num_rows;
    if (!n) return result;
    const bool adaptive = k == 0;
    std::vector<int> order(n);
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b) { return data.columns[0][a] < data.columns[0][b]; });
    std::vector<long long> values;
    std::vector<int> counts, starts;
    for (int i=0;i<n;++i) {
        const auto value = data.columns[0][order[i]];
        if (values.empty() || values.back() != value) { values.push_back(value); counts.push_back(0); starts.push_back(i); }
        ++counts.back();
    }
    starts.push_back(n);
    const int u = static_cast<int>(values.size());
    result.dp_states = result.distinct_values = u;
    k = std::min(k,u);
    const auto sorted_at = Clock::now();
    result.sort_ns = ns(start,sorted_at);
    result.result.cluster_assignment.resize(n);
    auto append_cluster = [&](int l,int end,int ref) {
        const int id = static_cast<int>(result.result.medoids.size());
        result.result.medoids.push_back({values[ref]});
        result.result.cluster_sizes.push_back(starts[end]-starts[l]);
        result.reference_cost += paperScalarCost(values[ref]);
        for (int i=l;i<end;++i) result.residual_cost += static_cast<long long>(counts[i])*paperScalarCost(values[i],values[ref]);
        for (int i=starts[l];i<starts[end];++i) result.result.cluster_assignment[order[i]] = id;
    };
    if (!adaptive && k == u) {

        for (int i=0;i<u;++i) append_cluster(i,i+1,i);
        result.objective = n;
        result.result = sortClusteringBySize(result.result.medoids,result.result.cluster_assignment,result.result.cluster_sizes);
        result.dp_ns = ns(sorted_at,Clock::now());
        return result;
    }
    std::vector<size_t> offsets(u+1,0);
    for (int l=0;l<u;++l) offsets[l+1] = offsets[l]+u-l;
    constexpr uint32_t inf = std::numeric_limits<uint32_t>::max()/4;
    std::vector<uint32_t> costs(offsets[u],inf), prefix;
    std::vector<uint16_t> references;

        prefix.resize(static_cast<size_t>(u)*(u+1),0);
        for (int r=0;r<u;++r) {
            auto* row = prefix.data()+static_cast<size_t>(r)*(u+1);
            for (int i=0;i<u;++i) row[i+1] = row[i]+counts[i]*paperScalarCost(values[i],values[r]);
        }
        const auto prepared_at = Clock::now();
        result.prefix_ns = ns(sorted_at,prepared_at);

            references.resize(offsets[u]);
            for (int l=0;l<u;++l) {
                costs[offsets[l]] = counts[l]+(adaptive ? paperScalarCost(values[l]) : 0);
                references[offsets[l]] = static_cast<uint16_t>(l);
            }
            result.reference_candidates = u;

            for (int length=2;length<=u;++length) {
                for (int l=0;l+length<=u;++l) {
                    const int end = l+length;
                    const auto index = offsets[l]+length-1;
                    const int lower = references[index-1];
                    const int upper = references[offsets[l+1]+length-2];
                    if (lower > upper) throw std::runtime_error("Invalid monotone reference bounds.");
                    uint32_t best = inf;
                    int chosen = lower;
                    for (int r=lower;r<=upper;++r) {
                        const auto* row = prefix.data()+static_cast<size_t>(r)*(u+1);
                        const auto candidate = row[end]-row[l]+(adaptive ? paperScalarCost(values[r]) : 0);
                        if (candidate < best) { best = candidate; chosen = r; }
                    }
                    result.reference_candidates += upper-lower+1;
                    costs[index] = best;
                    references[index] = static_cast<uint16_t>(chosen);
                }
            }

    result.search_ns = ns(prepared_at,Clock::now());
    const auto searched_at = Clock::now();
    std::vector<uint32_t> previous(u+1,inf),current(u+1,inf);
    std::vector<int> splits(static_cast<size_t>(adaptive ? 1 : k+1)*(u+1),-1);
    previous[0] = 0;
    if (adaptive) {

        for (int l=0;l<u;++l) {
            const auto base = previous[l]; const auto* row = costs.data()+offsets[l];
            for (int end=l+1;end<=u;++end) {
                const auto candidate = base+row[end-l-1];
                if (candidate < previous[end]) { previous[end] = candidate; splits[end] = l; }
            }
        }
    } else {
        for (int groups=1;groups<=k;++groups) {
            std::fill(current.begin(),current.end(),inf);
            auto* split_row = splits.data()+static_cast<size_t>(groups)*(u+1);
            for (int l=groups-1;l<u;++l) {
                const auto base = previous[l];
                if (base == inf) continue;
                const auto* row = costs.data()+offsets[l];
                for (int end=l+1;end<=u;++end) {
                    const auto candidate = base+row[end-l-1];
                    if (candidate < current[end]) { current[end] = candidate; split_row[end] = l; }
                }
            }
            previous.swap(current);
        }
    }
    result.objective = previous[u];
    std::vector<std::pair<int,int>> intervals;
    for (int end=u,groups=k;end>0;--groups) {
        const int l = splits[(adaptive ? 0 : static_cast<size_t>(groups)*(u+1))+end];
        if (l < 0 || l >= end) throw std::runtime_error("Optimal backtracking failed.");
        intervals.emplace_back(l,end); end = l;
    }
    std::reverse(intervals.begin(),intervals.end());
    for (const auto& interval : intervals) {
        const int l = interval.first,end = interval.second;
        const int ref = references[offsets[l]+end-l-1];
        append_cluster(l,end,ref);
    }
    if (result.objective != result.residual_cost+(adaptive ? result.reference_cost : 0))
        throw std::runtime_error("Optimal reconstructed cost mismatch.");
    result.result = sortClusteringBySize(result.result.medoids,result.result.cluster_assignment,result.result.cluster_sizes);
    result.dp_ns = ns(searched_at,Clock::now());
    return result;
}
}
