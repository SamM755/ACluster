#include "cluster_1d.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <unordered_map>
#include <unordered_set>

namespace ClusterLogic {
namespace {

int cost(long long x, long long b=0) {
    const uint64_t d = x>=b ? uint64_t(x)-uint64_t(b) : uint64_t(b)-uint64_t(x);
    return d ? 65-__builtin_clzll(d) : 2;
}
struct Groups {
    std::vector<long long> values;
    std::vector<int> count, row;
    explicit Groups(const SoAData<long long>& data) {
        if (data.dim!=1 || data.columns.size()!=1 || data.num_rows<0 ||
            data.columns[0].size()!=static_cast<size_t>(data.num_rows))
            throw std::invalid_argument("Scalar heuristics require valid 1D input.");
        std::unordered_map<long long,int> index;
        row.reserve(data.num_rows);
        for (const auto x : data.columns[0]) {
            const auto inserted = index.emplace(x,static_cast<int>(values.size()));
            const int id = inserted.first->second;
            if (inserted.second) { values.push_back(x); count.push_back(0); }
            ++count[id]; row.push_back(id);
        }
    }
};
ClusteringResult finish(const Groups& g, const std::vector<long long>& references,
    std::vector<int>& assignments) {
    std::vector<std::vector<long long>> medoids;
    for (auto x : references) medoids.push_back({x});
    std::vector<long long> sizes(references.size(),0);
    for (int owner : assignments) ++sizes.at(owner);
    return sortClusteringBySize(medoids,assignments,sizes);
}

std::vector<long long> initialize(const Groups& g, int k, uint32_t seed) {
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> first(0,static_cast<int>(g.row.size())-1);
    std::vector<long long> medoids{g.values[g.row[first(gen)]]};
    std::unordered_set<long long> selected{medoids.front()};
    std::vector<long long> distance(g.values.size());
    for (size_t i=0;i<distance.size();++i) distance[i]=std::abs(g.values[i]-medoids[0]);
    for (int m=1;m<k;++m) {
        long long total=0;
        for (size_t i=0;i<distance.size();++i) total+=distance[i]*g.count[i];
        int chosen=0;
        if (total) {
            std::uniform_real_distribution<double> draw(0.0,static_cast<double>(total));
            const long long target=static_cast<long long>(draw(gen));
            long long prefix=0;
            for (;chosen<static_cast<int>(g.row.size());++chosen) {
                prefix+=distance[g.row[chosen]];
                if (prefix>=target) break;
            }
            if (chosen==static_cast<int>(g.row.size())) --chosen;
        }
        while (selected.count(g.values[g.row[chosen]])) chosen=(chosen+1)%g.row.size();
        const auto next=g.values[g.row[chosen]];
        medoids.push_back(next); selected.insert(next);
        if (total) for (size_t i=0;i<distance.size();++i)
            distance[i]=std::min(distance[i],std::abs(g.values[i]-next));
    }
    return medoids;
}
}

KMedoidResult kMedoid1D(const SoAData<long long>& data, int k, int max_iter,
    double tolerance, uint32_t seed) {
    Groups g(data);
    if (g.row.empty()) return {{},0};
    if (k<1 || max_iter<1) throw std::invalid_argument("Invalid scalar KCluster parameters.");
    k=std::min(k,static_cast<int>(g.values.size()));
    auto medoids=initialize(g,k,seed);
    std::vector<int> owner(g.values.size()),assignments(g.row.size());
    long long previous=std::numeric_limits<long long>::max(),final_cost=0;
    for (int iteration=0;iteration<max_iter;++iteration) {
        long long total=0;
        for (size_t v=0;v<g.values.size();++v) {
            int best=std::numeric_limits<int>::max(),leader=-1;
            for (int m=0;m<k;++m) {
                const int c=cost(g.values[v],medoids[m]);
                if (c<best) { best=c; leader=m; }
                else if (c==best && g.values[v]==medoids[m]) { leader=m; break; }
            }
            owner[v]=leader; total+=static_cast<long long>(best)*g.count[v];
        }
        final_cost=total;
        if (iteration>0 && std::abs(static_cast<double>(previous)-total)<tolerance) break;
        previous=total;
        std::vector<std::vector<int>> members(k);

        for (size_t v=0;v<g.values.size();++v) members[owner[v]].push_back(static_cast<int>(v));
        auto next=medoids;
        for (int m=0;m<k;++m) {
            long long best=std::numeric_limits<long long>::max();
            for (int candidate : members[m]) {
                long long sum=-2;
                for (int v : members[m]) sum+=static_cast<long long>(g.count[v])*cost(g.values[candidate],g.values[v]);
                if (sum<best) { best=sum; next[m]=g.values[candidate]; }
            }
        }
        if (next==medoids) break;
        medoids=std::move(next);
    }
    for (size_t i=0;i<g.row.size();++i) assignments[i]=owner[g.row[i]];
    return {finish(g,medoids,assignments),final_cost};
}

AClusterResult adaptiveGreedy1D(const SoAData<long long>& data) {
    Groups g(data);
    const int n=static_cast<int>(g.row.size()),u=static_cast<int>(g.values.size());
    if (!n) return {{},0};

    std::vector<long long> medoids{g.values[g.row[0]]};
    std::vector<std::vector<int>> members{{0}};
    std::vector<int> existing(u,-1),checked(u,0),nearest(u,-1),best_cost(u,66);
    std::vector<int> assignments(n,-1),row_cost(n,2);
    existing[g.row[0]]=0; assignments[0]=0;
    std::vector<std::pair<int,int>> moves;
    for (int i=1;i<n;++i) {
        const int v=g.row[i]; const auto x=g.values[v];
        if (existing[v]>=0) {
            const int leader=existing[v]; assignments[i]=leader; members[leader].push_back(i);
            continue;
        }

        for (int m=checked[v];m<static_cast<int>(medoids.size());++m) {
            const int c=cost(x,medoids[m]);
            if (c<best_cost[v]) { best_cost[v]=c; nearest[v]=m; }
        }
        checked[v]=static_cast<int>(medoids.size());
        const int leader=nearest[v];
        long long savings=best_cost[v];
        moves.clear();
        for (int row : members[leader]) {
            const int c=cost(g.values[g.row[row]],x);
            if (c<row_cost[row]) { savings+=row_cost[row]-c; moves.emplace_back(row,c); }
        }
        if (savings>cost(x)) {
            const int next=static_cast<int>(medoids.size());
            medoids.push_back(x); existing[v]=next;
            assignments[i]=next; members.push_back({i});
            auto& old=members[leader];
            size_t move=0,keep=0;
            for (int row : old) {
                if (move<moves.size() && moves[move].first==row) {
                    assignments[row]=next; row_cost[row]=moves[move].second;
                    members[next].push_back(row); ++move;
                } else old[keep++]=row;
            }
            old.resize(keep);
        } else {
            assignments[i]=leader; row_cost[i]=best_cost[v]; members[leader].push_back(i);
        }
    }
    long long total=0;
    for (int c : row_cost) total+=c;
    for (auto reference : medoids) total+=cost(reference);
    return {finish(g,medoids,assignments),total};
}
}
