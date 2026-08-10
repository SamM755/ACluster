#include "pipeline.h"
#include "optimal_1d.h"
namespace {
long long cost(long long a, long long b=0) {
    uint64_t d=a>=b?uint64_t(a)-uint64_t(b):uint64_t(b)-uint64_t(a);
    long long n=1;
    while(d){++n;d>>=1;}
    return n;
}
long long oracle(std::vector<long long> x, int k) {
    std::sort(x.begin(),x.end());
    const int n=x.size();
    const long long inf=1LL<<50;
    std::vector<std::vector<long long>> dp(n+1,std::vector<long long>(n+1,inf));
    dp[0][0]=0;
    for(int g=1;g<=n;++g) for(int end=g;end<=n;++end) for(int start=g-1;start<end;++start){
        long long best=inf;
        for(int r=start;r<end;++r){
            long long v=k?0:cost(x[r]);
            for(int i=start;i<end;++i)v+=cost(x[i],x[r]);
            best=std::min(best,v);
        }
        dp[g][end]=std::min(dp[g][end],dp[g-1][start]+best);
    }
    long long best=inf;
    for(int g=1;g<=n;++g) best=std::min(best,dp[g][n]);
    return k?dp[k][n]:best;
}
}
void test_optimal_1d() {
    std::mt19937 rng(20261001);
    for(int trial=0;trial<300;++trial){
        const int n=1+rng()%9;
        SoAData<long long> data;
        data.dim=1;data.num_rows=n;data.columns.resize(1);
        for(int i=0;i<n;++i)data.columns[0].push_back(int(rng()%33)-16);
        const int distinct=std::set<long long>(data.columns[0].begin(),data.columns[0].end()).size();
        for(int k=0;k<=distinct;++k){
            auto actual=ClusterLogic::optimal1D(data,k);
            if(actual.objective!=oracle(data.columns[0],k))throw std::runtime_error("Optimal oracle mismatch");
            long long residual=0,reference=0;
            for(const auto& r:actual.result.medoids)reference+=cost(r[0]);
            for(int i=0;i<n;++i)residual+=cost(data.columns[0][i],actual.result.medoids.at(actual.result.cluster_assignment.at(i))[0]);
            if(actual.objective!=residual+(k?0:reference))throw std::runtime_error("Optimal reconstruction mismatch");
        }
    }
    std::cout << "Optimal exhaustive oracle checks passed\n";
}
