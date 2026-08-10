#include "experiment_hook.h"
#include <algorithm>
#include <limits>
#include <numeric>
#include <random>
#include <iostream>
using I=long long;
using U=unsigned long long;
using Wide=unsigned __int128;
int experimentBudget=-2;
bool experimentVerify=false;
ExperimentStats experimentStats;
namespace {
U dist(I a,I b){return a>=b?U(a)-U(b):U(b)-U(a);}
int bits(U x){return x?65-__builtin_clzll(x):2;}
I rcost(const std::vector<I>& v){I c=0;for(I x:v)c+=bits(dist(x,0));return c;}
I cost(const SoAData<I>& data,int row,const std::vector<I>& ref){
 I c=0;for(int d=0;d<data.dim;++d)c+=bits(dist(data.columns[d][row],ref[d]));return c;
}
struct Bound {I cost=0;Wide l1=0;};
bool lessBound(const Bound&a,const Bound&b){return a.cost!=b.cost?a.cost<b.cost:a.l1<b.l1;}
struct Node {int lo=0,hi=0,left=-1,right=-1,parent=-1,live=0;};
class ReferenceIndex {
 const std::vector<std::vector<I>>& refs;
 const std::vector<unsigned char>& live;
 int dim;
 std::vector<int> ids,leaf;
 std::vector<Node> nodes;
 std::vector<I> low,high;
 int build(int begin,int end,int parent){
  int id=nodes.size();nodes.push_back({begin,end,-1,-1,parent,end-begin});
  low.resize(size_t(id+1)*dim);high.resize(size_t(id+1)*dim);
  int axis=0;U range=0;
  for(int d=0;d<dim;++d){
   I lo=refs[ids[begin]][d],hi=lo;
   for(int i=begin+1;i<end;++i){lo=std::min(lo,refs[ids[i]][d]);hi=std::max(hi,refs[ids[i]][d]);}
   low[size_t(id)*dim+d]=lo;high[size_t(id)*dim+d]=hi;
   if(dist(lo,hi)>range){range=dist(lo,hi);axis=d;}
  }
  if(end-begin<=8){for(int i=begin;i<end;++i)leaf[ids[i]]=id;return id;}
  int mid=(begin+end)/2;
  std::nth_element(ids.begin()+begin,ids.begin()+mid,ids.begin()+end,[&](int a,int b){
   return refs[a][axis]!=refs[b][axis]?refs[a][axis]<refs[b][axis]:refs[a]<refs[b];
  });
  int a=build(begin,mid,id),b=build(mid,end,id);nodes[id].left=a;nodes[id].right=b;return id;
 }
 Bound bound(int node,const SoAData<I>& data,int row)const {
  Bound b;
  for(int d=0;d<dim;++d){
   I x=data.columns[d][row],lo=low[size_t(node)*dim+d],hi=high[size_t(node)*dim+d];
   U delta=x<lo?dist(x,lo):x>hi?dist(x,hi):0;b.cost+=bits(delta);b.l1+=delta;
  }
  return b;
 }
 void consider(int r,const SoAData<I>& data,int row,int exclude,int& winner,Bound& best,bool count=true){
  if(r<0||!live[r]||r==exclude)return;
  if(count)++experimentStats.distance_evaluations;
  Bound b;
  for(int d=0;d<dim;++d){
   U delta=dist(data.columns[d][row],refs[r][d]);b.cost+=bits(delta);b.l1+=delta;
   if(b.cost>best.cost)return;
  }
  if(winner<0||lessBound(b,best)||(!lessBound(best,b)&&refs[r]<refs[winner])){best=b;winner=r;}
 }
 void visit(int id,const SoAData<I>& data,int row,int exclude,int& winner,Bound& best){
  if(nodes[id].live==0)return;
  auto bd=bound(id,data,row);if(lessBound(best,bd))return;
  const auto node=nodes[id];
  if(node.left<0){for(int i=node.lo;i<node.hi;++i)consider(ids[i],data,row,exclude,winner,best);return;}
  auto a=bound(node.left,data,row),b=bound(node.right,data,row);
  int first=node.left,second=node.right;if(lessBound(b,a))std::swap(first,second);
  visit(first,data,row,exclude,winner,best);visit(second,data,row,exclude,winner,best);
 }
public:
 ReferenceIndex(const std::vector<std::vector<I>>& r,const std::vector<unsigned char>& l):
  refs(r),live(l),dim(r[0].size()),ids(r.size()),leaf(r.size()){
  std::iota(ids.begin(),ids.end(),0);nodes.reserve((r.size()+7)/8*4);
  low.reserve(nodes.capacity()*dim);high.reserve(nodes.capacity()*dim);build(0,r.size(),-1);
 }
 int nearest(const SoAData<I>& data,int row,int exclude,int initial){
  ++experimentStats.queries;int winner=-1;
  Bound best{std::numeric_limits<I>::max(),~Wide(0)};
  consider(initial,data,row,exclude,winner,best);visit(0,data,row,exclude,winner,best);
  if(winner<0)throw std::runtime_error("No alternative reference");
  if(experimentVerify){
   int direct=-1;Bound best2{std::numeric_limits<I>::max(),~Wide(0)};
   for(int r=0;r<int(refs.size());++r)consider(r,data,row,exclude,direct,best2,false);
   if(direct!=winner)throw std::runtime_error("Multidimensional index differs from exhaustive bit-cost search");
  }
  return winner;
 }
 void remove(int r){for(int id=leaf[r];id>=0;id=nodes[id].parent)--nodes[id].live;}
 size_t bytes()const{return sizeof(*this)+(ids.capacity()+leaf.capacity())*sizeof(int)+
  nodes.capacity()*sizeof(Node)+(low.capacity()+high.capacity())*sizeof(I);}
};
}
ClusteringResult experimentSelect(const SoAData<I>& data){
 experimentStats={};auto& st=experimentStats;
 auto baseline=ClusterLogic::adaptiveGreedyBasePointSelection_with_diag(data,data.dim);
 ClusteringResult result=std::move(baseline.result);
 int n=data.num_rows,k=result.medoids.size();st.initial_k=k;st.before_bits=baseline.final_cost;
 if(n==0)return result;
 if(experimentBudget>=0){
  auto& refs=result.medoids;auto& owner=result.cluster_assignment;
  std::vector<int> head(k,-1),next(n,-1),order(k),compact(k,-1);
  std::vector<unsigned char> live(k,1);
  std::vector<std::pair<int,int>> alternatives;alternatives.reserve(n);
  ReferenceIndex index(refs,live);
  auto workspace=[&](){return index.bytes()+(head.capacity()+next.capacity()+order.capacity()+compact.capacity())*sizeof(int)
    +live.capacity()+alternatives.capacity()*sizeof(std::pair<int,int>);};
  st.workspace_bytes=workspace();
  std::fill(result.cluster_sizes.begin(),result.cluster_sizes.end(),0);
  auto total=[&](){I value=0;for(int r=0;r<k;++r)if(live[r])value+=rcost(refs[r]);
   for(int i=0;i<n;++i)value+=cost(data,i,refs[owner[i]]);return value;};
  if(experimentVerify&&total()!=st.before_bits)throw std::runtime_error("Baseline cost mismatch");
  for(int row=0;row<n;++row){
   int to=index.nearest(data,row,-1,owner[row]);
   if(cost(data,row,refs[to])>cost(data,row,refs[owner[row]]))throw std::runtime_error("Reassignment increases cost");
   st.reassigned+=owner[row]!=to;owner[row]=to;next[row]=head[to];head[to]=row;++result.cluster_sizes[to];
  }
  std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),[&](int a,int b){
   return result.cluster_sizes[a]!=result.cluster_sizes[b]?result.cluster_sizes[a]<result.cluster_sizes[b]:refs[a]<refs[b];
  });
  int remaining=k;I budget=I(experimentBudget)*n;
  for(int r:order){
   if(remaining==1)break;
   if(st.deletion_queries+result.cluster_sizes[r]>budget){++st.skipped;continue;}
   ++st.checked;I gain=rcost(refs[r]);alternatives.clear();
   for(int row=head[r];row>=0;row=next[row]){
    int to=index.nearest(data,row,r,-1);++st.deletion_queries;
    alternatives.emplace_back(row,to);gain-=cost(data,row,refs[to])-cost(data,row,refs[r]);
   }
   if(gain<=0)continue;
   I before=experimentVerify?total():0;live[r]=0;index.remove(r);
   for(auto item:alternatives){
    int row=item.first,to=item.second;owner[row]=to;next[row]=head[to];head[to]=row;++result.cluster_sizes[to];
   }
   head[r]=-1;result.cluster_sizes[r]=0;--remaining;++st.deleted;
   if(experimentVerify&&before-total()!=gain)throw std::runtime_error("Deletion audit failed");
  }
  std::vector<std::vector<I>> kept;std::vector<I> sizes;
  kept.reserve(remaining);sizes.reserve(remaining);
  for(int r=0;r<k;++r)if(live[r]){
   compact[r]=kept.size();kept.push_back(refs[r]);sizes.push_back(result.cluster_sizes[r]);
  }
  I finalCapacity=kept.capacity()*sizeof(std::vector<I>)+sizes.capacity()*sizeof(I);
  for(const auto& ref:kept)finalCapacity+=ref.capacity()*sizeof(I);
  st.workspace_bytes=workspace()+finalCapacity;
  for(int& r:owner){if(!live[r])throw std::runtime_error("Deleted assignment");r=compact[r];}
  result=ClusterLogic::sortClusteringBySize(kept,owner,sizes);
  if(st.deletion_queries>budget)throw std::runtime_error("Query budget exceeded");
 }
 st.k=result.medoids.size();
 for(const auto& ref:result.medoids)for(I x:ref){
  st.reference_bits+=bits(dist(x,0));st.paper_reference_bits+=x?bits(dist(x,0)):1;
 }
 for(int row=0;row<n;++row)for(int d=0;d<data.dim;++d){
  U delta=dist(data.columns[d][row],result.medoids[result.cluster_assignment[row]][d]);
  st.residual_bits+=bits(delta);st.paper_residual_bits+=delta?bits(delta):1;
 }
 st.after_bits=st.reference_bits+st.residual_bits;
 if(st.after_bits>st.before_bits||st.k+st.deleted!=st.initial_k)throw std::runtime_error("Cost or reference accounting failed");
 return result;
}
void experimentSelfTest(){
 std::mt19937 rng(20261005);experimentVerify=true;long long deleted=0;
 for(int test=0;test<240;++test){
  int n=1+test%79,dim=std::vector<int>{1,2,4,23,49}[test%5];
  SoAData<I> data;data.dim=dim;data.num_rows=n;data.columns.assign(dim,std::vector<I>(n));
  for(int i=0;i<n;++i)for(int d=0;d<dim;++d)data.columns[d][i]=I(rng()%(test%2?200:7))-3;
  auto expected=ClusterLogic::adaptiveGreedyBasePointSelection_with_diag(data,dim);
  experimentBudget=-2;auto base=experimentSelect(data);
  if(base.medoids!=expected.result.medoids||base.cluster_assignment!=expected.result.cluster_assignment||
   base.cluster_sizes!=expected.result.cluster_sizes||experimentStats.after_bits!=expected.final_cost)
    throw std::runtime_error("Baseline decisions differ");
  for(int b:{0,1,4}){
   experimentBudget=b;auto out=experimentSelect(data);deleted+=experimentStats.deleted;
   if(std::accumulate(out.cluster_sizes.begin(),out.cluster_sizes.end(),I(0))!=n)throw std::runtime_error("Lost rows");
  }
 }
 experimentVerify=false;std::cout<<"PASS 240 baseline cases and 720 multidimensional cleanup cases; "<<deleted<<" deletions audited\n";
}
