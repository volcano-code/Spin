#include "spin2014/SourceAdaptiveDomain.h"
#include "source_surface_test_utils.h"
#include "manufacturing/SlicerValidation.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <type_traits>

using namespace spin2014;
using namespace surface_test;

// Allows an executable red test against the pre-interface production header.
template<class T,class=void> struct HasRestore:std::false_type {};
template<class T> struct HasRestore<T,std::void_t<decltype(&T::restorePartition)>>:std::true_type {};
std::shared_ptr<const SourceSurfaceOctree> source(double outerHi=1.) {
    auto outer=cube(0,outerHi),design=cube(.125,.875);
    SourceSurfaceOctreeConfig c;c.maximumDepth=4;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
}
std::vector<SourcePartitionKey> keys(const SourceAdaptiveDomain& d) {
    std::vector<SourcePartitionKey> result;for(const auto& l:d.leaves())result.push_back(l);return result;
}
// Deliberately forge a checksum consistent with an invalid input. Integrity
// hashes cannot replace coverage validation or authenticate a checkpoint.
std::string forgedHash(const SourceSurfaceOctree& s,const std::vector<SourcePartitionKey>& k,std::size_t budget=10000) {
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<std::hexfloat<<"source-adaptive-partition-v1 "<<s.partitionSha256()<<" 2 4 "<<budget<<' ';
    for(const auto& l:k) {
        out<<l.level<<':'<<l.coordinates[0]<<','<<l.coordinates[1]<<','<<l.coordinates[2]<<';';
        const auto moments=s.cell(l.level,l.coordinates).moments.values();
        for(int j=0;j<10;++j)out<<moments[j]<<',';
    }
    return manufacturing::computeContentSha256(out.str());
}
void same(const SourceAdaptiveDomain& a,const SourceAdaptiveDomain& b,const Eigen::VectorXd& beta) {
    require(a.partitionHash()==b.partitionHash() && keys(a)==keys(b),"partition identity/order changed");
    require(a.edges()==b.edges(),"restored adjacency changed");
    require(Eigen::MatrixXd(a.uniformLaplacian())==Eigen::MatrixXd(b.uniformLaplacian()),"restored Laplacian changed");
    require(a.removedMoments(beta).values()==b.removedMoments(beta).values(),"restored material moments changed");
    for(std::size_t i=0;i<a.leaves().size();++i)
        require(a.leaves()[i].moments.values()==b.leaves()[i].moments.values(),"per-cell moments changed");
}
template<class Domain> void restorationTests() {
    if constexpr(!HasRestore<Domain>::value)throw std::runtime_error("partition restoration API missing");
    else {
        auto s=source();Domain original(s,2,4,10000);
        auto restore=[&](const Domain& d){return Domain::restorePartition(s,2,4,10000,keys(d),d.partitionHash());};
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(original.leaves().size(),.37);
        auto copy=restore(original);same(original,copy,beta);
        // Independent analytic conservative grid: closed boundary contacts are
        // excluded, leaving [3/16,13/16]^3 at source geometry depth four.
        close(copy.removedMoments(Eigen::VectorXd::Ones(beta.size())).values(),
              boxMoments({3./16,3./16,3./16},{13./16,13./16,13./16}),"analytic restored volume/moments");
        for(Eigen::Index i=0;i<beta.size();++i)beta[i]=i%3?0.:.37;
        original.splitFractionalLeaves(beta,.05);
        copy=restore(original);same(original,copy,beta);
        // Continue the actual split/merge flow from the restored mixed partition.
        auto continuedBeta=beta;
        const auto left=original.splitFractionalLeaves(beta,.05);
        const auto right=copy.splitFractionalLeaves(continuedBeta,.05);
        require(left.afterHash==right.afterHash && beta==continuedBeta,"post-restore split diverged");
        same(original,copy,beta);
        beta.setZero();continuedBeta.setZero();
        original.mergeBinarySiblings(beta,.05);copy.mergeBinarySiblings(continuedBeta,.05);
        same(original,copy,beta);
        require(original.leaves().size()==1 && original.leaves()[0].level==0,"recursive merge not exercised");
        copy=restore(original);same(original,copy,beta);
        // Hash binds source surface identity even when the carveable region is
        // unchanged. It also binds original initial level and resource policy.
        const auto good=keys(original);const auto hash=original.partitionHash();
        auto reject=[&](const std::vector<SourcePartitionKey>& k,const std::string& h){rejects([&]{Domain::restorePartition(s,2,4,10000,k,h);},"invalid restored partition accepted");};
        require(forgedHash(*s,good)==hash,"test hash encoder disagrees with existing format");
        reject({},forgedHash(*s,{}));reject(good,std::string(64,'0'));reject(good,"");
        auto bad=good;bad.push_back(good.front());reject(bad,hash);
        bad=good;bad.push_back({1,{0,0,0}});reject(bad,forgedHash(*s,bad)); // parent + descendant
        bad=good;bad.front().level=-1;reject(bad,hash);
        bad=good;bad.front().level=5;reject(bad,hash);
        bad=good;bad.front().coordinates[0]=UINT32_MAX;reject(bad,hash);
        rejects([&]{Domain::restorePartition(source(.95),2,4,10000,good,hash);},"wrong source accepted");
        rejects([&]{Domain::restorePartition(s,1,4,10000,good,hash);},"changed initial level accepted");
        rejects([&]{Domain::restorePartition(s,2,3,10000,good,hash);},"changed maximum level accepted");
        rejects([&]{Domain::restorePartition(s,2,4,10001,good,hash);},"changed resource identity accepted");
        rejects([&]{Domain::restorePartition(nullptr,2,4,10000,good,hash);},"null source accepted");
        rejects([&]{Domain::restorePartition(s,2,4,0,good,hash);},"zero budget accepted");
        Domain many(s,2,4,10000);auto manyKeys=keys(many);
        rejects([&]{Domain::restorePartition(s,2,4,1,manyKeys,forgedHash(*s,manyKeys,1));},"restore leaf budget ignored");
        bad=manyKeys;bad.erase(bad.begin());reject(bad,forgedHash(*s,bad));
        bad=manyKeys;std::reverse(bad.begin(),bad.end());reject(bad,many.partitionHash());
        // One shifted leaf lies inside the source root but outside carving.
        bad=manyKeys;bad.front()={4,{0,0,0}};std::sort(bad.begin(),bad.end());reject(bad,forgedHash(*s,bad));
        same(original,copy,beta); // failed independent loads cannot mutate live state
        // A truly empty conservative domain can be restored, but is NOT a
        // feasibility witness or nonempty optimization result.
        auto outer=cube(0,1),design=cube(.125,.875);SourceSurfaceOctreeConfig c;c.maximumDepth=0;
        auto emptySource=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
        Domain empty(emptySource,0,0,1);
        auto emptyCopy=Domain::restorePartition(emptySource,0,0,1,{},empty.partitionHash());
        require(emptyCopy.leaves().empty(),"empty source was filled during restoration");
        std::cout<<"partition restoration, analytic moments, continued transitions and rejection checks passed\n";
    }
}
int main(){try{restorationTests<SourceAdaptiveDomain>();return 0;}
catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
