#include "spin2014/SourceAdaptiveDomain.h"
#include "manufacturing/SlicerValidation.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace spin2014 {
bool SourcePartitionKey::operator<(const SourcePartitionKey& b) const {return std::tie(level,coordinates)<std::tie(b.level,b.coordinates);}
bool SourcePartitionKey::operator==(const SourcePartitionKey& b) const {return level==b.level && coordinates==b.coordinates;}
namespace {
SourcePartitionKey child(const SourcePartitionKey& p,int k) {
    auto c=p;++c.level;for(int a=0;a<3;++a)c.coordinates[a]=2*c.coordinates[a]+((k>>a)&1);return c;
}
SourcePartitionKey parent(SourcePartitionKey k) {--k.level;for(auto& c:k.coordinates)c/=2;return k;}
SourceAdaptiveLeaf query(const SourceSurfaceOctree& s,const SourcePartitionKey& k) {
    SourceAdaptiveLeaf r;static_cast<SourcePartitionKey&>(r)=k;r.moments=s.cell(k.level,k.coordinates).moments;
    if(!r.moments.isFinite() || r.moments.mass()<0)throw std::invalid_argument("nonfinite/negative source contribution");return r;
}
void validateBeta(const Eigen::VectorXd& b,std::size_t n) {
    if(b.size()!=static_cast<Eigen::Index>(n)||!b.allFinite()||(b.array()<0).any()||(b.array()>1).any())throw std::invalid_argument("source beta must be finite in [0,1]");
}
void validateEpsilon(double e) {if(!std::isfinite(e)||e<=0||e>=.5)throw std::invalid_argument("source epsilon must be in (0,.5)");}
}
SourceAdaptiveDomain::SourceAdaptiveDomain(std::shared_ptr<const SourceSurfaceOctree> s,int initial,int maximum,std::size_t budget,EmptyPartitionTag)
    :source_(std::move(s)),initialLevel_(initial),maximumLevel_(maximum),maximumLeaves_(budget) {
    if(!source_||initial<0||maximum<initial||maximum>source_->config().maximumDepth||budget==0||budget>static_cast<std::size_t>(std::numeric_limits<int>::max()/6))throw std::invalid_argument("invalid source optimization levels/budget");
}
SourceAdaptiveDomain::SourceAdaptiveDomain(std::shared_ptr<const SourceSurfaceOctree> s,int initial,int maximum,std::size_t budget)
    :SourceAdaptiveDomain(std::move(s),initial,maximum,budget,EmptyPartitionTag{}) {
    std::function<void(SourcePartitionKey)> visit=[&](SourcePartitionKey k) {
        auto q=query(*source_,k);if(q.moments.mass()==0)return;
        if(k.level==initial){if(leaves_.size()>=budget)throw std::length_error("RESOURCE_LIMIT: initial source leaves");leaves_.push_back(q);return;}
        for(int i=0;i<8;++i)visit(child(k,i));
    };
    visit({});rebuild();
}
SourceAdaptiveDomain SourceAdaptiveDomain::restorePartition(
    std::shared_ptr<const SourceSurfaceOctree> s,int initial,int maximum,std::size_t budget,
    const std::vector<SourcePartitionKey>& keys,const std::string& expectedHash) {
    SourceAdaptiveDomain result(std::move(s),initial,maximum,budget,EmptyPartitionTag{});
    if(expectedHash.size()!=64 || !std::all_of(expectedHash.begin(),expectedHash.end(),
        [](char c){return (c>='0'&&c<='9') || (c>='a'&&c<='f');}))
        throw std::invalid_argument("invalid restored partition hash");
    if(keys.size()>budget)throw std::length_error("RESOURCE_LIMIT: restored source leaves");
    // Check before allocation, geometry queries, integer shifts, or sorting.
    // Sorting keys alone could silently assign another leaf's saved beta.
    for(std::size_t i=0;i<keys.size();++i) {
        const auto& key=keys[i];
        if(key.level<0 || key.level>maximum || (i && !(keys[i-1]<key)))
            throw std::invalid_argument("invalid/noncanonical restored source keys");
        const auto extent=std::uint32_t(1)<<key.level;
        for(auto coordinate:key.coordinates)if(coordinate>=extent)
            throw std::invalid_argument("restored source coordinate outside root");
    }
    result.leaves_.reserve(keys.size());
    for(const auto& key:keys) {
        auto leaf=query(*result.source_,key);
        if(leaf.moments.mass()<=0)throw std::invalid_argument("empty restored source leaf");
        result.leaves_.push_back(std::move(leaf));
    }
    result.rebuild(); // Requery moments, validate complete nonoverlapping cover.
    if(result.partitionHash()!=expectedHash)
        throw std::invalid_argument("restored partition identity mismatch");
    return result;
}
void SourceAdaptiveDomain::validateCoverage() const {
    std::set<SourcePartitionKey> keys;for(const auto& l:leaves_) {
        if(l.level<0||l.level>maximumLevel_||!keys.insert(l).second)throw std::logic_error("invalid/duplicate source leaf");
        const auto q=query(*source_,l);if(q.moments.mass()<=0||q.moments.values()!=l.moments.values())throw std::logic_error("invalid source leaf contribution");
    }
    std::size_t visited=0;
    std::function<void(SourcePartitionKey)> cover=[&](SourcePartitionKey k) {
        if(keys.count(k)){++visited;return;}
        if(query(*source_,k).moments.mass()==0)return;
        if(k.level>=maximumLevel_)throw std::logic_error("incomplete source coverage");
        for(int i=0;i<8;++i)cover(child(k,i));
    };
    cover({});if(visited!=keys.size())throw std::logic_error("overlapping source partition");
}
void SourceAdaptiveDomain::rebuild() {
    std::sort(leaves_.begin(),leaves_.end(),[](const auto& a,const auto& b){return a<b;});
    validateCoverage();
    // A dyadic face is (axis, finest integer plane, level, u, v, side).
    // Each fine face probes its opposite face and ancestors: O(N depth log N),
    // six stored keys per leaf, without expansion into finest face tiles.
    using Face=std::tuple<int,std::uint32_t,int,std::uint32_t,std::uint32_t,int>;
    std::map<Face,std::size_t> faces;
    const int depth=source_->config().maximumDepth;
    for(std::size_t i=0;i<leaves_.size();++i){const auto& l=leaves_[i];for(int a=0;a<3;++a)for(int side=0;side<2;++side)
        faces.emplace(Face{a,(l.coordinates[a]+side)<<(depth-l.level),l.level,l.coordinates[(a+1)%3],l.coordinates[(a+2)%3],side},i);}
    std::set<std::pair<std::size_t,std::size_t>> pairs;
    for(const auto& item:faces){auto [axis,plane,level,u,v,side]=item.first;
        for(int candidate=level;candidate>=0;--candidate){auto found=faces.find(Face{axis,plane,candidate,u,v,1-side});
            if(found!=faces.end()){if(item.second!=found->second)pairs.emplace(std::min(item.second,found->second),std::max(item.second,found->second));break;}u/=2;v/=2;}}
    edges_.assign(pairs.begin(),pairs.end());
    std::vector<Eigen::Triplet<double>> triplets;triplets.reserve(4*edges_.size());
    for(auto e:edges_){triplets.emplace_back(e.first,e.first,1);triplets.emplace_back(e.second,e.second,1);triplets.emplace_back(e.first,e.second,-1);triplets.emplace_back(e.second,e.first,-1);}
    laplacian_.resize(leaves_.size(),leaves_.size());laplacian_.setFromTriplets(triplets.begin(),triplets.end());
    std::ostringstream out;out.imbue(std::locale::classic());out<<std::hexfloat<<"source-adaptive-partition-v1 "<<source_->partitionSha256()<<' '<<initialLevel_<<' '<<maximumLevel_<<' '<<maximumLeaves_<<' ';
    for(const auto& l:leaves_){out<<l.level<<':'<<l.coordinates[0]<<','<<l.coordinates[1]<<','<<l.coordinates[2]<<';';for(int j=0;j<10;++j)out<<l.moments.values()(j)<<',';}
    hash_=manufacturing::computeContentSha256(out.str());
}
physics::RawMassMoments10 SourceAdaptiveDomain::removedMoments(const Eigen::VectorXd& b) const {
    validateBeta(b,leaves_.size());physics::RawMassMoments10 r;for(std::size_t i=0;i<leaves_.size();++i)r+=leaves_[i].moments.scaledDensity(b(i));
    if(!r.isFinite())throw std::overflow_error("nonfinite removed moments");return r;
}
SourcePartitionTransition SourceAdaptiveDomain::commit(SourceAdaptiveDomain& candidate,Eigen::VectorXd& beta,Eigen::VectorXd& next,std::vector<SourcePartitionChange> changes) {
    candidate.rebuild();SourcePartitionTransition t;t.beforeHash=hash_;t.afterHash=candidate.hash_;t.changes=std::move(changes);
    t.removedBefore=removedMoments(beta);t.removedAfter=candidate.removedMoments(next);
    // Every potentially allocating/throwing operation finishes before these swaps.
    leaves_.swap(candidate.leaves_);edges_.swap(candidate.edges_);laplacian_.swap(candidate.laplacian_);hash_.swap(candidate.hash_);beta.swap(next);return t;
}
SourcePartitionTransition SourceAdaptiveDomain::splitFractionalLeaves(Eigen::VectorXd& beta,double epsilon) {
    validateBeta(beta,leaves_.size());validateEpsilon(epsilon);
    auto candidate=*this;candidate.leaves_.clear();std::map<SourcePartitionKey,double> values;std::vector<SourcePartitionChange> changes;
    auto add=[&](const SourceAdaptiveLeaf& l,double b){if(candidate.leaves_.size()>=maximumLeaves_)throw std::length_error("RESOURCE_LIMIT: adaptive source leaves");candidate.leaves_.push_back(l);values.emplace(l,b);};
    for(std::size_t i=0;i<leaves_.size();++i){const auto& l=leaves_[i];const double b=beta(i);
        if(b>=epsilon&&b<=1-epsilon&&l.level<maximumLevel_){SourcePartitionChange change;change.parent=l;
            for(int k=0;k<8;++k){auto c=query(*source_,child(l,k));if(c.moments.mass()>0){add(c,b);change.children.push_back(c);}}changes.push_back(std::move(change));
        }else add(l,b);}
    std::sort(candidate.leaves_.begin(),candidate.leaves_.end(),[](const auto&a,const auto&b){return a<b;});
    Eigen::VectorXd next(candidate.leaves_.size());for(std::size_t i=0;i<candidate.leaves_.size();++i)next(i)=values.at(candidate.leaves_[i]);
    return commit(candidate,beta,next,std::move(changes));
}
SourcePartitionTransition SourceAdaptiveDomain::mergeBinarySiblings(Eigen::VectorXd& beta,double epsilon) {
    validateBeta(beta,leaves_.size());validateEpsilon(epsilon);
    std::map<SourcePartitionKey,double> values;for(std::size_t i=0;i<leaves_.size();++i)values.emplace(leaves_[i],beta(i));
    std::vector<SourcePartitionChange> changes;
    for(int level=maximumLevel_;level>0;--level){std::set<SourcePartitionKey> parents;for(auto& item:values)if(item.first.level==level)parents.insert(parent(item.first));
        for(auto p:parents){SourcePartitionChange change;change.parent=p;bool complete=true,first=true;double common=0;
            for(int k=0;k<8;++k){auto c=child(p,k);if(query(*source_,c).moments.mass()==0)continue;
                auto it=values.find(c);if(it==values.end()){complete=false;break;}
                if(first){common=it->second;first=false;}else if(common!=it->second){complete=false;break;}
                change.children.push_back(c);}
            if(complete&&!first&&(common<epsilon||common>1-epsilon)){for(auto c:change.children)values.erase(c);values.emplace(p,common);changes.push_back(std::move(change));}
        }}
    auto candidate=*this;candidate.leaves_.clear();Eigen::VectorXd next(values.size());std::size_t i=0;
    for(auto item:values){candidate.leaves_.push_back(query(*source_,item.first));next(i++)=item.second;}
    return commit(candidate,beta,next,std::move(changes));
}
} // namespace spin2014
