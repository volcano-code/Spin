#include "detail/RankTwoNormalMesh.h"
#include "verification/Sha256.h"
#include <algorithm>
#include <locale>
#include <map>
#include <set>
#include <sstream>
#include <utility>

namespace {
using Z=CGAL::Gmpz;
using I3=std::array<Z,3>;
using Triangle=std::array<std::size_t,3>;
using Edge=std::array<std::size_t,2>;
using Cone=std::array<I3,3>;
[[noreturn]] void invalid(const std::string& witness) {
    throw spin::RankTwoSourceError("SOURCE_TOPOLOGY_INVALID",witness);
}
bool zero(const I3& a){return a[0]==0 && a[1]==0 && a[2]==0;}
Z dot(const I3& a,const I3& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
I3 cross(const I3& a,const I3& b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
I3 primitive(I3 a){
    Z divisor(0);
    for(const auto& x:a)mpz_gcd(divisor.mpz(),divisor.mpz(),x.mpz());
    if(divisor!=0)for(auto& x:a)x/=divisor;
    return a; // divisor is positive: never reverse the ray
}
I3 integerRay(const std::array<CGAL::Gmpq,3>& q){
    Z denominator(1);I3 result;
    for(const auto& x:q)mpz_lcm(denominator.mpz(),denominator.mpz(),x.denominator().mpz());
    for(int i=0;i<3;++i)result[i]=q[i].numerator()*(denominator/q[i].denominator());
    return primitive(std::move(result));
}
Edge edge(std::size_t a,std::size_t b){return {std::min(a,b),std::max(a,b)};}
std::string hash(const std::string& text){
    verification::Sha256 sha;
    sha.update(reinterpret_cast<const unsigned char*>(text.data()),text.size());return sha.finish();
}
std::string authority(const spin::PreparedRankTwoNormalMesh::Rays& rays,
                      const spin::PreparedRankTwoNormalMesh::Faces& faces){
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<"rank-two-normal-authority-v1\nRAYS "<<rays.size()<<'\n';
    for(const auto& v:rays){for(const auto& q:v)out<<q.numerator()<<'/'<<q.denominator()<<' ';out<<'\n';}
    out<<"FACES "<<faces.size()<<'\n';
    for(const auto& f:faces)out<<f[0]<<' '<<f[1]<<' '<<f[2]<<'\n';
    return out.str();
}

void checkPair(const Cone& a,const Cone& b,const Triangle& fa,const Triangle& fb,
               std::size_t ia,std::size_t ib){
    std::array<I3,6> halfplanes{{a[0],a[1],a[2],b[0],b[1],b[2]}};
    std::array<bool,3> shared{};
    for(int i=0;i<3;++i)shared[i]=std::find(fb.begin(),fb.end(),fa[i])!=fb.end();
    // K=C_a intersect C_b is pointed. The section of K by
    // (a_0+a_1+a_2).x=1 is compact. Every vertex of that section has
    // >=2 independent active facets, so every extreme ray of K occurs
    // among the +/- cross products below, including lower-dimensional K.
    // K equals the prescribed common face iff all these rays have zero
    // barycentric coordinates at the non-shared vertices of a. The reverse
    // inclusion is automatic because shared input rays are identical.
    for(int p=0;p<6;++p)for(int q=p+1;q<6;++q){
        auto r=cross(halfplanes[p],halfplanes[q]);if(zero(r))continue;
        for(int direction=0;direction<2;++direction){
            if(direction)for(auto& x:r)x=-x;
            bool inside=true;
            for(const auto& h:halfplanes)if(dot(h,r)<0){inside=false;break;}
            if(!inside)continue;
            for(int i=0;i<3;++i)if(!shared[i] && dot(a[i],r)!=0){
                r=primitive(std::move(r));std::ostringstream witness;witness.imbue(std::locale::classic());
                witness<<"cone intersection outside common face: "<<ia<<' '<<ib
                    <<" ray "<<r[0]<<' '<<r[1]<<' '<<r[2]<<" nonshared "<<fa[i];
                invalid(witness.str());
            }
        }
    }
}
} // namespace

namespace spin {
struct PreparedRankTwoNormalMesh::Impl {
    Rays rays;Faces faces;std::vector<Edge> edges;
    std::size_t checkedPairs=0;
    std::string inputHash,record;
};
PreparedRankTwoNormalMesh::PreparedRankTwoNormalMesh(std::shared_ptr<const Impl> p):impl_(std::move(p)){}
PreparedRankTwoNormalMesh PreparedRankTwoNormalMesh::fromExact(
        const Rays& rays,const Faces& faces,const RankTwoSourceLimits& limits){
    if(rays.size()>limits.maximumVertices || faces.size()>limits.maximumArrangementCells)
        throw RankTwoSourceError("RESOURCE_LIMIT","normal mesh input budget exceeded");
    if(rays.empty() || faces.empty())invalid("empty normal mesh");
    std::vector<I3> integers;std::set<I3> distinct;
    for(std::size_t i=0;i<rays.size();++i){
        auto r=integerRay(rays[i]);
        if(zero(r))invalid("zero ray at vertex "+std::to_string(i));
        if(!distinct.insert(r).second)invalid("duplicate positive ray at vertex "+std::to_string(i));
        integers.push_back(std::move(r));
    }
    struct Incidence {std::size_t count=0;int orientation=0;};
    std::map<Edge,Incidence> incidences;
    std::set<Triangle> distinctFaces;
    std::vector<Cone> cones;
    std::vector<std::map<std::size_t,std::set<std::size_t>>> links(rays.size());
    for(std::size_t fi=0;fi<faces.size();++fi){
        const auto& f=faces[fi];
        for(const auto v:f)if(v>=rays.size())invalid("invalid vertex index in face "+std::to_string(fi));
        auto sorted=f;std::sort(sorted.begin(),sorted.end());
        if(sorted[0]==sorted[1] || sorted[1]==sorted[2])invalid("repeated vertex in face "+std::to_string(fi));
        if(!distinctFaces.insert(sorted).second)invalid("duplicate face "+std::to_string(fi));
        const auto& a=integers[f[0]];const auto& b=integers[f[1]];const auto& c=integers[f[2]];
        if(dot(a,cross(b,c))<=0)invalid("singular or reversed orientation in face "+std::to_string(fi));
        cones.push_back({primitive(cross(b,c)),primitive(cross(c,a)),primitive(cross(a,b))});
        for(int i=0;i<3;++i){
            const auto v=f[i],next=f[(i+1)%3],last=f[(i+2)%3];
            auto& incidence=incidences[edge(v,next)];
            ++incidence.count;incidence.orientation+=v<next?1:-1;
            if(!links[v][next].insert(last).second || !links[v][last].insert(next).second)
                invalid("duplicate vertex link edge at "+std::to_string(v));
        }
    }
    if(incidences.size()>limits.maximumArrangementCells-faces.size()
        || rays.size()>limits.maximumArrangementCells-faces.size()-incidences.size())
        throw RankTwoSourceError("RESOURCE_LIMIT","normal mesh cell budget exceeded");
    std::vector<std::vector<std::size_t>> adjacency(rays.size());
    for(const auto& entry:incidences){
        if(entry.second.count!=2 || entry.second.orientation!=0)
            invalid("edge incidence or orientation at "+std::to_string(entry.first[0])+" "+std::to_string(entry.first[1]));
        adjacency[entry.first[0]].push_back(entry.first[1]);adjacency[entry.first[1]].push_back(entry.first[0]);
    }
    for(std::size_t v=0;v<links.size();++v){
        const auto& link=links[v];
        if(link.empty())invalid("unused vertex "+std::to_string(v));
        if(link.size()<3)invalid("vertex link is not a simple cycle at "+std::to_string(v));
        for(const auto& p:link)if(p.second.size()!=2)invalid("vertex link degree at "+std::to_string(v));
        std::set<std::size_t> visited;std::vector<std::size_t> stack{link.begin()->first};
        while(!stack.empty()){
            const auto x=stack.back();stack.pop_back();if(!visited.insert(x).second)continue;
            for(const auto y:link.at(x))stack.push_back(y);
        }
        if(visited.size()!=link.size())invalid("disconnected vertex link at "+std::to_string(v));
    }
    std::set<std::size_t> visited;std::vector<std::size_t> stack{0};
    while(!stack.empty()){
        const auto x=stack.back();stack.pop_back();if(!visited.insert(x).second)continue;
        for(const auto y:adjacency[x])stack.push_back(y);
    }
    if(visited.size()!=rays.size())invalid("disconnected normal mesh");
    if(rays.size()+faces.size()!=incidences.size()+2)invalid("Euler characteristic is not two");
    auto result=std::make_shared<Impl>();result->rays=rays;result->faces=faces;
    for(const auto& e:incidences)result->edges.push_back(e.first);
    for(std::size_t i=0;i<faces.size();++i)for(std::size_t j=i+1;j<faces.size();++j){
        checkPair(cones[i],cones[j],faces[i],faces[j],i,j);++result->checkedPairs;
    }
    // The verified simple links and edge incidences give a closed connected
    // 2-manifold. Every face maps homeomorphically to its spherical cone;
    // pairwise intersection checks give an injection. Invariance of domain
    // plus compactness then gives all of S^2, not just a locally valid cover.
    const auto input=authority(rays,faces);result->inputHash=hash(input);
    std::ostringstream out;out.imbue(std::locale::classic());
    out<<"rank-two-normal-mesh-admission-v1\nINPUT "<<result->inputHash<<'\n'<<input
        <<"ALGORITHM exact-six-halfplane-extreme-rays-v1\nEDGES "<<result->edges.size()<<'\n';
    for(const auto& e:result->edges)out<<e[0]<<' '<<e[1]<<'\n';
    out<<"CHECKED_FACE_PAIRS "<<result->checkedPairs<<"\nSCOPE embedded-oriented-normal-sphere-only\n";
    result->record=out.str();return PreparedRankTwoNormalMesh(std::move(result));
}
PreparedRankTwoNormalMesh PreparedRankTwoNormalMesh::reload(
        const Rays& rays,const Faces& faces,const std::string& saved,const RankTwoSourceLimits& limits){
    const auto result=fromExact(rays,faces,limits);
    if(result.canonicalRecord()!=saved)throw RankTwoSourceError("INPUT_IDENTITY_MISMATCH","normal mesh replay mismatch");
    return result;
}
const PreparedRankTwoNormalMesh::Rays& PreparedRankTwoNormalMesh::rays() const{return impl_->rays;}
const PreparedRankTwoNormalMesh::Faces& PreparedRankTwoNormalMesh::faces() const{return impl_->faces;}
const std::vector<Edge>& PreparedRankTwoNormalMesh::edges() const{return impl_->edges;}
RankTwoNormalFace PreparedRankTwoNormalMesh::face(std::size_t i) const{
    RankTwoNormalFace f;f.vertexIds=impl_->faces.at(i);
    for(int j=0;j<3;++j)f.rays[j]=impl_->rays.at(f.vertexIds[j]);return f;
}
std::size_t PreparedRankTwoNormalMesh::checkedFacePairs() const{return impl_->checkedPairs;}
const std::string& PreparedRankTwoNormalMesh::inputHash() const{return impl_->inputHash;}
const std::string& PreparedRankTwoNormalMesh::canonicalRecord() const{return impl_->record;}
} // namespace spin
