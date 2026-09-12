#include "spin2014/SourceConformingBoundary.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <map>
#include <set>
#include <cstring>
#include <iostream>
#include <limits>
using namespace spin2014;
using namespace surface_test;
using Boundary=SourceConformingBoundary;
using Key=Boundary::Key;
template<class F> void resource(F f) {
    bool caught=false;try {f();}catch(const std::length_error& e){caught=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;}
    require(caught,"RESOURCE_LIMIT required");
}
SourceCageBbw solve(const Mesh& outer,Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1) {
    SourceSurfaceCageConfig cc;cc.cageLevel=2;cc.rootOrigin=origin;cc.rootSide=side;
    SourceCageGridConfig gc;gc.fineLevel=3;
    return SourceCageBbw(SourceCageFem(SourceCageGrid(SourceSurfaceCage(outer.v,outer.f,cc),gc),{}),{});
}
std::shared_ptr<const SourceSurfaceOctree> oracle(const Mesh& outer,const Mesh& design,int depth=5,
    Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1) {
    SourceSurfaceOctreeConfig c;c.maximumDepth=depth;c.rootOrigin=origin;c.rootSide=side;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
}
// Independent signed simplex quadrature through a local tetrahedron origin.
Ten triangleMoments(const Eigen::MatrixXd& points,Boundary::Triangle t,int orientation,Eigen::Vector3d origin) {
    if(orientation<0)std::swap(t[1],t[2]);
    Eigen::MatrixXd v(4,3);v.row(0)=origin.transpose();
    for(int i=0;i<3;++i)v.row(i+1)=points.row(t[i]);
    Eigen::Vector3d a=v.row(1).transpose()-origin,b=v.row(2).transpose()-origin,c=v.row(3).transpose()-origin;
    return (a.dot(b.cross(c))<0?-1.:1.)*tetraMoments(v);
}
Ten checkBoundaries(const Boundary& b,const Eigen::MatrixXd& points,bool rest) {
    std::vector<int> uses(b.triangles().size(),0),balance(uses.size(),0);
    Ten sum=Ten::Zero(),weighted=Ten::Zero();
    for(std::size_t owner=0;owner<b.primitiveBoundaries().size();++owner) {
        const auto& primitive=b.binding().primitives().at(owner);
        std::map<std::array<std::size_t,2>,std::pair<int,int>> edges;
        Ten moments=Ten::Zero();const Eigen::Vector3d center=(primitive.bounds.lower+primitive.bounds.upper)/2;
        for(const auto& ref:b.primitiveBoundaries()[owner]) {
            require(ref.orientation==1||ref.orientation==-1,"valid reference orientation");
            auto tri=b.triangles().at(ref.triangle);++uses[ref.triangle];balance[ref.triangle]+=ref.orientation;
            moments+=triangleMoments(points,tri,ref.orientation,center);
            if(ref.orientation<0)std::swap(tri[1],tri[2]);
            if(rest) {
                const Eigen::Vector3d a=points.row(tri[0]).transpose(),b=points.row(tri[1]).transpose(),c=points.row(tri[2]).transpose();
                require((b-a).cross(c-a).dot(a-center)>0,"every primitive rest triangle outward and nondegenerate");
            }
            for(int k=0;k<3;++k) {
                auto a=tri[k],c=tri[(k+1)%3];const int sign=a<c?1:-1;if(a>c)std::swap(a,c);
                auto& e=edges[{a,c}];++e.first;e.second+=sign;
            }
        }
        for(const auto& e:edges)require(e.second.first==2&&e.second.second==0,"indexed primitive closure and orientation");
        if(rest)close(moments,boxMoments(primitive.bounds.lower,primitive.bounds.upper),"all ten per-primitive signed tetra moments");
        sum+=moments;weighted+=(double(primitive.betaIndex+1)/(b.binding().domain().leaves().size()+1))*moments;
    }
    Ten external=Ten::Zero();
    for(const auto& p:b.patches())for(std::size_t i=p.triangleBegin;i<p.triangleBegin+p.triangleCount;++i) {
        const bool shared=p.negativeOwner!=Boundary::noOwner&&p.positiveOwner!=Boundary::noOwner;
        require(uses.at(i)==(shared?2:1),"every actual triangle owned exactly once per occupied side");
        require(balance.at(i)==(shared?0:(p.negativeOwner!=Boundary::noOwner?1:-1)),"shared triangle cancels exactly by index");
        external+=triangleMoments(points,b.triangles()[i],1,Eigen::Vector3d::Zero())*balance[i];
    }
    close(sum,external,"deformed primitive reintegration equals independently summed shared surfaces");
    if(rest) {
        Eigen::VectorXd beta(b.binding().domain().leaves().size());for(Eigen::Index i=0;i<beta.size();++i)beta[i]=double(i+1)/(beta.size()+1);
        close(weighted,b.binding().domain().removedMoments(beta).values(),"nonuniform beta owner moments");
    }
    return sum;
}
// Tiny-fixture oracle uses included unit cells, not the production overlay.
void checkCoverage(const Boundary& b) {
    const int d=b.binding().source().config().maximumDepth;
    std::map<Key,std::size_t> cells;
    for(std::size_t i=0;i<b.binding().primitives().size();++i) {
        const auto& p=b.binding().primitives()[i].sourceKey;const auto step=1u<<(d-p.level);
        for(unsigned x=0;x<step;++x)for(unsigned y=0;y<step;++y)for(unsigned z=0;z<step;++z)
            require(cells.emplace(Key{p.coordinates[0]*step+x,p.coordinates[1]*step+y,p.coordinates[2]*step+z},i).second,"nonoverlapping included unit cells");
    }
    using Face=std::array<std::size_t,6>; // axis, plane, u, v, owner, upper
    std::set<Face> expected,emitted;
    for(const auto& c:cells)for(int axis=0;axis<3;++axis)for(int upper=0;upper<2;++upper) {
        Key neighbor=c.first;neighbor[axis]+=upper?1:unsigned(-1);const auto n=cells.find(neighbor);
        if(n!=cells.end()&&n->second==c.second)continue;
        expected.insert({std::size_t(axis),c.first[axis]+std::size_t(upper),c.first[(axis+1)%3],c.first[(axis+2)%3],c.second,std::size_t(upper)});
    }
    for(const auto& p:b.patches()) {
        const auto step=1u<<(d-p.level);
        for(unsigned u=0;u<step;++u)for(unsigned v=0;v<step;++v)for(int upper=0;upper<2;++upper) {
            const auto owner=upper?p.negativeOwner:p.positiveOwner;if(owner==Boundary::noOwner)continue;
            require(emitted.insert({std::size_t(p.axis),p.plane,p.coordinates[0]*step+u,p.coordinates[1]*step+v,owner,std::size_t(upper)}).second,"patch coverage never repeats");
        }
    }
    require(expected==emitted,"independent included-unit-face coverage and orientations");
}
void checkPoints(const Boundary& b) {
    require(b.points().topRows(b.outerPointCount())==b.binding().source().outerVertices(),"original outer unchanged");
    const auto& keys=b.materialPointKeys();require(std::is_sorted(keys.begin(),keys.end()),"sorted keys");
    require(std::set<Key>(keys.begin(),keys.end()).size()==keys.size(),"unique keys");
    for(auto key:b.binding().materialCornerKeys()) {for(auto& x:key)x*=2;require(std::binary_search(keys.begin(),keys.end(),key),"every original primitive corner retained");}
    for(std::size_t i=0;i<keys.size();++i)for(int a=0;a<3;++a)
        require(b.points()(b.outerPointCount()+i,a)==std::fma(std::ldexp(double(keys[i][a]),-b.latticeLevel()),b.binding().source().config().rootSide,b.binding().source().config().rootOrigin[a]),"canonical FMA points");
}
int main() {
    try {
        const auto outer=cube(.3,.7),design=cube(.34,.66);const auto src=oracle(outer,design);
        SourceAdaptiveDomain domain(src,2,5);const auto bbw=solve(outer);require(bbw.success(),"actual BBW successful");
        const SourceCageBinding binding(domain,bbw,{});const auto before=binding.outputSha256();
        const Boundary b(binding,{}),again(binding,{});
        require(binding.primitives().size()==496&&domain.leaves().size()==8&&binding.handleCount()==26,"unaltered actual 496-primitive fixture");
        checkPoints(b);checkCoverage(b);const auto total=checkBoundaries(b,b.points(),true);
        require(std::abs(total[0]-std::pow(10./32,3))<1e-14,"literal source volume");
        require(b.points()==again.points()&&b.materialPointKeys()==again.materialPointKeys()&&b.triangles()==again.triangles()&&b.outputSha256()==again.outputSha256()&&b.configSha256()==again.configSha256(),"deterministic complete topology identity");
        require(std::memcmp(b.points().data(),again.points().data(),std::size_t(b.points().size())*sizeof(double))==0,"byte-identical repeated point arrays");
        for(std::size_t i=0;i<b.patches().size();++i) {
            const auto& p=b.patches()[i];const auto& q=again.patches()[i];
            require(p.axis==q.axis&&p.plane==q.plane&&p.level==q.level&&p.coordinates==q.coordinates&&p.negativeOwner==q.negativeOwner&&p.positiveOwner==q.positiveOwner&&p.triangleBegin==q.triangleBegin&&p.triangleCount==q.triangleCount,"deterministic patch topology");
        }
        for(std::size_t i=0;i<b.primitiveBoundaries().size();++i) {
            const auto& refs=b.primitiveBoundaries()[i];const auto& repeated=again.primitiveBoundaries()[i];require(refs.size()==repeated.size(),"deterministic primitive reference count");
            for(std::size_t j=0;j<refs.size();++j)require(refs[j].triangle==repeated[j].triangle&&refs[j].orientation==repeated[j].orientation,"deterministic oriented references");
        }
        bool coarseFine=false,hanging=false;
        for(const auto& p:b.patches()) {
            hanging|=p.triangleCount>4;
            if(p.negativeOwner!=Boundary::noOwner&&p.positiveOwner!=Boundary::noOwner)
                coarseFine|=binding.primitives()[p.negativeOwner].sourceKey.level!=binding.primitives()[p.positiveOwner].sourceKey.level;
        }
        require(coarseFine&&hanging,"actual coarse/fine face and inserted edge breakpoint");
        Eigen::MatrixXd moved=b.points();for(Eigen::Index i=0;i<moved.rows();++i)moved(i,2)+=moved(i,0)*moved(i,1)/16;
        checkBoundaries(b,moved,false);
        // Naive 0..1 coarse diagonal interpolates xy/16 at center as 1/32;
        // refined shared center is 1/64. This is a synthetic counterexample.
        require((0.+1./16)/2-(.5*.5/16)==1./64,"naive coarse/fine center mismatch");
        const auto& grid=bbw.fem().grid();
        for(Eigen::Index begin=0;begin<moved.rows();begin+=11) {
            const Eigen::Index count=std::min(Eigen::Index(11),moved.rows()-begin);
            const Eigen::MatrixXd weights=grid.interpolate(bbw.finalWeights(),b.points().middleRows(begin,count));
            require(weights.cols()==26,"all actual BBW handles used");
            for(Eigen::Index i=0;i<count;++i) {
                Eigen::Vector3d q=Eigen::Vector3d::Zero(),p=b.points().row(begin+i).transpose();
                for(int h=0;h<26;++h)q+=weights(i,h)*((1.+.002*std::sin(h+1))*p+Eigen::Vector3d(.003*std::sin(h),.002*std::cos(2*h),.004*std::sin(3*h)));
                moved.row(begin+i)=q.transpose();
            }
        }
        checkBoundaries(b,moved,false);require((moved-b.points()).norm()>1e-4,"actual nonuniform BBW motion");
        for(int budget=0;budget<4;++budget)for(auto value:{std::size_t(0),std::numeric_limits<std::size_t>::max(),std::size_t(1)}) {
            SourceConformingBoundaryConfig c;auto* field=budget==0?&c.maximumFaceTreeNodes:budget==1?&c.maximumPatches:budget==2?&c.maximumPoints:&c.maximumTriangles;*field=value;
            resource([&]{Boundary rejected(binding,c);});
        }
        for(int budget=0;budget<4;++budget) {
            SourceConformingBoundaryConfig c;auto* field=budget==0?&c.maximumFaceTreeNodes:budget==1?&c.maximumPatches:budget==2?&c.maximumPoints:&c.maximumTriangles;++*field;
            Boundary changed(binding,c);require(changed.configSha256()!=b.configSha256()&&changed.outputSha256()!=b.outputSha256(),"each budget binds identities");
        }
        SourceConformingBoundaryConfig exactCaps;exactCaps.maximumPoints=b.points().rows();exactCaps.maximumPatches=b.patches().size();exactCaps.maximumTriangles=b.triangles().size();
        const Boundary exactCapBoundary(binding,exactCaps);require(exactCapBoundary.triangles()==b.triangles(),"inclusive exact caps preserve output topology");
        for(int cap=0;cap<3;++cap) {
            auto c=exactCaps;if(cap==0)--c.maximumPoints;if(cap==1)--c.maximumPatches;if(cap==2)--c.maximumTriangles;
            resource([&]{Boundary rejected(binding,c);});
        }
        Eigen::VectorXd beta=Eigen::VectorXd::Constant(8,.5);domain.splitFractionalLeaves(beta,1e-6);
        const SourceCageBinding splitBinding(domain,bbw,{});const Boundary split(splitBinding,{});checkCoverage(split);close(checkBoundaries(split,split.points(),true),total,"split geometric union");
        require(split.outputSha256()!=b.outputSha256()&&split.binding().domain().leaves().size()>8,"new binding updates owners");
        beta.setOnes();domain.mergeBinarySiblings(beta,1e-6);const Boundary merged(SourceCageBinding(domain,bbw,{}),{});
        checkCoverage(merged);close(checkBoundaries(merged,merged.points(),true),total,"merge geometric union");
        require(binding.outputSha256()==before&&b.binding().domain().leaves().size()==8,"upstream immutable snapshot");
        const Eigen::Vector3d shift(1.25,-.75,.125);const double scale=2;
        auto otherOuter=outer;otherOuter.v=(outer.v*scale).rowwise()+shift.transpose();
        const auto otherDesign=box((Eigen::Vector3d(.34,.36,.38)*scale+shift).eval(),(Eigen::Vector3d(.66,.64,.62)*scale+shift).eval());
        const auto transformedBbw=solve(otherOuter,shift,scale);SourceAdaptiveDomain td(oracle(otherOuter,otherDesign,5,shift,scale),2,5);
        const Boundary transformed(SourceCageBinding(td,transformedBbw,{}),{});checkPoints(transformed);checkCoverage(transformed);checkBoundaries(transformed,transformed.points(),true);
        // True depth 9 source, tiny design; no 512^3 allocation or solve.
        SourceAdaptiveDomain deep(oracle(outer,cube(.491,.509),9),2,9);
        const Boundary sparse(SourceCageBinding(deep,bbw,{}),{});require(sparse.latticeLevel()==10&&!sparse.triangles().empty()&&sparse.points().rows()<10000,"sparse true D9 material");
        checkPoints(sparse);checkCoverage(sparse);checkBoundaries(sparse,sparse.points(),true);
        SourceAdaptiveDomain emptyDomain(oracle(outer,cube(.499,.501)),2,5);
        const Boundary empty(SourceCageBinding(emptyDomain,bbw,{}),{});
        require(empty.binding().primitives().empty()&&empty.materialPointKeys().empty()&&empty.patches().empty()&&empty.triangles().empty()&&empty.primitiveBoundaries().empty(),"empty source never invents material");
        require(empty.points()==outer.v,"empty source preserves original outer");
        // Source D5 remains represented (step=2*ULP), but D6 center spacing
        // equals the root endpoint ULP and must conservatively reject.
        const Eigen::Vector3d precisionOrigin(std::ldexp(1.,40),0,0);const double precisionSide=std::ldexp(1.,-6);
        auto precisionOuter=outer,precisionDesign=design;
        precisionOuter.v=(outer.v*precisionSide).rowwise()+precisionOrigin.transpose();
        precisionDesign.v=(design.v*precisionSide).rowwise()+precisionOrigin.transpose();
        SourceAdaptiveDomain precisionDomain(oracle(precisionOuter,precisionDesign,5,precisionOrigin,precisionSide),2,5);
        const auto precisionBbw=solve(precisionOuter,precisionOrigin,precisionSide);
        const SourceCageBinding precisionBinding(precisionDomain,precisionBbw,{});
        bool unsupported=false;try{Boundary rejected(precisionBinding,{});}catch(const std::invalid_argument& e){unsupported=std::string(e.what()).find("D+1 spacing")!=std::string::npos;}
        require(unsupported,"actual binding accepted upstream but unsupported D+1 frame rejected");
        std::cout<<"FACT primitives="<<binding.primitives().size()<<" patches="<<b.patches().size()<<" points="<<b.points().rows()<<" triangles="<<b.triangles().size()<<" config="<<b.configSha256()<<" output="<<b.outputSha256()<<'\n';
        std::cout<<"PASS SYNTHETIC_REGRESSION shared source boundary/rest moments/nonaffine/full26handleBBW motion\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
}
