#include "spin2014/SourceCageFem.h"
#include "source_surface_test_utils.h"
#include <Eigen/Geometry>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <limits>
using namespace spin2014;
using namespace surface_test;
using Key=SourceCageGrid::Key;

SourceCageGrid grid(const Mesh& m,int cageLevel=3,int fineLevel=4,
                    Eigen::Vector3d origin=Eigen::Vector3d::Zero(),double side=1) {
    SourceSurfaceCageConfig c; c.cageLevel=cageLevel; c.rootOrigin=origin; c.rootSide=side;
    SourceCageGridConfig f; f.fineLevel=fineLevel;
    return {SourceSurfaceCage(m.v,m.f,c),f};
}
void near(double a,double b,double tol=2e-12) {
    require(std::abs(a-b)<=tol*std::max(1.,std::abs(b)),"independent analytic value mismatch");
}
int node(const SourceCageGrid& g,Key k) {
    const auto it=std::lower_bound(g.nodeKeys().begin(),g.nodeKeys().end(),k);
    require(it!=g.nodeKeys().end() && *it==k,"fixture node exists");
    return int(it-g.nodeKeys().begin());
}
const SourceCageFem::BoundaryRow& boundary(const SourceCageFem& f,int index) {
    const auto& rows=f.boundaryConstraints();
    const auto it=std::lower_bound(rows.begin(),rows.end(),index,[](const auto& r,int n){return r.nodeIndex<n;});
    require(it!=rows.end() && it->nodeIndex==index,"complete boundary row exists"); return *it;
}
// A missing boundary row, wrong tet split/orientation, stiffness sign or mass
// scaling must fail these independent topology and variational checks.
void topology(const SourceCageFem& f) {
    const auto& g=f.grid();
    require(f.tetrahedra().size()==6*g.cells().size(),"six tetrahedra per actual cell");
    std::map<std::array<int,3>,int> faces;
    double volume=0;
    for(std::size_t t=0;t<f.tetrahedra().size();++t) {
        const auto& tet=f.tetrahedra()[t]; Eigen::Matrix3d e;
        for(int a=0;a<3;++a) e.col(a)=(g.nodePositions().row(tet[a+1])-g.nodePositions().row(tet[0])).transpose();
        require(e.determinant()>0 && f.tetrahedronVolumes()[t]>0,"strict positive orientation and volume");
        near(f.tetrahedronVolumes()[t],e.determinant()/6); volume+=f.tetrahedronVolumes()[t];
        for(int skip=0;skip<4;++skip) {
            std::array<int,3> face{}; int j=0;
            for(int a=0;a<4;++a) if(a!=skip) face[j++]=tet[a];
            std::sort(face.begin(),face.end()); ++faces[face];
        }
    }
    double unionVolume=0;
    for(const auto& c:g.cellNodes()) unionVolume+=(g.nodePositions().row(c[7])-g.nodePositions().row(c[0])).prod();
    near(volume,unionVolume); near(f.lumpedMass().sum(),unionVolume);
    std::set<int> boundaryNodes;
    for(const auto& item:faces) {
        require(item.second==1 || item.second==2,"conforming interior faces have multiplicity two");
        if(item.second==1) for(int n:item.first) boundaryNodes.insert(n);
    }
    require(boundaryNodes.size()==f.boundaryConstraints().size(),"every and only tetrahedral union boundary node constrained");
    for(const auto& row:f.boundaryConstraints()) require(boundaryNodes.count(row.nodeIndex),"interior unconstrained");
    require(f.lumpedMass().allFinite() && f.lumpedMass().minCoeff()>0,"positive finite lumped mass");
    const SourceCageFem::SparseMatrix delta=f.stiffness()-SourceCageFem::SparseMatrix(f.stiffness().transpose());
    near(delta.norm(),0); near((f.stiffness()*Eigen::VectorXd::Ones(g.nodePositions().rows())).norm(),0);
    Eigen::MatrixXd ones=Eigen::MatrixXd::Ones(g.nodePositions().rows(),1);
    near(f.applyQ(ones).norm(),0,1e-9); near(f.energy(ones),0);
    near(f.energy(Eigen::MatrixXd::Zero(ones.rows(),1)),0);
}
void constraints(const SourceCageFem& f) {
    const auto& g=f.grid(); const auto& c=g.cage();
    int previous=-1;
    for(const auto& row:f.boundaryConstraints()) {
        require(row.nodeIndex>previous,"sorted unique boundary rows"); previous=row.nodeIndex;
        require(!row.entries.empty() && row.entries.size()<=3,"sparse triangle support");
        double sum=0; Eigen::Matrix<long double,3,1> residual=Eigen::Matrix<long double,3,1>::Zero(); int last=-1;
        for(const auto& e:row.entries) {
            require(e.handleIndex>last && std::size_t(e.handleIndex)<c.handleVertexIndices().size(),"sorted actual handle support");
            last=e.handleIndex; require(e.weight>0 && std::isfinite(e.weight),"nonzero finite nonnegative weights");
            sum+=e.weight;
            residual+=static_cast<long double>(e.weight)*(c.boundaryVertices().row(c.handleVertexIndices()[e.handleIndex]).transpose().cast<long double>()-g.nodePositions().row(row.nodeIndex).transpose().cast<long double>());
        }
        near(sum,1); near(static_cast<double>(residual.norm()),0);
        bool supported=false;
        for(int t=0;t<c.boundaryFaces().rows();++t) {
            bool all=true;
            for(const auto& e:row.entries) all &= (c.boundaryFaces().row(t).array()==c.handleVertexIndices()[e.handleIndex]).any();
            if(all) {supported=true; break;}
        }
        require(supported,"each row references corners of one supplied coarse triangle");
    }
    for(std::size_t h=0;h<g.handleNodeIndices().size();++h) {
        const auto& r=boundary(f,g.handleNodeIndices()[h]);
        require(r.entries.size()==1 && r.entries[0].handleIndex==int(h) && r.entries[0].weight==1,"every handle one-hot");
    }
}
// Independent brick formula: integrate the six path simplices symbolically.
// Reference right simplex: gradients (-1,-1,-1),(1,0,0),(0,1,0),
// (0,0,1), volume1/6, K=gradientGram/6, nodal barycentric mass1/24.
// Mapping these gradients onto each brick path and summing gives the literal
// edge formula below. Cube edges incident to000/111 carry -1/3, others -1/6.
// Diagonals carry zero. This does not inspect production tets or assembly.
std::pair<Eigen::MatrixXd,Eigen::VectorXd> reference(const SourceCageGrid& g,SourceCageFemMassRule rule) {
    const int n=g.nodePositions().rows();
    Eigen::MatrixXd k=Eigen::MatrixXd::Zero(n,n); Eigen::VectorXd m=Eigen::VectorXd::Zero(n);
    const int edges[12][4]={{0,1,0,2},{0,2,1,2},{0,4,2,2},{3,7,2,2},{5,7,1,2},{6,7,0,2},
                           {1,3,1,1},{1,5,2,1},{2,3,0,1},{2,6,2,1},{4,5,0,1},{4,6,1,1}};
    for(const auto& c:g.cellNodes()) {
        const Eigen::Vector3d d=(g.nodePositions().row(c[7])-g.nodePositions().row(c[0])).transpose();
        const double v=d.prod();
        for(int b=0;b<8;++b) m(c[b])+=v*(rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI?1./8:((b==0 || b==7)?1./4:1./12));
        for(const auto& e:edges) {
            double w=e[3]*v/(6*d[e[2]]*d[e[2]]); int a=c[e[0]],b=c[e[1]];
            k(a,a)+=w; k(b,b)+=w; k(a,b)-=w; k(b,a)-=w;
        }
    }
    return {k,m};
}
template<class F> void resource(F run) {
    bool ok=false; try {run();} catch(const std::length_error& e) {ok=std::string(e.what()).find("RESOURCE_LIMIT")!=std::string::npos;}
    require(ok,"specific RESOURCE_LIMIT length_error required");
}
template<class E,class F> void throws(F run) {
    bool ok=false; try {run();} catch(const E&) {ok=true;} require(ok,"specific exception type required");
}
Mesh tunnel() {
    Mesh m; m.v.resize(16,3); m.f.resize(32,3);
    const double rings[2][4][2]={{{.1,.1},{.9,.1},{.9,.9},{.1,.9}},{{.3,.3},{.7,.3},{.7,.7},{.3,.7}}};
    for(int z=0;z<2;++z) for(int r=0;r<2;++r) for(int j=0;j<4;++j)
        m.v.row(8*z+4*r+j)<<rings[r][j][0],rings[r][j][1],z?.9:.1;
    int row=0; auto quad=[&](int a,int b,int c,int d){m.f.row(row++)<<a,b,c; m.f.row(row++)<<a,c,d;};
    for(int j=0;j<4;++j) {int k=(j+1)%4; quad(j,k,k+8,j+8); quad(j+4,j+12,k+12,k+4); quad(j,j+4,k+4,k); quad(j+8,k+8,k+12,j+12);}
    return m;
}
int main() {
    int failures=0;
    auto test=[&](const char* name,const std::function<void()>& run) {
        try {run(); std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& e) {++failures; std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';}
    };
    const auto g=grid(cube(.25,.75)); const SourceCageFem f(g,{});
    test("full_selected_union_fem_and_complete_boundary_nodes",[&] {
        require(g.cells().size()==1728 && g.nodeKeys().size()==2197,"fixture cells and nodes");
        require(g.cage().boundaryFaces().rows()==432 && g.handleNodeIndices().size()==218,"fixture coarse cage");
        require(f.tetrahedra().size()==10368 && f.boundaryConstraints().size()==866,"complete 866 boundary nodes, not just 218 handles");
        require(g.nodeKeys().size()-f.boundaryConstraints().size()==1331,"1331 unconstrained interior nodes");
        topology(f); constraints(f);
    });
    for(auto rule:{SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,SourceCageFemMassRule::BARYCENTRIC_TET}) {
    SourceCageFemConfig massConfig; massConfig.massRule=rule;
    test(rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI?"voronoi_independent_brick_and_energy":"barycentric_independent_brick_and_energy",[&] {
        // All actual production cells assembled against independent brick formula.
        const auto tiny=grid(cube(.3,.4),2,3); const SourceCageFem tf(tiny,massConfig);
        const auto r=reference(tiny,rule);
        near((Eigen::MatrixXd(tf.stiffness())-r.first).norm(),0);
        near((tf.lumpedMass()-r.second).norm(),0);
        // At outer cube corner only one brick contributes: 1/4, not Voronoi1/8.
        near(tf.lumpedMass()(0)/std::pow(.125,3),rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI?.125:.25);
        near(tf.lumpedMass()(node(tiny,{2,2,4}))/std::pow(.125,3),rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI?.125:1./12);
        if(rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI) for(int i=0;i<27;++i) {
            // Independent Cartesian dual widths, clipped to [.25,.5]^3.
            double dual=1;
            for(int a=0;a<3;++a) dual*=tiny.nodeKeys()[i][a]==3?.125:.0625;
            near(tf.lumpedMass()(i),dual);
        }
        topology(tf); constraints(tf);
        Eigen::MatrixXd x(tiny.nodePositions().rows(),2),direction(x.rows(),2);
        for(int i=0;i<x.rows();++i) for(int j=0;j<2;++j) {x(i,j)=std::sin(.3+i*.7+j); direction(i,j)=std::cos(i*.4-j);}
        const Eigen::MatrixXd q=r.first.transpose()*r.second.cwiseInverse().asDiagonal()*r.first;
        near((tf.applyQ(x)-q*x).norm(),0,1e-9);
        const double e=.5*(x.array()*(q*x).array()).sum(); near(tf.energy(x),e);
        require(tf.energy(x)>0,"PSD nontrivial energy");
        const double derivative=(direction.array()*tf.applyQ(x).array()).sum();
        for(double h:{1e-3,1e-4,1e-5}) near((tf.energy(x+h*direction)-tf.energy(x-h*direction))/(2*h),derivative,1e-8);
        near(tf.energy(x),tf.energy(x.leftCols(1))+tf.energy(x.rightCols(1)));
        Eigen::MatrixXd affine=tiny.nodePositions(); const Eigen::MatrixXd residual=tf.stiffness()*affine;
        near(residual.row(node(tiny,{3,3,3})).norm(),0);
        for(int axis=0;axis<3;++axis) near(affine.col(axis).dot(residual.col(axis)),std::pow(.25,3));
    });
    test(rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI?"voronoi_actual_coordinate_translation_and_scale_covariance":"barycentric_actual_coordinate_translation_and_scale_covariance",[&] {
        const auto base=grid(cube(.3,.4),2,3); const SourceCageFem a(base,massConfig);
        for(double s:{1.,1.7,4.}) {
            const Eigen::Vector3d origin(-.3,-.2,-.1);
            auto mesh=cube(.3,.4); for(int i=0;i<mesh.v.rows();++i) for(int j=0;j<3;++j) mesh.v(i,j)=std::fma(mesh.v(i,j),s,origin[j]);
            const auto shifted=grid(mesh,2,3,origin,s); const SourceCageFem b(shifted,massConfig);
            require(base.nodeKeys()==shifted.nodeKeys(),"same integer geometry under transform");
            near((Eigen::MatrixXd(b.stiffness())-s*Eigen::MatrixXd(a.stiffness())).norm(),0);
            near((b.lumpedMass()-s*s*s*a.lumpedMass()).norm(),0);
            Eigen::MatrixXd x(base.nodePositions().rows(),1); for(int i=0;i<x.rows();++i)x(i)=i%5;
            near((b.applyQ(x)-a.applyQ(x)/s).norm(),0,1e-9); near(b.energy(x),a.energy(x)/s);
            const auto r=reference(shifted,rule);
            near((Eigen::MatrixXd(b.stiffness())-r.first).norm(),0); near((b.lumpedMass()-r.second).norm(),0);
            topology(b); constraints(b);
        }
    });
    }
    test("large_shifted_nondyadic_root_requires_actual_unequal_intervals",[&] {
        const Eigen::Vector3d origin(std::ldexp(1.,45),std::ldexp(1.,44),-std::ldexp(1.,45));
        auto mesh=cube(.3,.4); for(int i=0;i<8;++i)for(int a=0;a<3;++a)mesh.v(i,a)=std::fma(mesh.v(i,a),1.7,origin[a]);
        const auto actual=grid(mesh,2,3,origin,1.7);
        require(actual.cells().size()==8 && actual.nodeKeys().size()==27,"large shifted fixture is valid upstream");
        require(actual.nodePositions()(node(actual,{3,2,2}),0)-actual.nodePositions()(node(actual,{2,2,2}),0)!=
                actual.nodePositions()(node(actual,{4,2,2}),0)-actual.nodePositions()(node(actual,{3,2,2}),0),"actual fine intervals unequal, not nominal step");
        for(auto rule:{SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,SourceCageFemMassRule::BARYCENTRIC_TET}) {
            SourceCageFemConfig c; c.massRule=rule; const SourceCageFem tf(actual,c); const auto r=reference(actual,rule);
            near((Eigen::MatrixXd(tf.stiffness())-r.first).norm(),0); near((tf.lumpedMass()-r.second).norm(),0);
            topology(tf); constraints(tf);
            if(rule==SourceCageFemMassRule::RESTRICTED_GRID_VORONOI)for(int i=0;i<27;++i) {
                long double dual=1;
                for(int a=0;a<3;++a) {
                    const auto index=actual.nodeKeys()[i][a]; Key lo=actual.nodeKeys()[i],hi=lo;
                    if(index>2)--lo[a]; if(index<4)++hi[a];
                    dual*=(static_cast<long double>(actual.nodePositions()(node(actual,hi),a))-actual.nodePositions()(node(actual,lo),a))/2;
                }
                near(tf.lumpedMass()(i),static_cast<double>(dual));
            }
        }
    });
    test("non_handle_face_edge_diagonal_rows_and_offnode_transfer_mismatch",[&] {
        for(const Key key: {Key{2,2,3},Key{2,3,3}}) {
            const auto& row=boundary(f,node(g,key)); require(row.entries.size()==2,"edge/diagonal midpoint has two coarse handles");
            near(row.entries[0].weight,.5); near(row.entries[1].weight,.5);
        }
        const auto fine=grid(cube(.3,.4),2,4); const SourceCageFem ff(fine,{}); constraints(ff);
        const auto& row=boundary(ff,node(fine,{4,5,6}));
        require(row.entries.size()==3,"off-diagonal face node has three corner support");
        std::vector<double> weights; for(const auto& e:row.entries) weights.push_back(e.weight);
        std::sort(weights.begin(),weights.end()); near(weights[0],.25); near(weights[1],.25); near(weights[2],.5);
        // Locate the actual coarse low-x face's (y high,z low) handle. Its
        // triangular weight is max(u-v,0), while trilinear gives 1/8 here.
        const auto tiny=grid(cube(.3,.4),2,3); const SourceCageFem tf(tiny,{});
        int handle=-1; for(std::size_t h=0;h<tiny.handleNodeIndices().size();++h)
            if(tiny.nodeKeys()[tiny.handleNodeIndices()[h]]==Key{2,4,2})handle=int(h);
        require(handle>=0,"analytic corner handle found");
        Eigen::MatrixXd field=Eigen::MatrixXd::Zero(tiny.nodePositions().rows(),1);
        for(const auto& row:tf.boundaryConstraints()) for(const auto& e:row.entries) if(e.handleIndex==handle)field(row.nodeIndex)=e.weight;
        Eigen::MatrixXd query(1,3); query<<.25,.3125,.3125;
        const double mismatch=tiny.interpolate(field,query)(0,0); near(mismatch,.125);
        std::cout<<"off-node triangular weight=0, trilinear transfer="<<mismatch<<" (NODAL constraint only)\n";
    });
    test("cavity_tunnel_and_disconnected_component_boundaries_and_nullspaces",[&] {
        for(auto rule:{SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,SourceCageFemMassRule::BARYCENTRIC_TET}) {
        SourceCageFemConfig c; c.massRule=rule;
        const SourceCageFem hollow(grid(join(cube(.1,.9),cube(.3,.7),true)),c);
        const SourceCageFem ring(grid(tunnel()),c);
        const SourceCageFem islands(grid(join(cube(.03,.09),cube(.91,.97))),c);
        topology(hollow); constraints(hollow); topology(ring); constraints(ring); topology(islands); constraints(islands);
        boundary(hollow,node(hollow.grid(),{6,8,8})); boundary(ring,node(ring.grid(),{6,8,8}));
        Eigen::MatrixXd x(islands.grid().nodePositions().rows(),2);
        for(int i=0;i<x.rows();++i) {const bool left=islands.grid().nodePositions()(i,0)<.5; x(i,0)=left?1:0; x(i,1)=left?0:1;}
        near(islands.applyQ(x).norm(),0,1e-9); near(islands.energy(x),0);
        require(hollow.boundaryConstraints().size()!=ring.boundaryConstraints().size(),"cavity and tunnel constraints remain distinct");
        }
    });
    test("sparse_level20_no_root_expansion_resource_and_evaluation_validation",[&] {
        const auto tiny=grid(cube(100.25/524288,100.75/524288),19,20);
        SourceCageFemConfig c; c.maximumTetrahedra=48; c.maximumAssemblyEntries=768;
        c.maximumStiffnessEntries=200; c.maximumBoundaryEntries=80; c.maximumEvaluationValues=54;
        const SourceCageFem tf(tiny,c); require(tf.tetrahedra().size()==48 && tf.boundaryConstraints().size()==26,"genuine sparse level20, 27nodes");
        for(int which=0;which<5;++which) for(std::size_t v:{std::size_t(0),std::numeric_limits<std::size_t>::max()}) {
            auto d=c; if(which==0)d.maximumTetrahedra=v; if(which==1)d.maximumAssemblyEntries=v; if(which==2)d.maximumStiffnessEntries=v; if(which==3)d.maximumBoundaryEntries=v; if(which==4)d.maximumEvaluationValues=v;
            resource([&]{SourceCageFem bad(tiny,d);});
        }
        for(int which=0;which<4;++which) {
            auto d=c; if(which==0)d.maximumTetrahedra=47; if(which==1)d.maximumAssemblyEntries=767;
            if(which==2)d.maximumStiffnessEntries=1; if(which==3)d.maximumBoundaryEntries=1;
            resource([&]{SourceCageFem bad(tiny,d);});
        }
        Eigen::MatrixXd field=Eigen::MatrixXd::Ones(27,1); tf.applyQ(field);
        resource([&]{tf.applyQ(Eigen::MatrixXd::Ones(27,2));}); resource([&]{tf.energy(Eigen::MatrixXd::Ones(27,2));});
        throws<std::invalid_argument>([&]{tf.applyQ(Eigen::MatrixXd(26,1));});
        throws<std::invalid_argument>([&]{tf.applyQ(Eigen::MatrixXd(27,0));});
        for(double v:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
            field(0)=v; throws<std::invalid_argument>([&]{tf.applyQ(field);}); throws<std::invalid_argument>([&]{tf.energy(field);});
        }
        field.setZero(); field(0)=std::numeric_limits<double>::max();
        throws<std::overflow_error>([&]{tf.applyQ(field);}); throws<std::overflow_error>([&]{tf.energy(field);});
    });
    test("unrepresentable_fem_contributions_reject_after_valid_upstream_grid",[&] {
        const double s=5e-108;
        const auto tiny=grid(cube(.25*s,.75*s),0,1,Eigen::Vector3d::Zero(),s);
        require(tiny.cells().size()==8 && tiny.nodeKeys().size()==27,"valid upstream subnormal-volume geometry");
        for(auto rule:{SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,SourceCageFemMassRule::BARYCENTRIC_TET}) {
            SourceCageFemConfig c; c.massRule=rule;
            throws<std::overflow_error>([&]{SourceCageFem unsupported(tiny,c);});
        }
    });
    test("immutable_inputs_full_configuration_output_hashes_and_field_semantics",[&] {
        const auto tiny=grid(cube(.3,.4),2,3); const SourceCageFem a(tiny,{}),b(tiny,{});
        require(a.configSha256().size()==64 && a.outputSha256().size()==64 && a.configSha256()==b.configSha256() && a.outputSha256()==b.outputSha256(),"deterministic hashes");
        require(a.sourceSha256()==tiny.sourceSha256() && a.gridConfigSha256()==tiny.configSha256() && a.gridOutputSha256()==tiny.outputSha256(),"immutable upstream identities");
        require(a.tetrahedra()==b.tetrahedra() && a.tetrahedronVolumes()==b.tetrahedronVolumes(),"deterministic tetrahedral arrays");
        require(a.lumpedMass()==b.lumpedMass() && (a.stiffness()-b.stiffness()).norm()==0,"deterministic computational coefficients");
        for(int which=0;which<5;++which) {
            SourceCageFemConfig c; if(which==0)++c.maximumTetrahedra; if(which==1)++c.maximumAssemblyEntries; if(which==2)++c.maximumStiffnessEntries; if(which==3)++c.maximumBoundaryEntries; if(which==4)++c.maximumEvaluationValues;
            const SourceCageFem d(tiny,c); require(a.configSha256()!=d.configSha256() && a.outputSha256()!=d.outputSha256(),"each resource parameter binds config and output");
        }
        const SourceCageFem d(grid(cube(.31,.39),2,3),{});
        require(a.grid().nodeKeys()==d.grid().nodeKeys() && a.outputSha256()!=d.outputSha256(),"different upstream source changes identity on same domain");
        Eigen::MatrixXd x=Eigen::MatrixXd::Zero(27,1); const auto hash=a.outputSha256(); const double e=a.energy(x); x(0)=1;
        require(a.energy(x)>e && a.outputSha256()==hash,"field changes energy, never immutable operator identity");
        SourceCageFemConfig c; c.massRule=SourceCageFemMassRule::BARYCENTRIC_TET; const SourceCageFem bary(tiny,c);
        require(a.configSha256()!=bary.configSha256() && a.outputSha256()!=bary.outputSha256(),"mass rule bound to hashes");
        require((a.stiffness()-bary.stiffness()).norm()==0 && (a.lumpedMass()-bary.lumpedMass()).norm()>0,"mass alternatives share K and differ in M");
        require((a.applyQ(x)-bary.applyQ(x)).norm()>0,"Q differs under explicit mass alternatives");
        require(a.boundaryConstraints().size()==bary.boundaryConstraints().size(),"mass alternatives preserve complete boundary size");
        for(std::size_t i=0;i<a.boundaryConstraints().size();++i) {
            const auto& first=a.boundaryConstraints()[i]; const auto& second=bary.boundaryConstraints()[i];
            require(first.nodeIndex==second.nodeIndex && first.entries.size()==second.entries.size(),"mass alternatives preserve sparse boundary support");
            for(std::size_t j=0;j<first.entries.size();++j)
                require(first.entries[j].handleIndex==second.entries[j].handleIndex && first.entries[j].weight==second.entries[j].weight,"mass alternatives preserve exact boundary coefficients");
        }
        c.massRule=static_cast<SourceCageFemMassRule>(100); throws<std::invalid_argument>([&]{SourceCageFem invalid(tiny,c);});
    });
    std::cout<<"source cage FEM SYNTHETIC_REGRESSION: "<<failures<<" failed groups; no BBW solve or physical mass authority\n";
    return failures?1:0;
}
