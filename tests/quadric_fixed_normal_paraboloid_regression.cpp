#include "spin/QuadricSourceMotion.h"
#include "spin/detail/ProjectedDomainClearance.h"
#include "verification/ExactBoundaryAudit.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace {
void require(bool ok,const std::string& why) {if(!ok)throw std::runtime_error(why);}
template<class T,class=void> struct HasFixedNormal:std::false_type {};
template<class T> struct HasFixedNormal<T,std::void_t<decltype(T{}.fixedNormalParaboloid)>>:std::true_type {};
template<class T> void enable(T& settings) {
    if constexpr(HasFixedNormal<T>::value)settings.fixedNormalParaboloid=true;
    else throw std::runtime_error("explicit fixed-normal paraboloid path is missing");
}
Eigen::MatrixXd box(double radius) {
    Eigen::MatrixXd v(8,3);
    for(int i=0;i<8;++i)for(int j=0;j<3;++j)v(i,j)=(i&(1<<j))?radius:-radius;
    return v;
}
void projectedCoverage() {
    const std::vector<Eigen::Vector3d> concave={{0,0,0},{3,0,0},{3,1,0},{1,1,0},{1,3,0},{0,3,0}};
    Eigen::MatrixXd point(1,3);point<<2,.5,4;
    require(spin::detail::ProjectedDomainClearance(point,Eigen::Vector3d::UnitZ(),1e-10).inside(concave),
            "single-point design bound lost coverage compatibility");
    Eigen::MatrixXd line(2,3);line<<1.9,.4,-1,2.1,.6,1;
    const spin::detail::ProjectedDomainClearance lineDomain(line,Eigen::Vector3d::UnitZ(),1e-10);
    require(lineDomain.inside(concave),"two-point design bound lost coverage compatibility");
    require(lineDomain.separated({{0,0,0},{1,0,0},{0,1,0}}),"separated cap rejected for line-shaped projected design");
    require(!lineDomain.separated({{1,.1,0},{3,.1,0},{2,.9,0}}),"cap containing projected design line accepted");
    Eigen::MatrixXd domain=box(.1);domain.rowwise()+=Eigen::RowVector3d(2,.5,0);
    spin::detail::ProjectedDomainClearance inside(domain,Eigen::Vector3d::UnitZ(),1e-10);
    require(inside.inside(concave),"inside L-shaped arm is incorrectly required to lie in polygon kernel");
    // Both end clusters are inside the L, but their convex hull crosses its
    // missing corner: checking only design vertices is insufficient.
    Eigen::MatrixXd crossing(8,3);
    crossing<<2.7,.2,0,2.8,.2,0,2.8,.3,0,2.7,.3,0,
              .2,2.7,0,.3,2.7,0,.3,2.8,0,.2,2.8,0;
    require(!spin::detail::ProjectedDomainClearance(crossing,Eigen::Vector3d::UnitZ(),1e-10).inside(concave),
            "design hull crossing a concavity was accepted from vertex-only containment");
    domain=box(.5);domain.rowwise()+=Eigen::RowVector3d(2,.5,0);
    require(!spin::detail::ProjectedDomainClearance(domain,Eigen::Vector3d::UnitZ(),1e-10).inside(concave),
            "design touching the finite rim was accepted");
    const std::vector<Eigen::Vector3d> bowtie={{0,0,0},{3,3,0},{0,3,0},{3,0,0}};
    require(!inside.inside(bowtie),"self-crossing projected source boundary was accepted");
    require(inside.separated({{1,2,0},{1,3,0},{1,2,5}}),"finite lateral segment mistaken for its infinite supporting line");
    require(!inside.separated({{1.5,.5,0},{2.5,.5,0},{1.5,.5,5}}),"projected cap crosses the design hull");
    require(!inside.separated({{2,.5,0},{2,.5,5},{2,.5,9}}),"cap projecting to an interior point was accepted");
    domain=box(.1);domain.rowwise()+=Eigen::RowVector3d(2,.5,0);
    spin::detail::ProjectedDomainClearance margin(domain,Eigen::Vector3d(0,0,7),.01);
    require(!margin.separated({{2.105,.2,0},{2.105,.8,0},{2.105,.2,5}}),"positive but sub-margin clearance was accepted");
    require(margin.separated({{2.12,.2,0},{2.12,.8,0},{2.12,.2,5}}),"clearance depends on arbitrary axis length");
}
void localGraphCertificate() {
    const std::vector<std::array<int,3>> faces={{0,1,2}};
    const std::vector<int> rim={0,1,2};
    const std::vector<Eigen::Vector3d> sloped={{-100,-100,-90},{100,-100,110},{0,100,10}};
    for(bool oblique:{false,true}) {
        Eigen::Matrix3d rotation=Eigen::Matrix3d::Identity();
        if(oblique)rotation=Eigen::AngleAxisd(.73,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
        const Eigen::Vector3d shift(3,-5,2);
        Eigen::MatrixXd domain=(box(.25)*rotation.transpose()).eval();domain.rowwise()+=shift.transpose();
        const Eigen::Vector3d axis=7*rotation*Eigen::Vector3d::UnitZ();
        spin::detail::ProjectedDomainClearance check(domain,axis,.001);
        const auto transform=[&](std::vector<Eigen::Vector3d> v) {
            for(auto& x:v)x=rotation*x+shift;return v;
        };
        require(check.graphAbove(transform(sloped),faces,rim),"remote low vertex rejected a locally distant graph");
        auto bad=sloped;for(auto& x:bad)x.z()-=20;
        require(!check.graphAbove(transform(bad),faces,rim),"wrong-side far graph was accepted");
        bad=sloped;for(auto& x:bad)x.z()-=10;
        require(!check.graphAbove(transform(bad),faces,rim),"far graph crossing the design was accepted");
        bad=sloped;for(auto& x:bad)x.z()=.2505;
        require(!check.graphAbove(transform(bad),faces,rim),"sub-margin local height was accepted");
        for(auto& x:bad)x.z()=.252;
        require(check.graphAbove(transform(bad),faces,rim),"local clearance depends on axis length or orientation");
    }
    spin::detail::ProjectedDomainClearance check(box(.25),Eigen::Vector3d::UnitZ(),1e-10);
    const std::vector<Eigen::Vector3d> square={{-1,-1,3},{1,-1,3},{1,1,3},{-1,1,3}};
    const std::vector<int> outer={0,1,2,3};
    require(check.graphAbove(square,{{0,1,2},{0,2,3}},outer),"coherent far graph rejected");
    require(!check.graphAbove(square,{{0,1,2}},outer),"incomplete far graph coverage accepted");
    require(!check.graphAbove(square,{{0,1,2},{0,3,2}},outer),"folded projected graph accepted");
    require(!check.graphAbove(square,{{0,1,2},{0,2,3},{0,1,2}},outer),"duplicate far face accepted");
    require(!check.graphAbove(square,{{0,1,9}},outer),"out-of-range far graph index accepted");
    require(!check.graphAbove(square,{},outer),"empty far graph accepted");
}
spin::QuadricCoefficients10 retract(const spin::QuadricCoefficients10& q,
    const Eigen::VectorXd& direction,double step) {
    auto p=(q.values()+step*direction).eval();
    Eigen::Matrix3d A;
    A<<p(4),p(9),p(8),p(9),p(5),p(7),p(8),p(7),p(6);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> e(A);
    int zero=0;e.eigenvalues().cwiseAbs().minCoeff(&zero);
    auto d=e.eigenvalues().eval();d(zero)=0;
    A=e.eigenvectors()*d.asDiagonal()*e.eigenvectors().transpose();
    p(4)=A(0,0);p(5)=A(1,1);p(6)=A(2,2);
    p(7)=A(1,2);p(8)=A(0,2);p(9)=A(0,1);
    return spin::QuadricCoefficients10::fromUnnormalized(p);
}
void geometry(bool saddle,bool inverted) {
    using Q=spin::QuadricCoefficients10;
    Q::Vector p;p<<.125,.1,-.08,-1,1,(saddle?-2:2),0,0,0,0;
    const Q q=Q::fromUnnormalized(inverted?(-p).eval():p);
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions=2;settings.truncationRadius=16;
    settings.maximumRadiusDoublings=0;
    const auto legacy=spin::triangulateQuadric(q,box(.2),settings);
    enable(settings);
    const auto t=spin::triangulateQuadric(q,box(.2),settings);
    require(t.hash!=legacy.hash,"new policy silently executed the legacy grid");
    require(t.effectiveIcosphereSubdivisions==2,"fixed normal resolution changed");
    require(t.sourceNormalVertexByVertex.size()==size_t(t.closedBoundary.vertices.rows()),"normal identities missing");
    require(t.fixedNormals.rows()>0 && t.interfaceVertexCount>3,"normal patch missing");
    const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
        t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
    require(audit.facts().boundaryValid,"fixed-normal boundary failed exact audit");
    for(int i=0;i<t.interfaceVertexCount;++i) {
        const int id=t.sourceNormalVertexByVertex[i];
        require(id>=0 && id<t.fixedNormals.rows(),"invalid finite normal identity");
        const Eigen::Vector3d normal=t.fixedNormals.row(id),x=t.closedBoundary.vertices.row(i);
        // Hand-derived map for this fixture; independent of production map.
        const double k=-1/normal.z();
        const double a=(k*normal.x()-.1)/2;
        const double b=(k*normal.y()+.08)/(saddle?-4:4);
        const double z=.125+.1*a-.08*b+a*a+(saddle?-2:2)*b*b;
        require((x-Eigen::Vector3d(a,b,z)).norm()<1e-10,"finite vertex is not the analytic fixed-normal point");
        require(std::abs(q.evaluate(x))<1e-10,"source vertex left paraboloid");
    }
    for(int f=t.interfaceFaceCount;f<t.closedBoundary.faces.rows();++f) {
        bool separated=false;
        for(int a=0;a<3;++a) {
            bool above=true,below=true;
            for(int k=0;k<3;++k) {
                double x=t.closedBoundary.vertices(t.closedBoundary.faces(f,k),a);
                above&=x>.2;below&=x<-.2;
            }
            separated|=above||below;
        }
        require(separated,"closure changes material within the design AABB");
    }
    require(spin::triangulateQuadric(q,box(.2),settings).hash==t.hash,"source is nondeterministic");
    bool rejected=false;
    try {(void)spin::triangulateQuadric(q,box(100),settings);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"insufficient finite patch coverage accepted");
    for(int bad=0;bad<3;++bad) {
        auto invalid=settings;
        if(bad==0)invalid.algebraicFiniteMap=true;
        if(bad==1)invalid.paperProjectiveHyperboloid=true;
        if(bad==2)invalid.maximumIcosphereSubdivisions=3;
        rejected=false;
        try{(void)spin::triangulateQuadric(q,box(.2),invalid);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"incompatible fixed-normal profile silently accepted");
    }
    const auto full=spin::evaluateQuadricSourceMotion(q,box(.2),settings);
    require(!full.defined,"reduced finite patch was promoted to full coefficient source motion");
    const auto motion=spin::evaluateQuadricStratumSourceMotion(q,box(.2),settings);
    require(motion.defined,"fixed-normal reduced motion unavailable: "+motion.diagnostic);
    require(motion.triangulation.hash==t.hash,"motion constructed a different source");
    require(motion.coefficientBasis.rows()==10 && motion.coefficientBasis.cols()==8,
            "fixed-normal source substituted a full-dimensional derivative");
    double firstError=0,lastError=0;
    for(double h:{1e-5,3e-6,1e-6}) {
        double e2=0,d2=0;
        for(int j=0;j<8;++j) {
            const auto plus=spin::triangulateQuadric(retract(q,motion.coefficientBasis.col(j),h),box(.2),settings);
            const auto minus=spin::triangulateQuadric(retract(q,motion.coefficientBasis.col(j),-h),box(.2),settings);
            require(plus.closedBoundary.faces==t.closedBoundary.faces && minus.closedBoundary.faces==t.closedBoundary.faces &&
                    plus.sourceNormalVertexByVertex==t.sourceNormalVertexByVertex &&
                    minus.sourceNormalVertexByVertex==t.sourceNormalVertexByVertex &&
                    plus.truncationRadius==t.truncationRadius && minus.truncationRadius==t.truncationRadius,
                    "source topology or cap branch changed in differential stencil");
            const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
            for(int v=0;v<fd.rows();++v) {
                const Eigen::Vector3d d=motion.vertexVelocity.at(v).col(j);
                const double error=(fd.row(v).transpose()-d).norm();
                require(d.norm()<1e-8?error<1e-7:error/d.norm()<1e-6,
                        "per-vertex fixed-normal derivative failed");
                e2+=error*error;d2+=d.squaredNorm();
            }
        }
        const double relative=std::sqrt(e2/d2);
        require(relative<1e-6,"complete eight-direction source derivative failed");
        if(h==1e-5)firstError=relative;
        lastError=relative;
        std::cout<<"saddle="<<saddle<<" inverted="<<inverted<<" h="<<h<<" relative="<<relative<<'\n';
    }
    require(firstError<1e-8 || lastError<firstError*.3,"source finite differences do not refine");
    std::cout<<"saddle="<<saddle<<" inverted="<<inverted<<" vertices="
             <<t.closedBoundary.vertices.rows()<<" source="<<t.hash<<" motion="<<motion.scientificHash<<'\n';
}
void distantClosure() {
    using Q=spin::QuadricCoefficients10;
    Q::Vector p;p<<1000,.1,-.08,-1,1,2,0,0,0,0;
    const Q q=Q::fromUnnormalized(p);
    spin::QuadricTriangulationSettings s;
    enable(s);s.icosphereSubdivisions=2;s.truncationRadius=16;s.maximumRadiusDoublings=0;
    bool rejected=false;
    try{(void)spin::triangulateQuadric(q,box(.2),s);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"distant positive domain omitted by wrong-side disjoint closure");
    s.maximumRadiusDoublings=6;
    const auto t=spin::triangulateQuadric(q,box(.2),s);
    require(t.truncationRadius>16,"distant closure did not use its bounded radius policy");
    // Independent solid-angle containment. q is strictly positive on this box.
    for(int i=0;i<8;++i) {
        const Eigen::Vector3d x=box(.2).row(i);
        double angle=0;
        for(int f=0;f<t.closedBoundary.faces.rows();++f) {
            const Eigen::Vector3d a=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,0)).transpose()-x;
            const Eigen::Vector3d b=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,1)).transpose()-x;
            const Eigen::Vector3d c=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,2)).transpose()-x;
            angle+=2*std::atan2(a.dot(b.cross(c)),a.norm()*b.norm()*c.norm()+a.dot(b)*c.norm()+b.dot(c)*a.norm()+c.dot(a)*b.norm());
        }
        require(std::abs(angle/(4*std::acos(-1.0))-1)<1e-8,"bounded closure omits positive design material");
    }
    Q::Vector plane=Q::Vector::Zero();plane(3)=1;
    rejected=false;
    try{(void)spin::triangulateQuadric(Q::fromUnnormalized(plane),box(.2),s);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"non-paraboloid accepted by fixed-normal profile");
}
void savedLocalCap() {
    // Frozen seed-1's first rejected trial, not a new optimizer input.
    using Q=spin::QuadricCoefficients10;
    Q::Vector p;p<<.1571303233553719,.576412573974079,.1306184650714625,-.6322026023925956,
        0,0,0,.4740922569529424,-.03941857499305815,0;
    p.normalize();
    spin::QuadricTriangulationSettings s;enable(s);
    s.icosphereSubdivisions=2;s.truncationRadius=4;s.maximumRadiusDoublings=4;
    spin::QuadricTriangulation t;
    try {t=spin::triangulateQuadric(Q::fromUnnormalized(p),box(.25),s);}
    catch(const std::exception& e) {
        throw std::runtime_error(std::string("locally distant frozen cap was rejected: ")+e.what());
    }
    require(t.truncationRadius<=64,"saved cap bypassed its original radius budget");
    require(verification::auditExactTriangleBoundary(t.closedBoundary.vertices,t.closedBoundary.faces,
        verification::ExactBoundaryOrientation::OutwardPositive).facts().boundaryValid,
        "locally certified cap produces invalid boundary");
}
void savedSaddleWinding() {
    using Q=spin::QuadricCoefficients10;
    // First boundary-loop failure of frozen seed 0, evaluation 3. This is a
    // source-geometry regression, not a new optimizer configuration.
    Q::Vector p;
    p<<.22723201164453788,.5754976591179629,-.07095486715854094,.754327523868548,
        0,0,0,-.1836033899756779,.09702188404694938,0;
    p.normalize();
    const auto q=Q::fromUnnormalized(p);
    spin::QuadricTriangulationSettings s;
    enable(s);s.icosphereSubdivisions=2;s.truncationRadius=16;s.maximumRadiusDoublings=6;
    // Independent rational clipping of all 134 projected triangles found
    // disjoint interiors and covered area exactly equal to the design hull.
    // The prior halfspace-kernel rejection was a conservative false negative,
    // not evidence that this frozen domain lies outside the finite patch.
    const auto original=spin::triangulateQuadric(q,box(.25),s);
    require(verification::auditExactTriangleBoundary(original.closedBoundary.vertices,
        original.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive).facts().boundaryValid,
        "covered frozen domain has invalid source boundary");

    // A separate local geometry control: the uniquely determined transverse
    // stationary point and a fixed .02 half-width, chosen before evaluation.
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(q.quadraticMatrix());
    int zero=0;eigen.eigenvalues().cwiseAbs().minCoeff(&zero);
    Eigen::Vector3d axis=eigen.eigenvectors().col(zero),center=Eigen::Vector3d::Zero();
    if(axis.dot(q.linearVector())<0)axis=-axis;
    for(int j=0;j<3;++j)if(j!=zero) {
        const Eigen::Vector3d d=eigen.eigenvectors().col(j);
        center-=d*q.linearVector().dot(d)/(2*eigen.eigenvalues()(j));
    }
    center-=axis*q.evaluate(center)/axis.dot(q.linearVector());
    const Eigen::MatrixXd domain=(box(.02).rowwise()+center.transpose()).eval();
    const auto t=spin::triangulateQuadric(q,domain,s);
    const auto audit=verification::auditExactTriangleBoundary(t.closedBoundary.vertices,
        t.closedBoundary.faces,verification::ExactBoundaryOrientation::OutwardPositive);
    require(audit.facts().boundaryValid,"saved saddle local control has an invalid closed boundary");
    int centroidDisagreements=0;
    for(int f=0;f<t.interfaceFaceCount;++f) {
        const Eigen::Vector3d a=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,0));
        const Eigen::Vector3d b=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,1));
        const Eigen::Vector3d c=t.closedBoundary.vertices.row(t.closedBoundary.faces(f,2));
        const Eigen::Vector3d normal=(b-a).cross(c-a);
        require(normal.dot(axis)<0,"finite source is not a coherently oriented height graph");
        centroidDisagreements+=normal.dot(-q.gradient((a+b+c)/3))<0;
    }
    require(centroidDisagreements==2,"saved regression no longer exercises centroid-gradient disagreement");
    require(spin::triangulateQuadric(q,domain,s).hash==t.hash,"saved source winding is nondeterministic");
    const auto motion=spin::evaluateQuadricStratumSourceMotion(q,domain,s);
    require(motion.defined,"saved source derivative unavailable: "+motion.diagnostic);
    require(motion.coefficientBasis.cols()==8,"saved source derivative lost its true reduced dimension");
    bool differentialPassed=true;
    double coarseError=0,fineError=0;
    // Retain the coarse stencil as a truncation-error diagnostic. Three finer
    // stencils must each meet the unchanged 1e-6 derivative gate; no solver
    // tolerance, normal resolution or geometry parameter changes here.
    for(double h:{1e-5,3e-6,1e-6,3e-7}) {
        double error2=0,norm2=0;
        for(int j=0;j<8;++j) {
            const auto plus=spin::triangulateQuadric(retract(q,motion.coefficientBasis.col(j),h),domain,s);
            const auto minus=spin::triangulateQuadric(retract(q,motion.coefficientBasis.col(j),-h),domain,s);
            require(plus.closedBoundary.faces==t.closedBoundary.faces && minus.closedBoundary.faces==t.closedBoundary.faces &&
                    plus.sourceNormalVertexByVertex==t.sourceNormalVertexByVertex &&
                    minus.sourceNormalVertexByVertex==t.sourceNormalVertexByVertex &&
                    plus.truncationRadius==t.truncationRadius && minus.truncationRadius==t.truncationRadius,
                    "saved source derivative stencil changed branch");
            const Eigen::MatrixXd fd=(plus.closedBoundary.vertices-minus.closedBoundary.vertices)/(2*h);
            for(int v=0;v<fd.rows();++v) {
                const Eigen::Vector3d d=motion.vertexVelocity.at(v).col(j);
                error2+=(fd.row(v).transpose()-d).squaredNorm();norm2+=d.squaredNorm();
            }
        }
        const double relative=std::sqrt(error2/norm2);
        if(h==1e-5)coarseError=relative;
        else differentialPassed&=relative<1e-6;
        fineError=relative;
        std::cout<<"saved_saddle h="<<h<<" relative="<<relative<<'\n';
    }
    require(differentialPassed,"saved source eight-direction differential failed");
    require(fineError<coarseError*.05,"saved source differential does not refine with step size");
}
}
int main() {
    try {savedLocalCap();localGraphCertificate();projectedCoverage();savedSaddleWinding();distantClosure();for(bool saddle:{false,true})for(bool inverted:{false,true})geometry(saddle,inverted);}
    catch(const std::exception& error){std::cerr<<"[FAIL] "<<error.what()<<'\n';return 1;}
    std::cout<<"[PASS] bounded fixed-normal paraboloid geometry and eight-direction source motion\n";
}
