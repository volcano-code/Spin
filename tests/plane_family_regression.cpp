#include "spin/PlaneMomentEvaluator.h"
#include "quadric_certified_empty_cases.h"
#include <array>
#include <cmath>
#include <iomanip>

namespace {
using P=Eigen::Vector4d;
using J=spin::PlaneMomentEvaluation::Jacobian;
using certified_empty_tests::require;
using certified_empty_tests::mustReject;

const int powers[10][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{2,0,0},
                        {0,2,0},{0,0,2},{0,1,1},{1,0,1},{1,1,0}};
double integral(int k,double hi) {
    return (std::pow(hi,k+1)-std::pow(-1.,k+1))/(k+1);
}
void basic() {
    const auto r=certified_empty_tests::request();
    const P p(.2,0,0,-1);
    const auto e=spin::evaluatePlaneMaterialMoments(r,p);
    require(e.coefficientJacobian.has_value(),"plane-family 10x4 Jacobian unavailable");
    require(e.variableMaterialBoundaryTarget.faces.rows()>0,
        "plane-family result omitted the actual Boolean material boundary");
    const auto independent=physics::RawMassMoments10::fromMassMoments(
        computeSurfaceMassMoments(e.variableMaterialBoundaryTarget.vertices,
                                  e.variableMaterialBoundaryTarget.faces,r.materialDensity));
    require((independent.values()-e.variableMomentsTarget.values()).norm()<1e-10,
        "independent output-mesh integration disagrees with the plane evaluator moments");
    J exact;
    for(int m=0;m<10;++m) {
        const int i=powers[m][0],j=powers[m][1],k=powers[m][2];
        const double moment=1.7*integral(i,1)*integral(j,1)*integral(k,.2);
        require(std::abs(e.variableMomentsTarget.values()(m)-moment)<1e-10,"box slice analytic moment mismatch");
        exact(m,0)=1.7*integral(i,1)*integral(j,1)*std::pow(.2,k);
        exact(m,1)=1.7*integral(i+1,1)*integral(j,1)*std::pow(.2,k);
        exact(m,2)=1.7*integral(i,1)*integral(j+1,1)*std::pow(.2,k);
        exact(m,3)=1.7*integral(i,1)*integral(j,1)*std::pow(.2,k+1);
    }
    const double error=(*e.coefficientJacobian-exact).norm()/exact.norm();
    require(error<1e-10,"analytic 10x4 slice derivative mismatch");
    require((e.totalMomentsTarget.values()-e.variableMomentsTarget.values()-r.fixedShellMoments.values()).norm()<1e-10,"fixed moments lost");
    require(e.interfaceTriangleCount>0 && !e.capSeparationCertificate.empty() && e.scientificHash.size()==64,"missing plane evidence");
    std::cout<<"[PASS] analytic slice all10 moments and 40 derivatives error="<<error<<'\n';
}

void invariance() {
    auto r=certified_empty_tests::request();
    const P p(.17,.37,-.53,.91);
    const auto e=spin::evaluatePlaneMaterialMoments(r,p);
    const auto repeat=spin::evaluatePlaneMaterialMoments(r,p);
    require(e.coefficientJacobian.has_value(),"oblique derivative unavailable");
    require(e.scientificHash==repeat.scientificHash,"plane identity nondeterministic");
    for(double scale:{.125,3.,16.}) {
        const auto s=spin::evaluatePlaneMaterialMoments(r,scale*p);
        require(s.coefficientJacobian.has_value(),"scaled derivative unavailable");
        require((s.variableMomentsTarget.values()-e.variableMomentsTarget.values()).norm()<1e-10,"positive scale changes shape");
        require((scale**s.coefficientJacobian-*e.coefficientJacobian).norm()<1e-9,"missing inverse-scale derivative chain");
        require(s.scientificHash!=e.scientificHash,"input scale not bound to derivative identity");
    }
    require((*e.coefficientJacobian*p).norm()<1e-10,"radial scale-null direction violated");
    const auto neg=spin::evaluatePlaneMaterialMoments(r,-p);
    const auto full=spin::evaluatePlaneMaterialMoments(r,P(10,0,0,-1));
    const auto empty=spin::evaluatePlaneMaterialMoments(r,P(-10,0,0,-1));
    require(neg.coefficientJacobian && full.coefficientJacobian && empty.coefficientJacobian,"complement/full/empty unavailable");
    require((e.variableMomentsTarget.values()+neg.variableMomentsTarget.values()-full.variableMomentsTarget.values()).norm()<1e-9,"complement moments violated");
    require((*e.coefficientJacobian-*neg.coefficientJacobian).norm()<1e-9,"complement derivative sign wrong");
    require(full.coefficientJacobian->isZero(0) && empty.coefficientJacobian->isZero(0),"full/empty derivative not exactly zero");
    require(empty.variableMomentsTarget.values().isZero(0),"empty value fabricated");
    r.materialDensity*=2;
    const auto den=spin::evaluatePlaneMaterialMoments(r,p);
    require(den.coefficientJacobian && (*den.coefficientJacobian-2**e.coefficientJacobian).norm()<1e-9,"density factor wrong");
    require(den.scientificHash!=e.scientificHash,"density missing from identity");
    std::cout<<"[PASS] scaling normalization radial-null complement empty full density hashes\n";
}

void differences() {
    const auto r=certified_empty_tests::request();
    const std::array<P,6> cases={{P(.2,0,0,-1),P(-.31,1,0,0),P(.23,0,-1,0),
        P(.17,.37,-.53,.91),P(-.27,-.61,.29,.83),P(.13,.2,.3,1)}};
    for(std::size_t c=0;c<cases.size();++c) {
        const P p=cases[c];const auto e=spin::evaluatePlaneMaterialMoments(r,p);
        require(e.coefficientJacobian.has_value(),"FD center derivative unavailable");
        double first=0,last=0;
        for(double h:{1e-3,3e-4,1e-4,3e-5}) {
            J fd;
            for(int i=0;i<4;++i) {
                P a=p,b=p;a(i)+=h;b(i)-=h;
                fd.col(i)=(spin::evaluatePlaneMaterialMoments(r,a).variableMomentsTarget.values()-
                    spin::evaluatePlaneMaterialMoments(r,b).variableMomentsTarget.values())/(2*h);
            }
            const double error=(fd-*e.coefficientJacobian).norm()/std::max(1.,e.coefficientJacobian->norm());
            if(h==1e-3)first=error;last=error;
            std::cout<<"FD case="<<c<<" h="<<h<<" relative="<<error<<'\n';
        }
        require(last<1e-6 && (first<1e-9 || last<first*.1),"multi-step finite difference did not converge");
        const P n=p.normalized();P u(1,.3,-.2,.1);u-=n*n.dot(u);u.normalize();
        const auto center=spin::evaluatePlaneMaterialMoments(r,n);
        const double h=1e-4;
        const auto a=spin::evaluatePlaneMaterialMoments(r,std::cos(h)*n+std::sin(h)*u);
        const auto b=spin::evaluatePlaneMaterialMoments(r,std::cos(h)*n-std::sin(h)*u);
        const auto fd=(a.variableMomentsTarget.values()-b.variableMomentsTarget.values())/(2*h);
        const auto jvp=*center.coefficientJacobian*u;
        require((fd-jvp).norm()/std::max(1.,jvp.norm())<1e-6,"tangent JVP chain mismatch");
    }
    std::cout<<"[PASS] six planes multi-step finite differences and tangent JVP\n";
}

void invalid() {
    const auto r=certified_empty_tests::request();
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P::Zero());},"zero polynomial accepted");
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P(1,0,0,0));},"zero normal accepted");
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P(NAN,0,0,1));},"NaN accepted");
    auto bad=r;bad.materialDensity=-1;
    mustReject([&]{spin::evaluatePlaneMaterialMoments(bad,P(.2,0,0,-1));},"negative density accepted");
    bad=r;bad.designDomain.faces.conservativeResize(11,3);
    mustReject([&]{spin::evaluatePlaneMaterialMoments(bad,P(.2,0,0,-1));},"open mesh accepted");
    spin::PlaneMomentEvaluatorSettings s;s.truncationRadius=.5;
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P(.2,0,0,-1),s);},"cap interference accepted as plane-family region");
    s.truncationRadius=-1;
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P(.2,0,0,-1),s);},"negative truncation radius accepted");
    s.truncationRadius=1e308;
    mustReject([&]{spin::evaluatePlaneMaterialMoments(r,P(.2,0,0,-1),s);},"overflowing finite truncation radius accepted");
    for(double t:{1.,1.-1e-12}) {
        bool rejected=false;
        try { rejected=!spin::evaluatePlaneMaterialMoments(r,P(t,0,0,-1)).coefficientJacobian.has_value(); }
        catch(const std::runtime_error&){rejected=true;}
        require(rejected,"contact/near-contact granted derivative");
    }
    for(double bx:{0.0}) {
        const auto nearAxis=spin::evaluatePlaneMaterialMoments(r,P(.2,bx,0,1));
        require(nearAxis.coefficientJacobian.has_value(),
            "stable axis-adjacent plane chart withheld its derivative");
        const double h=1e-5;
        Eigen::Matrix<double,10,4> fd;
        const P coefficients(.2,bx,0,1);
        for(int column=0;column<4;++column) {
            P plus=coefficients,minus=coefficients;
            plus(column)+=h;minus(column)-=h;
            fd.col(column)=(spin::evaluatePlaneMaterialMoments(r,plus).variableMomentsTarget.values()-
                spin::evaluatePlaneMaterialMoments(r,minus).variableMomentsTarget.values())/(2*h);
        }
        require((fd-*nearAxis.coefficientJacobian).norm()/std::max(1.,fd.norm())<1e-6,
            "axis-adjacent plane differential failed independent finite differences");
    }
    bool chartSeamRejected=false;
    try {
        chartSeamRejected=!spin::evaluatePlaneMaterialMoments(
            r,P(.2,.9,std::sqrt(.19),0)).coefficientJacobian.has_value();
    } catch(const std::runtime_error&) { chartSeamRejected=true; }
    require(chartSeamRejected,"plane tangent-frame chart seam granted a derivative");
    std::cout<<"[PASS] invalid inputs cap interference contact and near-contact fail closed\n";
}

void largeTruncation() {
    const auto r=certified_empty_tests::request();
    const P p(1.8100001,.37,-.53,.91);
    spin::QuadricCoefficients10::Vector q=spin::QuadricCoefficients10::Vector::Zero();q.head<4>()=p;
    int nativeCuts=0;
    for(double radius:{1e8,1e10,1e12,1e14,1e16}) {
        spin::QuadricMomentEvaluatorSettings native;native.triangulation.truncationRadius=radius;
        native.emptyRegionPolicy=spin::QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1;
        try {
            const auto n=spin::evaluateQuadricMaterialMoments(r,spin::QuadricCoefficients10::fromUnnormalized(q),native);
            if(n.interfaceTriangleCount>0) {
                ++nativeCuts;bool safelyRejected=false;
                try {
                    spin::PlaneMomentEvaluatorSettings s;s.truncationRadius=radius;
                    const auto e=spin::evaluatePlaneMaterialMoments(r,p,s);
                    safelyRejected=!e.coefficientJacobian || !e.coefficientJacobian->isZero(0);
                } catch(const std::runtime_error&){safelyRejected=true;}
                require(safelyRejected,"ideal full sign granted zero derivative despite native partial cut");
            }
        } catch(const std::runtime_error& error) {
            if(std::string(error.what()).find("ideal full sign")!=std::string::npos)throw;
        }
    }
    require(nativeCuts>0,"large-radius adversarial fixture did not exercise native partial cut");
    std::cout<<"[PASS] large-radius ideal/native full-state distinction cases="<<nativeCuts<<'\n';
}

void transformedFrame() {
    const auto r=certified_empty_tests::request();auto moved=r;
    const Eigen::Matrix3d R=Eigen::AngleAxisd(.4,Eigen::Vector3d(1,2,3).normalized()).toRotationMatrix();
    const Eigen::Vector3d t(3,-2,1);
    for(auto* m:{&moved.designDomain,&moved.outerSurface})
        for(int i=0;i<m->vertices.rows();++i)m->vertices.row(i)=(R*m->vertices.row(i).transpose()+t).transpose();
    moved.targetFrame.axes=R;moved.targetFrame.origin=t;
    MassMoments f;f.mass=95.2;f.firstMoment=95.2*t;
    f.secondMoment=Eigen::Matrix3d::Identity()*(248./3.*1.7)+95.2*t*t.transpose();
    moved.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(f);
    const P p(.17,.37,-.53,.91);
    const auto a=spin::evaluatePlaneMaterialMoments(r,p),b=spin::evaluatePlaneMaterialMoments(moved,p);
    require(a.coefficientJacobian && b.coefficientJacobian,"transformed derivative unavailable");
    require((a.totalMomentsTarget.values()-b.totalMomentsTarget.values()).norm()<1e-9,"target-frame moments incorrect");
    require((*a.coefficientJacobian-*b.coefficientJacobian).norm()<1e-9,"target-frame derivative incorrect");
    require(a.scientificHash!=b.scientificHash,"request frame not hash-bound");
    std::cout<<"[PASS] nonidentity target frame and fixed-shell transform\n";
}

void otherDomains() {
    auto r=certified_empty_tests::request();
    r.designDomain.vertices.resize(4,3);r.designDomain.vertices<<0,0,0,1,0,0,0,1,0,0,0,1;
    r.designDomain.faces.resize(4,3);r.designDomain.faces<<0,2,1,0,1,3,1,2,3,2,0,3;
    const auto tetra=spin::evaluatePlaneMaterialMoments(r,P(.6,-1,-1,-1));
    require(tetra.coefficientJacobian.has_value(),"tetrahedral derivative unavailable");
    for(int m=0;m<10;++m) {
        const int i=powers[m][0],j=powers[m][1],k=powers[m][2],s=i+j+k;
        const double exact=1.7*std::tgamma(i+1)*std::tgamma(j+1)*std::tgamma(k+1)*std::pow(.6,s+3)/std::tgamma(s+4);
        require(std::abs(tetra.variableMomentsTarget.values()(m)-exact)<1e-10,"analytic tetra moment mismatch");
        Eigen::RowVector4d derivative;derivative<<(s+3)*exact/.6,(i+1)*exact,(j+1)*exact,(k+1)*exact;
        require((tetra.coefficientJacobian->row(m)-derivative).norm()<1e-10,"analytic tetra derivative mismatch");
    }
    r=certified_empty_tests::request();
    const auto a=certified_empty_tests::box(-1,-.2),b=certified_empty_tests::box(.2,1);
    r.designDomain.vertices.resize(16,3);r.designDomain.vertices<<a.vertices,b.vertices;
    r.designDomain.faces.resize(24,3);r.designDomain.faces<<a.faces,(b.faces.array()+8).matrix();
    const P p(.07,.3,-.2,.9);const auto e=spin::evaluatePlaneMaterialMoments(r,p);
    require(e.coefficientJacobian.has_value(),"disconnected closed-domain derivative unavailable");
    J fd;const double h=1e-5;
    for(int i=0;i<4;++i){P x=p,y=p;x(i)+=h;y(i)-=h;
        fd.col(i)=(spin::evaluatePlaneMaterialMoments(r,x).variableMomentsTarget.values()-
                   spin::evaluatePlaneMaterialMoments(r,y).variableMomentsTarget.values())/(2*h);}
    require((fd-*e.coefficientJacobian).norm()/std::max(1.,fd.norm())<1e-6,"disconnected-domain finite difference failed");
    std::cout<<"[PASS] analytic tetrahedron and disconnected polyhedral domain\n";
}
}

int main() {
    try {std::cout<<std::setprecision(17);basic();invariance();differences();invalid();largeTruncation();transformedFrame();otherDomains();return 0;}
    catch(const std::exception& e){std::cerr<<"[FAIL] plane-family: "<<e.what()<<'\n';return 1;}
}
