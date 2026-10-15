#include "spin2014/SourceBinaryMoments.h"
#include "spin2014/Spin2014GivensProblem.h"
#include "manufacturing/SlicerValidation.h"
#include <Eigen/Eigenvalues>
#include <functional>
#include <sstream>
#include <stdexcept>

namespace spin2014 {
SourceBinaryCuboids sourceBinaryCuboids(const SourceAdaptiveDomain& domain,const Eigen::VectorXd& beta,std::size_t cap) {
    if(cap==0||beta.size()!=static_cast<Eigen::Index>(domain.leaves().size())||!beta.allFinite()||((beta.array()!=0)&&(beta.array()!=1)).any())throw std::invalid_argument("binary source values/cuboid cap invalid");
    SourceBinaryCuboids r;r.sourceHash=domain.source().partitionSha256();r.partitionHash=domain.partitionHash();r.maximumCuboids=cap;
    const auto& source=domain.source();const auto& cfg=source.config();
    std::function<void(SourcePartitionKey)> visit=[&](SourcePartitionKey key){auto q=source.cell(key.level,key.coordinates);if(q.moments.mass()==0)return;
        if(q.kind==SourceSurfaceCellKind::Full){if(r.cuboids.size()>=cap)throw std::length_error("RESOURCE_LIMIT: binary output cuboids");
            SourceRemovedCuboid b;b.cell=key;
            const auto canonical=source.cellBounds(key.level,key.coordinates);
            b.lower=canonical.lower;b.upper=canonical.upper;
            if(!b.lower.allFinite()||!b.upper.allFinite()||(b.upper.array()<=b.lower.array()).any())throw std::overflow_error("invalid binary cuboid arithmetic");r.cuboids.push_back(b);return;}
        if(key.level>=cfg.maximumDepth)throw std::logic_error("nonfull terminal has positive carveable contribution");
        for(int k=0;k<8;++k){auto child=key;++child.level;for(int a=0;a<3;++a)child.coordinates[a]=2*key.coordinates[a]+((k>>a)&1);visit(child);}
    };
    for(std::size_t i=0;i<domain.leaves().size();++i)if(beta(i)==1)visit(domain.leaves()[i]);
    std::sort(r.cuboids.begin(),r.cuboids.end(),[](const auto&a,const auto&b){return a.cell<b.cell;});
    std::ostringstream o;o.imbue(std::locale::classic());o<<std::hexfloat<<"source-binary-cuboids-v1 "<<r.sourceHash<<' '<<r.partitionHash<<' '<<cap<<' ';
    for(const auto& b:r.cuboids){o<<b.cell.level<<':'<<b.cell.coordinates[0]<<','<<b.cell.coordinates[1]<<','<<b.cell.coordinates[2]<<';';for(int a=0;a<3;++a)o<<b.lower(a)<<','<<b.upper(a)<<';';}
    r.outputHash=manufacturing::computeContentSha256(o.str());return r;
}
SourceBinaryMomentCheck checkSourceBinaryMoments(const SourceSurfaceOctree& source,const SourceBinaryCuboids& output,
    const physics::RawMassMoments10& expected,Spin2014Objective objective,double phi,double density,double tolerance,
    const spin::SpinTargetFrame& targetFrame) {
    return checkSourceBinaryMoments(source,output,expected,objective,phi,density,density,tolerance,targetFrame);
}
SourceBinaryMomentCheck checkSourceBinaryMoments(const SourceSurfaceOctree& source,const SourceBinaryCuboids& output,
    const physics::RawMassMoments10& expected,Spin2014Objective objective,double phi,
    double boundaryDensity,double interiorDensity,double tolerance,const spin::SpinTargetFrame& targetFrame) {
    SourceBinaryMomentCheck r;
    if(!std::isfinite(boundaryDensity)||boundaryDensity<=0||!std::isfinite(interiorDensity)||interiorDensity<=0||
       !std::isfinite(tolerance)||tolerance<=0||!std::isfinite(phi))throw std::invalid_argument("invalid binary check controls");
    try {
        // This geometric integration intentionally uses only the original outer
        // and design triangles and emitted box coordinates, never oracle moment queries.
        MassMoments total=computeSurfaceMassMoments(source.outerVertices(),source.outerFaces(),boundaryDensity)
            -computeSurfaceMassMoments(source.designVertices(),source.designFaces(),boundaryDensity)
            +computeSurfaceMassMoments(source.designVertices(),source.designFaces(),interiorDensity);
        Eigen::MatrixXd vertices(8,3);Eigen::MatrixXi faces(12,3);
        faces<<0,2,1,1,2,3,4,5,6,5,7,6,0,1,4,1,5,4,2,6,3,3,6,7,0,4,2,2,4,6,1,3,5,3,7,5;
        for(const auto& b:output.cuboids){
            if(!b.lower.allFinite()||!b.upper.allFinite()||(b.upper.array()<=b.lower.array()).any())throw std::invalid_argument("invalid emitted cuboid");
            for(int i=0;i<8;++i)for(int a=0;a<3;++a)vertices(i,a)=(i&(1<<a))?b.upper(a):b.lower(a);
            total-=computeSurfaceMassMoments(vertices,faces,interiorDensity);
        }
        if(!total.isFinite()||!expected.isFinite())throw std::invalid_argument("nonfinite binary geometric moments");
        r.independentMoments=physics::RawMassMoments10::fromMassMoments(total);
        const double densityScale=std::max(boundaryDensity,interiorDensity);
        const double side=source.config().rootSide;double scale[3]={densityScale*std::pow(side,3),densityScale*std::pow(side,4),densityScale*std::pow(side,5)};
        for(double s:scale)if(!std::isfinite(s)||s<=0)throw std::invalid_argument("invalid fixed root scaling");
        r.normalizedMomentError=0;for(int i=0;i<10;++i)r.normalizedMomentError=std::max(r.normalizedMomentError,std::abs(r.independentMoments->values()(i)-expected.values()(i))/scale[i==0?0:(i<=3?1:2)]);
        if(total.mass<=0)throw std::invalid_argument("nonpositive binary mass");
        const Eigen::Matrix3d inertia=total.inertiaAboutCenter();if(!inertia.allFinite())throw std::invalid_argument("nonfinite central inertia");
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(inertia);if(eigen.info()!=Eigen::Success)throw std::invalid_argument("central inertia eigensolve failed");
        r.minimumCentralInertiaEigenvalue=eigen.eigenvalues().minCoeff();
        const auto targetMoments=r.independentMoments->transformed(
            targetFrame.axes.transpose(),-targetFrame.axes.transpose()*targetFrame.origin
        );
        Spin2014GivensProblem::Vector11 state;state.head<10>()=targetMoments.values();state(10)=phi;
        if(objective==Spin2014Objective::StaticLowering&&phi!=0)throw std::invalid_argument("static binary phi must be zero");
        r.constraints=Spin2014GivensProblem(objective).evaluate(state).constraints;
        r.rawConstraintResidual=r.constraints.lpNorm<Eigen::Infinity>();r.normalizedConstraintResidual=0;
        const int firstCount=objective==Spin2014Objective::YoYo?3:2;
        for(int i=0;i<r.constraints.size();++i)r.normalizedConstraintResidual=std::max(r.normalizedConstraintResidual,std::abs(r.constraints(i))/scale[i<firstCount?1:2]);
        r.passed=std::isfinite(r.normalizedMomentError)&&std::isfinite(r.normalizedConstraintResidual)&&r.normalizedMomentError<=1e-10&&r.normalizedConstraintResidual<=tolerance&&r.minimumCentralInertiaEigenvalue>=-128*std::numeric_limits<double>::epsilon()*scale[2];
        r.diagnostic=r.passed?"BINARY_MOMENT_VALIDATED: numerical cuboid-list check only":"BINARY_VALIDATION_FAILED: moment, central inertia or rounded paper constraint gate";
    }catch(const std::exception& e){r.passed=false;r.diagnostic=std::string("BINARY_VALIDATION_FAILED: ")+e.what();}
    return r;
}
}
