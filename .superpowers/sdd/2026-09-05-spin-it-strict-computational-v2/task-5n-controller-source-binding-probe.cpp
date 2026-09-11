#include "spin2014/SourceBinaryMoments.h"
#include "source_surface_test_utils.h"
#include <iostream>
#include <iomanip>
#include <map>
int main() {
    using namespace spin2014;
    auto outer=surface_test::cube(.3,.7), design=surface_test::cube(.34,.66);
    SourceSurfaceOctreeConfig config; config.maximumDepth=5;
    auto source=std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,config);
    SourceAdaptiveDomain domain(source,2,5);
    const auto primitives=sourceBinaryCuboids(domain,Eigen::VectorXd::Ones(domain.leaves().size()),100000);
    std::map<SourcePartitionKey,std::size_t> owners;
    for(std::size_t i=0;i<domain.leaves().size();++i)owners.emplace(domain.leaves()[i],i);
    Eigen::MatrixXd sums=Eigen::MatrixXd::Zero(10,domain.leaves().size());
    for(const auto& primitive:primitives.cuboids) {
        SourcePartitionKey ancestor=primitive.cell; int count=0; std::size_t owner=0;
        for(;;) {
            auto found=owners.find(ancestor); if(found!=owners.end()){++count;owner=found->second;}
            if(!ancestor.level)break;
            --ancestor.level;for(auto& k:ancestor.coordinates)k>>=1;
        }
        if(count!=1)throw std::runtime_error("nonunique owner");
        const auto l=primitive.lower,u=primitive.upper;
        const auto volume=(u-l).prod(); const Eigen::Vector3d mean=(u+l)/2;
        const Eigen::Vector3d second=(u.cwiseProduct(u)+u.cwiseProduct(l)+l.cwiseProduct(l))/3;
        Eigen::Matrix<double,10,1> moments;
        moments<<volume,volume*mean.x(),volume*mean.y(),volume*mean.z(),
            volume*second.x(),volume*second.y(),volume*second.z(),
            volume*mean.y()*mean.z(),volume*mean.x()*mean.z(),volume*mean.x()*mean.y();
        sums.col(owner)+=moments;
    }
    double maxError=0,boundingVolume=0;
    for(std::size_t i=0;i<domain.leaves().size();++i){
        maxError=std::max(maxError,(sums.col(i)-domain.leaves()[i].moments.values()).cwiseAbs().maxCoeff());
        boundingVolume+=std::pow(std::ldexp(1.,-domain.leaves()[i].level),3);
    }
    const double expectedVolume=std::pow(10./32,3);
    std::cout<<std::setprecision(17)<<"SYNTHETIC_REGRESSION source-binding preflight only\n"
        <<"source="<<source->sourceSha256()<<"\npartition="<<domain.partitionHash()
        <<"\nleaves="<<domain.leaves().size()<<" primitives="<<primitives.cuboids.size()
        <<"\nactual_volume="<<sums.row(0).sum()<<" expected_volume="<<expectedVolume
        <<" bounding_cell_volume="<<boundingVolume<<"\nmax_leaf_moment_error="<<maxError<<'\n';
    if(maxError>1e-12 || std::abs(sums.row(0).sum()-expectedVolume)>1e-12 ||
       primitives.cuboids.size()<=domain.leaves().size() || boundingVolume<=expectedVolume)
        throw std::runtime_error("source-binding audit mismatch");
}
