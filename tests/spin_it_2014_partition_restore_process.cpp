// Test-only file exchange across fresh processes. NOT a production checkpoint
// format: no saved pass history, provenance manifest, durable write or scheduler.
#include "spin2014/SourceAdaptiveSpinSolver.h"
#include "spin2014/SlqpIterationStore.h"
#include "verification/Sha256.h"
#include "source_surface_test_utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <iterator>

using namespace spin2014;
using namespace surface_test;
std::shared_ptr<const SourceSurfaceOctree> source() {
    const auto outer=box({-.5+1./8192,-.5,0},{.5+1./8192,.5,1});
    const auto design=box({.125,-.125,.625},{.375,.125,.875});
    SourceSurfaceOctreeConfig c;c.maximumDepth=5;c.rootOrigin={-1,-1,0};c.rootSide=2;
    return std::make_shared<const SourceSurfaceOctree>(outer.v,outer.f,design.v,design.f,c);
}
SlqpLbfgsResult solve(const SourceAdaptiveDomain& d,const Eigen::VectorXd& beta,double phi) {
    AdaptiveFillFunctional f(d,SpinIt2014MomentProblem::fromObjective(Spin2014Objective::StaticLowering,1,1),1,1);
    SlqpLbfgsSettings c;c.algorithm=SlqpAlgorithm::AuditedKktRestorationV1;
    c.directionPolicy=SlqpDirectionPolicy::FiniteBoxPrimalV2;
    c.maximumIterations=300;c.equalityTolerance=1e-10;c.projectedGradientTolerance=1e-8;
    return solveSlqpLbfgsPaper(f,beta,phi,c);
}
void vector(std::ostream& o,const Eigen::VectorXd& v) {
    o<<v.size()<<' ';for(Eigen::Index i=0;i<v.size();++i)o<<v[i]<<' ';o<<'\n';
}
std::string digest(const std::string& bytes) {
    verification::Sha256 h;
    h.update(reinterpret_cast<const unsigned char*>(bytes.data()),bytes.size());
    return h.finish();
}
int main(int argc,char** argv){try {
    require(argc==3,"mode and test-state path required");
    const std::string mode=argv[1];require(mode=="--produce" || mode=="--consume","bad mode");
    auto s=source();SourceAdaptiveDomain d(s,2,5,10000);
    Eigen::VectorXd beta;double phi=0;
    if(mode=="--produce") {
        auto first=solve(d,Eigen::VectorXd::Constant(d.leaves().size(),.5),0);
        require(first.status==SlqpStatus::Converged && first.trace.size()>1,"first stage not a nontrivial solve");
        beta=first.beta;phi=first.finalPhi;
        const auto split=d.splitFractionalLeaves(beta,.05);d.mergeBinarySiblings(beta,.05);
        require(!split.changes.empty(),"actual stage transition not exercised");
        std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
        out<<"PARTITION_TEST_ONLY_V2\n"<<d.partitionHash()<<'\n'<<d.leaves().size()<<' '<<phi<<'\n';
        for(std::size_t i=0;i<d.leaves().size();++i){const auto& k=d.leaves()[i];
            out<<k.level<<' '<<k.coordinates[0]<<' '<<k.coordinates[1]<<' '<<k.coordinates[2]<<' '<<beta[i]<<'\n';}
        auto file=std::fopen(argv[2],"wx");require(file,"exclusive test-state creation failed");
        const auto body=out.str();
        const auto bytes=body+"PARTITION_TEST_ONLY_END_V2 "+std::to_string(body.size())+" "+digest(body)+"\n";
        const bool wrote=std::fwrite(bytes.data(),1,bytes.size(),file)==bytes.size();
        const int closed=std::fclose(file);require(wrote && closed==0,"test-state write failed");
    } else {
        std::ifstream file(argv[2],std::ios::binary);require(file.good(),"bad test-state file");
        const std::string bytes(std::istreambuf_iterator<char>{file},std::istreambuf_iterator<char>{});
        const auto footer=bytes.rfind("PARTITION_TEST_ONLY_END_V2 ");
        require(footer!=std::string::npos && footer>0 && bytes[footer-1]=='\n',"missing test-state footer");
        const auto body=bytes.substr(0,footer);
        std::istringstream trailer(bytes.substr(footer));trailer.imbue(std::locale::classic());
        std::string marker,expectedHash;std::size_t expectedBytes=0;
        trailer>>marker>>expectedBytes>>expectedHash;
        require(!trailer.fail() && marker=="PARTITION_TEST_ONLY_END_V2" &&
            expectedBytes==body.size() && expectedHash==digest(body),"corrupt test-state framing");
        trailer>>std::ws;require(trailer.eof(),"trailing test-state footer data");
        std::istringstream in(body);in.imbue(std::locale::classic());std::string magic,hash;std::size_t n=0;
        in>>magic>>hash>>n>>phi;require(in.good() && magic=="PARTITION_TEST_ONLY_V2" && n>0 && n<=10000,"bad test-state header");
        std::vector<SourcePartitionKey> keys(n);beta.resize(n);
        for(std::size_t i=0;i<n;++i)in>>keys[i].level>>keys[i].coordinates[0]>>keys[i].coordinates[1]>>keys[i].coordinates[2]>>beta[i];
        require(!in.fail(),"truncated test-state");in>>std::ws;require(in.eof(),"trailing test-state data");
        require(phi==0,"static fixture phi changed");
        d=SourceAdaptiveDomain::restorePartition(s,2,5,10000,keys,hash);
        (void)d.removedMoments(beta); // existing finite [0,1] validation
    }
    // Compare the actual next SLQP, not just a reconstructed domain's getters.
    std::cout<<std::hexfloat<<"SYNTHETIC_SINGLE_BOUNDARY_NOT_FULL_RESUME\n"<<d.partitionHash()<<'\n'<<phi<<'\n';
    vector(std::cout,beta);
    auto next=solve(d,beta,phi);require(next.status==SlqpStatus::Converged,"next stage failed");
    require(next.trace.size()>1,"post-restoration solve no longer exercises actual iterations");
    std::cout<<next.traceHash<<'\n'<<next.finalPhi<<'\n';vector(std::cout,next.beta);
    vector(std::cout,next.evaluation.finalMoments.values());
    for(const auto& t:next.trace){const auto bytes=encodeSlqpIteration(t);std::cout<<bytes.size()<<':';std::cout.write(bytes.data(),bytes.size());}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
