// Actual per-column partial-run rejection, never mocked candidate status.
#include "spin2014/SourceCageBbw.h"
#include <iostream>
#include <stdexcept>
void require(bool ok,const char* why) {if(!ok) throw std::runtime_error(why);}
int main() {
    try {
        Eigen::MatrixXd v(8,3);
        for(int b=0;b<8;++b) for(int a=0;a<3;++a) v(b,a)=((b>>a)&1)?.7:.3;
        Eigen::MatrixXi f(12,3);
        f << 0,2,1,1,2,3,4,5,6,5,7,6,0,4,2,2,4,6,1,3,5,3,7,5,0,1,4,1,5,4,2,6,3,3,6,7;
        spin2014::SourceSurfaceCageConfig cc;cc.cageLevel=2;
        spin2014::SourceCageGridConfig gc;gc.fineLevel=3;
        const spin2014::SourceCageGrid grid(spin2014::SourceSurfaceCage(v,f,cc),gc);
        for(auto rule:{spin2014::SourceCageFemMassRule::RESTRICTED_GRID_VORONOI,
                       spin2014::SourceCageFemMassRule::BARYCENTRIC_TET}) {
            spin2014::SourceCageFemConfig fc;fc.massRule=rule;
            const spin2014::SourceCageFem fem(grid,fc);
            spin2014::SourceCageBbwConfig config;
            config.protocol=spin2014::SourceCageBbwProtocol::INDEPENDENT_NORMALIZED;
            config.iterationLimit=30;
            const spin2014::SourceCageBbw run(fem,config);
            require(!run.success(),"partial solve granted success");
            require(run.terminal()==spin2014::SourceCageBbwTerminal::ITERATION_LIMIT,"not actual iteration exhaustion");
            require(run.finalWeights().size()==0,"partial solve produced normalized field");
            require(!run.rawAssessment().available && !run.finalAssessment().available,"partial field got complete assessment");
            require(run.attempts().size()>1 && run.attempts().size()<grid.handleNodeIndices().size(),"not a partially completed solve");
            const auto& first=run.attempts().front();const auto& last=run.attempts().back();
            require(first.runInvoked && first.candidateAvailable && first.assessment.optimalityPassed,"first column never validated");
            require(last.runInvoked && last.qpIterations>=30,"last actual backend did not exhaust iteration budget");
            std::cout << "PASS partial run mass=" << int(rule) << " calls=" << run.attempts().size()
                      << " firstIterations=" << first.qpIterations << " lastIterations=" << last.qpIterations
                      << " rawAssessmentAvailable=" << run.rawAssessment().available
                      << " finalValues=" << run.finalWeights().size() << " terminal=" << run.terminalReason()
                      << " scientificSha256=" << run.outputSha256() << '\n';
        }
        return 0;
    } catch(const std::exception& e) {std::cerr << "FAIL " << e.what() << '\n';return 1;}
}
