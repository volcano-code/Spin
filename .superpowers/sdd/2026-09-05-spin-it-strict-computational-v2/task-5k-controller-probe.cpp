// Actual native BBW export for independent verification; no release authority.
#include "spin2014/SourceCageBbw.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

template<class Matrix> void matrix(const Matrix& x) {
    std::cout << '[';
    for(Eigen::Index i=0;i<x.rows();++i) {
        if(i) std::cout << ',';
        std::cout << '[';
        for(Eigen::Index j=0;j<x.cols();++j) {if(j) std::cout << ',';std::cout << x(i,j);}
        std::cout << ']';
    }
    std::cout << ']';
}
template<class Range> void values(const Range& x) {
    std::cout << '['; bool first=true;
    for(const auto& item:x) {if(!first) std::cout << ',';first=false;std::cout << item;}
    std::cout << ']';
}
void assessment(const spin2014::SourceCageBbwAssessment& x) {
    std::cout << "{\"available\":" << x.available << ",\"primal_passed\":" << x.primalPassed
              << ",\"optimality_passed\":" << x.optimalityPassed
              << ",\"energy\":" << x.energy << ",\"energy_scale\":" << x.energyScale
              << ",\"boundary\":" << x.boundaryResidual << ",\"bounds\":" << x.boundsResidual
              << ",\"partition\":" << x.partitionResidual << ",\"stationarity\":" << x.stationarityResidual
              << ",\"dual\":" << x.dualResidual << ",\"complementarity\":" << x.complementarityResidual
              << ",\"projected_gradient\":" << x.projectedGradientResidual << '}';
}
int main(int argc,char** argv) {
    try {
        if(argc!=4) throw std::invalid_argument("usage: probe cells_per_axis mass_rule protocol");
        const int cells=std::stoi(argv[1]);
        if(cells!=2 && cells!=4) throw std::invalid_argument("only locked analytic fixtures");
        const std::string rule=argv[2],protocol=argv[3];
        spin2014::SourceCageFemConfig femConfig;
        if(rule=="BARYCENTRIC_TET") femConfig.massRule=spin2014::SourceCageFemMassRule::BARYCENTRIC_TET;
        else if(rule!="RESTRICTED_GRID_VORONOI") throw std::invalid_argument("mass rule");
        spin2014::SourceCageBbwConfig solveConfig;
        if(protocol=="INDEPENDENT_NORMALIZED") solveConfig.protocol=spin2014::SourceCageBbwProtocol::INDEPENDENT_NORMALIZED;
        else if(protocol!="COUPLED_PARTITION") throw std::invalid_argument("protocol");
        Eigen::MatrixXd vertices(8,3);
        for(int b=0;b<8;++b) for(int a=0;a<3;++a) vertices(b,a)=((b>>a)&1)?(cells==2?.4:.7):.3;
        Eigen::MatrixXi faces(12,3);
        faces << 0,2,1,1,2,3,4,5,6,5,7,6,0,4,2,2,4,6,1,3,5,3,7,5,0,1,4,1,5,4,2,6,3,3,6,7;
        spin2014::SourceSurfaceCageConfig cageConfig;cageConfig.cageLevel=2;
        spin2014::SourceCageGridConfig gridConfig;gridConfig.fineLevel=3;
        const spin2014::SourceCageGrid grid(spin2014::SourceSurfaceCage(vertices,faces,cageConfig),gridConfig);
        const spin2014::SourceCageFem fem(grid,femConfig);
        const spin2014::SourceCageBbw result(fem,solveConfig);
        if(!result.success()) throw std::runtime_error(result.terminalReason());
        std::cout << std::boolalpha << std::setprecision(17)
                  << "{\"scope\":\"NATIVE_SYNTHETIC_BBW_NOT_PAPER_CASE\",\"cells_per_axis\":" << cells
                  << ",\"rule\":\"" << rule << "\",\"protocol\":\"" << protocol
                  << "\",\"config_sha256\":\"" << result.configSha256()
                  << "\",\"output_sha256\":\"" << result.outputSha256()
                  << "\",\"fem_output_sha256\":\"" << fem.outputSha256()
                  << "\",\"success\":" << result.success() << ",\"terminal\":" << int(result.terminal())
                  << ",\"energy_scale\":" << result.energyScale()
                  << ",\"normalization_change\":" << result.maximumNormalizationChange()
                  << ",\"nodes\":";matrix(grid.nodePositions());
        std::cout << ",\"handle_nodes\":";values(grid.handleNodeIndices());
        std::cout << ",\"raw_weights\":";matrix(result.rawWeights());
        std::cout << ",\"final_weights\":";matrix(result.finalWeights());
        std::cout << ",\"raw_assessment\":";assessment(result.rawAssessment());
        std::cout << ",\"final_assessment\":";assessment(result.finalAssessment());
        std::cout << ",\"attempts\":[";bool first=true;
        for(const auto& attempt:result.attempts()) {
            if(!first) std::cout << ',';first=false;
            std::cout << "{\"handle\":" << attempt.handle << ",\"return_status\":" << attempt.returnStatus
                      << ",\"pass_model_status\":" << attempt.passModelStatus
                      << ",\"run_invoked\":" << attempt.runInvoked
                      << ",\"diagnostics_truncated\":" << attempt.diagnosticsTruncated
                      << ",\"model_status\":" << attempt.modelStatus << ",\"qp_iterations\":" << attempt.qpIterations
                      << ",\"candidate_available\":" << attempt.candidateAvailable << ",\"assessment\":";
            assessment(attempt.assessment);std::cout << '}';
        }
        std::cout << "]}\n";
        return std::cout.good()?0:1;
    } catch(const std::exception& e) {std::cerr << e.what() << '\n';return 1;}
}
