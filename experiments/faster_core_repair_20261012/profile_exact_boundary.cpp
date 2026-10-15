#include "verification/ExactBoundaryAudit.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>

// Diagnostic executable only: the ordinary exact auditor is linked unchanged.
int main(int argc,char** argv) {
    try {
        if(argc!=3)throw std::invalid_argument("mesh tape and repetition count required");
        std::ifstream input(argv[1]);
        int vertices=0,faces=0;
        if(!(input>>vertices>>faces) || vertices<=0 || faces<=0)
            throw std::invalid_argument("invalid mesh tape size");
        std::vector<std::array<std::string,3>> points(vertices);
        Eigen::MatrixXi triangles(faces,3);
        for(auto& point:points)for(auto& value:point)
            if(!(input>>value))throw std::invalid_argument("incomplete exact point tape");
        for(int i=0;i<faces;++i)for(int j=0;j<3;++j)
            if(!(input>>triangles(i,j)))throw std::invalid_argument("incomplete face tape");
        const int repetitions=std::stoi(argv[2]);
        if(repetitions<1 || repetitions>20)throw std::invalid_argument("invalid repetitions");
        for(int iteration=0;iteration<repetitions;++iteration) {
            const auto start=std::chrono::steady_clock::now();
            const auto audit=verification::auditExactTriangleBoundaryRational(points,triangles,
                verification::ExactBoundaryOrientation::AnyNonzero);
            const auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now()-start).count();
            const auto& facts=audit.facts();
            if(!facts.boundaryValid)throw std::runtime_error("saved actual boundary rejected");
            std::cout<<elapsed<<' '<<facts.inputDigest<<' '<<facts.backend<<' '
                     <<facts.exactPairTests<<' '<<facts.componentSignedVolumesExact.size();
            for(const auto& volume:facts.componentSignedVolumesExact)std::cout<<' '<<volume;
            std::cout<<'\n'<<std::flush;
        }
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
