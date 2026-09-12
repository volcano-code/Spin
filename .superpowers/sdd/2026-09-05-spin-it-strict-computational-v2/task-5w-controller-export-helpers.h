#pragma once
// Controller-only raw export. No admission flags or integration in this helper.
#include "spin2014/SourceCagePaperFunctional.h"
#include <iomanip>
#include <iostream>

template<class T> void arrayJson(const T& a) {
    std::cout << '[';
    for (std::size_t i=0;i<std::size_t(a.size());++i) {
        if(i)std::cout<<',';
        std::cout<<a[i];
    }
    std::cout << ']';
}
template<class T> void matrixJson(const T& a) {
    std::cout << '[';
    for(Eigen::Index i=0;i<a.rows();++i) {
        if(i)std::cout<<',';
        const auto row=a.row(i).eval();arrayJson(row);
    }
    std::cout << ']';
}
void dumpSource(const spin2014::SourceCagePaperFunctional& f,const Eigen::VectorXd& x) {
    const auto e=f.evaluate(x);
    const auto& b=f.source().boundary();
    const auto& bind=b.binding();
    const auto& grid=bind.bbw().fem().grid();
    const auto& cage=grid.cage();
    std::cout<<"{\"state\":";arrayJson(x);
    std::cout<<",\"betas\":"<<f.betaCount()<<",\"handles\":"<<f.handleCount()
        <<",\"angle\":"<<f.hasAngle()<<",\"handle_offset\":"<<f.handleOffset()
        <<",\"density\":"<<f.source().density()
        <<",\"partition_hash\":"<<std::quoted(bind.domain().partitionHash())
        <<",\"bbw_hash\":"<<std::quoted(bind.bbw().outputSha256())
        <<",\"source_hash\":"<<std::quoted(bind.source().sourceSha256())
        <<",\"config_hash\":"<<std::quoted(f.configSha256())
        <<",\"upstream_hash\":"<<std::quoted(f.upstreamSha256())
        <<",\"evaluation_hash\":"<<std::quoted(e.outputSha256)
        <<",\"moments\":";arrayJson(e.world.finalMoments.values());
    std::cout<<",\"beta_jacobian\":";matrixJson(e.world.betaJacobian);
    std::cout<<",\"rest\":";matrixJson(b.points());
    std::cout<<",\"points\":";matrixJson(e.world.deformedPoints);
    std::cout<<",\"outer_faces\":";matrixJson(bind.source().outerFaces());
    std::cout<<",\"centers\":[";
    for(std::size_t h=0;h<f.handleCount();++h) {
        if(h)std::cout<<',';
        const auto p=cage.boundaryVertices().row(cage.handleVertexIndices()[h]).eval();arrayJson(p);
    }
    std::cout<<"],\"grid_nodes\":";matrixJson(grid.nodePositions());
    std::cout<<",\"weights\":";matrixJson(bind.bbw().finalWeights());
    std::cout<<",\"grid_cells\":[";
    for(std::size_t i=0;i<grid.cellNodes().size();++i){if(i)std::cout<<',';arrayJson(grid.cellNodes()[i]);}
    std::cout<<"],\"leaves\":[";
    const auto& leaves=bind.domain().leaves();
    for(std::size_t i=0;i<leaves.size();++i) {
        if(i)std::cout<<',';
        std::cout<<"{\"level\":"<<leaves[i].level<<",\"key\":";arrayJson(leaves[i].coordinates);
        std::cout<<",\"moments\":";arrayJson(leaves[i].moments.values());std::cout<<'}';
    }
    std::cout<<"],\"primitives\":[";
    for(std::size_t p=0;p<bind.primitives().size();++p) {
        if(p)std::cout<<',';
        const auto& prim=bind.primitives()[p];
        std::cout<<"{\"owner\":"<<prim.betaIndex<<",\"level\":"<<prim.sourceKey.level<<",\"key\":";arrayJson(prim.sourceKey.coordinates);
        std::cout<<",\"faces\":[";
        const auto& refs=b.primitiveBoundaries()[p];
        for(std::size_t k=0;k<refs.size();++k) {
            if(k)std::cout<<',';
            auto t=b.triangles()[refs[k].triangle];if(refs[k].orientation<0)std::swap(t[1],t[2]);arrayJson(t);
        }
        std::cout<<"]}";
    }
    std::cout<<"]}";
}
