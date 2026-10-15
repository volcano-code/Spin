#pragma once
#include "voxel/MaterialOccupancy.h"
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace dual_geometry_export {
namespace fs=std::filesystem;
inline std::ofstream file(const fs::path& p) {
    std::ofstream o;o.exceptions(std::ios::failbit|std::ios::badbit);o.open(p);
    o.imbue(std::locale::classic());o<<std::setprecision(17);return o;
}
template<class T> void vector(std::ostream& o,const T& v) {
    o<<'[';for(Eigen::Index i=0;i<v.size();++i){if(i)o<<',';o<<v[i];}o<<']';
}
inline void binaryMask(std::ostream& o,const std::vector<unsigned char>& values) {
    o<<'[';for(std::size_t i=0;i<values.size();++i){if(i)o<<',';o<<unsigned(values[i]);}o<<']';
}
inline void matrix(std::ostream& o,const Eigen::Matrix3d& values) {
    o<<'[';
    for(int row=0;row<3;++row) {
        if(row)o<<',';o<<'[';
        for(int column=0;column<3;++column){if(column)o<<',';o<<values(row,column);}
        o<<']';
    }
    o<<']';
}
inline void exportCase(const fs::path& directory,const spin2014::DualMaterialSpinRequest& request,
                       const spin2014::DualMaterialSpinResult& result) {
    if(!fs::create_directory(directory))throw std::invalid_argument("material output directory exists");
    const auto& domain=*request.hollowing.domain;
    if(result.materialLabels.size()!=domain.cellCount())throw std::invalid_argument("missing material labels");
    std::array<std::size_t,2> counts{};
    std::size_t fixedShell=0,removedCarveable=0;
    for(auto label:result.materialLabels){if(label>2)throw std::invalid_argument("invalid material label");if(label)++counts[label-1];}
    for(std::size_t index=0;index<domain.cellCount();++index) {
        if(domain.material(index)&&!domain.carveable(index))++fixedShell;
        if(domain.material(index)&&domain.carveable(index)&&result.materialLabels[index]==0U)++removedCarveable;
    }
    auto o=file(directory/"manifest.json");
    o<<"{\"schema\":\"dual-material-geometry-v1\",\"classification\":\"SYNTHETIC_REGRESSION\","
     <<"\"strict_result\":\""<<(counts[1]?"METHOD_GEOMETRY_ONLY":"NO_RETAINED_DENSE_INSERT")
     <<"\",\"solver_status\":\""<<spin2014::adaptiveSpinStatusToken(result.hollowing.status)
     <<"\",\"solver_output_sha256\":\""<<result.outputHash<<"\",\"input_domain_sha256\":\""<<domain.artifactId()
     <<"\",\"mode\":\""<<spin::spinModeToken(request.hollowing.mode)
     <<"\",\"legacy_binary_tolerance\":"<<request.hollowing.binaryConstraintTolerance
     <<",\"binary_constraint_tolerance\":"<<request.hollowing.binaryConstraintTolerance
     <<",\"binary_constraint_inf\":"<<result.hollowing.binaryConstraintInfinityNorm
     <<",\"binary_constraints\":";vector(o,result.hollowing.binaryConstraints);
    o<<",\"target_frame\":{\"origin\":";vector(o,request.hollowing.targetFrame.origin);
    o<<",\"axes\":";matrix(o,request.hollowing.targetFrame.axes);o<<'}';
    o<<",\"partition\":{\"void_cells\":"<<(domain.cellCount()-counts[0]-counts[1])
     <<",\"boundary_cells\":"<<counts[0]<<",\"interior_cells\":"<<counts[1]
     <<",\"fixed_shell_cells\":"<<fixedShell<<",\"removed_carveable_cells\":"<<removedCarveable<<'}'
     <<",\"total_moments\":";vector(o,result.hollowing.finalMoments.values());
    o<<",\"materials\":[";
    for(unsigned label=1;label<=2;++label) {
        if(label>1)o<<',';
        const double density=label==1?request.boundaryDensity:request.interiorDensity;
        const auto moments=label==1?result.boundaryContribution.values():result.interiorContribution.values();
        const std::string filename=label==1?"low_density_shell.obj":"high_density_insert.obj";
        o<<"{\"label\":"<<label<<",\"density\":"<<density<<",\"cell_count\":"<<counts[label-1]<<",\"moments\":";vector(o,moments);
        o<<",\"mesh\":";
        if(!counts[label-1]) {o<<"null}";continue;}
        std::vector<unsigned char> mask(domain.cellCount(),0);
        for(std::size_t i=0;i<mask.size();++i)mask[i]=result.materialLabels[i]==label;
        const auto material=voxel::VoxelDomain::create(domain.grid(),mask);
        const auto surface=voxel::extractMaterialBoundaries(material);
        auto mesh=file(directory/filename);
        for(Eigen::Index i=0;i<surface.outerVertices.rows();++i)mesh<<"v "<<surface.outerVertices.row(i)<<'\n';
        for(Eigen::Index i=0;i<surface.innerVertices.rows();++i)mesh<<"v "<<surface.innerVertices.row(i)<<'\n';
        for(Eigen::Index i=0;i<surface.outerFaces.rows();++i)mesh<<"f "<<surface.outerFaces(i,0)+1<<' '<<surface.outerFaces(i,1)+1<<' '<<surface.outerFaces(i,2)+1<<'\n';
        const auto offset=surface.outerVertices.rows()+1;
        for(Eigen::Index i=0;i<surface.innerFaces.rows();++i)mesh<<"f "<<surface.innerFaces(i,0)+offset<<' '<<surface.innerFaces(i,1)+offset<<' '<<surface.innerFaces(i,2)+offset<<'\n';
        o<<'"'<<filename<<"\"}";
    }
    o<<"],\"labels\":[";
    for(std::size_t i=0;i<result.materialLabels.size();++i){if(i)o<<',';o<<unsigned(result.materialLabels[i]);}
    o<<"],\"input_material\":";binaryMask(o,domain.materialMask());
    o<<",\"input_carveable\":";binaryMask(o,domain.carveableMask());
    o<<",\"grid_dimensions\":";vector(o,domain.grid().dimensions);
    o<<",\"grid_minimum\":";vector(o,domain.grid().minimum);
    o<<",\"grid_maximum\":";vector(o,domain.grid().maximum);
    o<<",\"sealed\":false,\"original_paper_case\":false}\n";
}
}
