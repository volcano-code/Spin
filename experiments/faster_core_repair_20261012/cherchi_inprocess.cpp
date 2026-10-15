// Reuse the frozen probe's arrangement/implicit-point implementation exactly.
// Its CLI entry is renamed and never executed by this module.
#define main cherchi_preserved_probe_cli
#include "../cherchi_reference_probe/probe.cpp"
#undef main
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py=pybind11;

PYBIND11_MODULE(_cherchi_core_boolean,module) {
    module.def("boolean",[](const std::vector<std::array<double,3>>& a,
        const std::vector<std::array<uint,3>>& af,const std::vector<std::array<double,3>>& b,
        const std::vector<std::array<uint,3>>& bf,const std::string& operation) {
        if(a.size()<4 || b.size()<4 || af.size()<4 || bf.size()<4)
            throw py::value_error("closed triangle operands required");
        const BoolOp op=operation=="intersection"?INTERSECTION:operation=="difference"?SUBTRACTION:NONE;
        if(op==NONE)throw py::value_error("unsupported Boolean operation");
        std::vector<double> coords,outCoords;
        std::vector<uint> triangles,labels,outTris;
        const auto append=[&](const auto& vertices,const auto& faces,uint label,uint offset) {
            for(const auto& v:vertices)for(double x:v) {
                if(!std::isfinite(x))throw py::value_error("nonfinite input coordinate");
                coords.push_back(x);
            }
            for(const auto& face:faces) {
                for(uint index:face) {
                    if(index>=vertices.size())throw py::value_error("input face index out of range");
                    triangles.push_back(index+offset);
                }
                if(face[0]==face[1] || face[1]==face[2] || face[0]==face[2])
                    throw py::value_error("repeated input face index");
                labels.push_back(label);
            }
        };
        append(a,af,0,0);append(b,bf,1,static_cast<uint>(a.size()));
        point_arena arena;
        std::vector<genericPoint*> arrangement;
        std::vector<const genericPoint*> points;
        std::vector<std::bitset<NBIT>> outputLabels;
        const tbb::global_control singleThread(tbb::global_control::max_allowed_parallelism,1);
        {
            py::gil_scoped_release release;
            runBooleanPipelineRetainingArrangementPoints(coords,triangles,labels,op,arena,
                arrangement,points,outCoords,outTris,outputLabels);
        }
        if(arrangement.size()<5)throw std::runtime_error("missing scale sentinel");
        const double multiplier=arrangement.back()->toExplicit3D().X();
        alignExactPointsWithApproximateCoordinates(points,multiplier,outCoords);
        const bigrational scale(multiplier);
        std::vector<std::array<double,3>> vertices;
        std::vector<std::array<std::string,3>> exact;
        py::list records;
        for(std::size_t index=0;index<points.size();++index) {
            const genericPoint& point=*points[index];
            bigrational x,y,z;
            if(!point.getExactXYZCoordinates(x,y,z))throw std::runtime_error("exact point unavailable");
            vertices.push_back({nearestBinary64(x/scale),nearestBinary64(y/scale),nearestBinary64(z/scale)});
            exact.push_back({exactRationalString(x/scale),exactRationalString(y/scale),exactRationalString(z/scale)});
            std::vector<const explicitPoint3D*> operands;double parameter=0.;
            if(point.isExplicit3D())operands={&point.toExplicit3D()};
            else if(point.isLPI()) {
                const auto& p=point.toLPI();operands={&p.P(),&p.Q(),&p.R(),&p.S(),&p.T()};
            } else if(point.isTPI()) {
                const auto& p=point.toTPI();operands={&p.V1(),&p.V2(),&p.V3(),&p.W1(),&p.W2(),&p.W3(),&p.U1(),&p.U2(),&p.U3()};
            } else if(point.isLNC()) {
                const auto& p=point.toLNC();operands={&p.P(),&p.Q()};parameter=p.T();
            } else throw std::runtime_error("unsupported point construction");
            std::vector<std::array<std::string,3>> operandCoordinates;
            for(const auto* operand:operands) {
                if(!operand->getExactXYZCoordinates(x,y,z))throw std::runtime_error("operand exact point unavailable");
                operandCoordinates.push_back({exactRationalString(x/scale),exactRationalString(y/scale),exactRationalString(z/scale)});
            }
            py::dict record;record["vertex_index"]=index;record["point_type"]=pointTypeToken(point);
            record["parameter_t"]=parameter;record["operand_coordinates_exact"]=operandCoordinates;
            records.append(record);
        }
        std::vector<std::array<uint,3>> faces;
        for(std::size_t i=0;i<outTris.size();i+=3)faces.push_back({outTris[i],outTris[i+1],outTris[i+2]});
        py::dict result;result["vertices"]=vertices;result["faces"]=faces;result["exact_vertices"]=exact;
        result["construction_records"]=records;result["normalization_multiplier"]=multiplier;
        result["backend"]="cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0-inprocess-v1";
        return result;
    });
}
