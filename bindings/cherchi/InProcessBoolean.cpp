// In-process adapter for pinned Cherchi 7bd6c269.
// Geometry helpers extracted unchanged from the validated reference adapter.
#include "booleans.h"
#include "io_functions.h"
#include <tbb/global_control.h>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
namespace {
std::string canonicalCoordinate(const double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value == 0.0 ? 0.0 : value);
    return buffer;
}

double nearestBinary64(const bigrational& exact) {
    if (exact.sgn() == 0) return 0.0;
    const bignatural numerator = exact.get_num();
    const bignatural denominator = exact.get_den();
    int exponent = static_cast<int>(numerator.getNumSignificantBits()) -
                   static_cast<int>(denominator.getNumSignificantBits());
    if (exponent >= 0) {
        bignatural scaledDenominator = denominator;
        scaledDenominator <<= static_cast<uint32_t>(exponent);
        if (numerator < scaledDenominator) --exponent;
    } else {
        bignatural scaledNumerator = numerator;
        scaledNumerator <<= static_cast<uint32_t>(-exponent);
        if (scaledNumerator < denominator) --exponent;
    }

    const auto roundedScaledQuotient = [](bignatural scaledNumerator,
                                          bignatural scaledDenominator) {
        bignatural remainder;
        bignatural quotient = scaledNumerator.divide_by(scaledDenominator, remainder);
        bignatural twiceRemainder = remainder;
        if (!twiceRemainder.empty()) twiceRemainder <<= 1U;
        const bool exactlyHalfway =
            !(twiceRemainder > scaledDenominator) &&
            !(twiceRemainder < scaledDenominator);
        if (twiceRemainder > scaledDenominator ||
            (exactlyHalfway && quotient.getBit(0))) {
            quotient = quotient + bignatural(1U);
        }
        uint64_t value = 0;
        if (!quotient.toUint64(value))
            throw std::runtime_error("rounded Cherchi coordinate significand overflow");
        return value;
    };

    uint64_t magnitudeBits = 0;
    if (exponent > 1023) {
        magnitudeBits = UINT64_C(0x7ff0000000000000);
    } else if (exponent >= -1022) {
        bignatural scaledNumerator = numerator;
        bignatural scaledDenominator = denominator;
        const int shift = 52 - exponent;
        if (shift >= 0) scaledNumerator <<= static_cast<uint32_t>(shift);
        else scaledDenominator <<= static_cast<uint32_t>(-shift);
        uint64_t significand = roundedScaledQuotient(scaledNumerator, scaledDenominator);
        if (significand == (UINT64_C(1) << 53U)) {
            significand >>= 1U;
            ++exponent;
        }
        if (exponent > 1023) magnitudeBits = UINT64_C(0x7ff0000000000000);
        else magnitudeBits = (static_cast<uint64_t>(exponent + 1023) << 52U) |
                             (significand - (UINT64_C(1) << 52U));
    } else {
        bignatural scaledNumerator = numerator;
        scaledNumerator <<= 1074U;
        const uint64_t significand = roundedScaledQuotient(scaledNumerator, denominator);
        magnitudeBits = significand >= (UINT64_C(1) << 52U)
            ? (UINT64_C(1) << 52U) : significand;
    }
    if (exact.sgn() < 0) magnitudeBits |= UINT64_C(1) << 63U;
    double result = 0.0;
    std::memcpy(&result, &magnitudeBits, sizeof(result));
    return result;
}

void alignExactPointsWithApproximateCoordinates(
    std::vector<const genericPoint*>& outputPoints,
    const double multiplier,
    const std::vector<double>& outputCoords) {
    if (!std::isfinite(multiplier) || multiplier == 0.0 ||
        outputCoords.size() != outputPoints.size() * 3U) {
        throw std::runtime_error("invalid Cherchi point/coordinate alignment inputs");
    }
    // Both this adapter and computeFinalExplicitResult bind vertex indices by
    // first encounter in the selected triangle walk. Keep that binding:
    // approximate coordinates are not injective for implicit exact points.
    // Re-matching collision buckets by coordinates can swap distinct points
    // without changing their displayed coordinates and fold the exact mesh.
    for (std::size_t vertex = 0; vertex < outputPoints.size(); ++vertex) {
        const genericPoint* point = outputPoints[vertex];
        if (point == nullptr) throw std::runtime_error("null Cherchi output point");
        double x = 0.0, y = 0.0, z = 0.0;
        if (!point->getApproxXYZCoordinates(x, y, z))
            throw std::runtime_error("Cherchi approximate output point evaluation failed");
        const std::array<double, 3> approximate = {x / multiplier, y / multiplier, z / multiplier};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double output = outputCoords[3U * vertex + axis];
            if (!std::isfinite(approximate[axis]) || !std::isfinite(output) ||
                canonicalCoordinate(approximate[axis]) != canonicalCoordinate(output))
                throw std::runtime_error("Cherchi arrangement/output vertex index binding mismatch");
        }
    }
}

std::string exactRationalString(const bigrational& value) {
    if (value.sgn() == 0) return "0/1";
    return std::string(value.sgn() < 0 ? "-" : "") +
           "0b" + value.get_num().get_str() + "/0b" +
           value.get_den().get_str();
}

void runBooleanPipelineRetainingArrangementPoints(
    const std::vector<double>& inCoords,
    const std::vector<uint>& inTris,
    const std::vector<uint>& inLabels,
    const BoolOp op,
    point_arena& arena,
    std::vector<genericPoint*>& arrVerts,
    std::vector<const genericPoint*>& outputPoints,
    std::vector<double>& outCoords,
    std::vector<uint>& outTris,
    std::vector<std::bitset<NBIT>>& outLabels) {
    std::vector<uint> arrInTris, arrOutTris;
    std::vector<std::bitset<NBIT>> arrInLabels;
    std::vector<DuplTriInfo> duplicateTriangles;
    Labels labels;
    std::vector<phmap::flat_hash_set<uint>> patches;
    cinolib::Octree octree;

    customArrangementPipeline(inCoords, inTris, inLabels, arrInTris,
                              arrInLabels, arena, arrVerts, arrOutTris,
                              labels, octree, duplicateTriangles, true);
    FastTrimesh tm(arrVerts, arrOutTris, true);
    computeAllPatches(tm, labels, patches, true);
    addDuplicateTrisInfoInStructures(duplicateTriangles, arrInTris,
                                     arrInLabels, octree);
    const cinolib::vec3d maxCoords(
        octree.root->bbox.max.x() + 0.5,
        octree.root->bbox.max.y() + 0.5,
        octree.root->bbox.max.z() + 0.5);
    computeInsideOut(tm, patches, octree, arrVerts, arrInTris, arrInLabels,
                     maxCoords, labels);
    uint selectedTriangleCount = 0;
    if (op == INTERSECTION) selectedTriangleCount = boolIntersection(tm, labels);
    else if (op == UNION) selectedTriangleCount = boolUnion(tm, labels);
    else if (op == SUBTRACTION) selectedTriangleCount = boolSubtraction(tm, labels);
    else if (op == XOR) selectedTriangleCount = boolXOR(tm, labels);
    else throw std::runtime_error("unsupported Cherchi diagnostic Boolean operation");
    computeFinalExplicitResult(tm, labels, selectedTriangleCount, outCoords,
                               outTris, outLabels, true);

    std::vector<int> vertexIndex(tm.numVerts(), -1);
    outputPoints.assign(outCoords.size() / 3U, nullptr);
    std::size_t nextVertex = 0;
    std::size_t selectedTriangle = 0;
    if (outTris.size() != 3U * selectedTriangleCount)
        throw std::runtime_error("Cherchi output triangle count differs from selected arrangement");
    for (uint triangle = 0; triangle < tm.numTris(); ++triangle) {
        if (tm.triInfo(triangle) == 0) continue;
        const uint* vertices = tm.tri(triangle);
        for (uint corner = 0; corner < 3; ++corner) {
            const uint oldVertex = vertices[corner];
            if (vertexIndex[oldVertex] < 0) {
                vertexIndex[oldVertex] = static_cast<int>(nextVertex);
                outputPoints[nextVertex++] = tm.vert(oldVertex);
            }
            if (outTris[3U * selectedTriangle + corner] !=
                static_cast<uint>(vertexIndex[oldVertex]))
                throw std::runtime_error("Cherchi arrangement/output triangle index binding mismatch");
        }
        ++selectedTriangle;
    }
    if (nextVertex != outputPoints.size() || selectedTriangle != selectedTriangleCount)
        throw std::runtime_error("Cherchi exact point order differs from output mesh");
}

std::string pointTypeToken(const genericPoint& point) {
    if (point.isExplicit3D()) return "EXPLICIT3D";
    if (point.isLPI()) return "LPI";
    if (point.isTPI()) return "TPI";
    if (point.isLNC()) return "LNC";
    return "UNSUPPORTED";
}


}
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

