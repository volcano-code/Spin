// Headless adapter for the pinned Cherchi et al. Boolean component.
// SYNTHETIC_REGRESSION / BOOLEAN_COMPONENT_ONLY, not an author solver.
#include "booleans.h"
#include "io_functions.h"

#include <tbb/global_control.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

std::string canonicalCoordinate(const double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value == 0.0 ? 0.0 : value);
    return buffer;
}

std::string roundTripCoordinate(const double value) {
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

void replaceWithNearestExactCoordinates(
    const std::vector<const genericPoint*>& outputPoints,
    const double multiplier,
    std::vector<double>& outputCoords) {
    if (!std::isfinite(multiplier) || multiplier == 0.0 ||
        outputCoords.size() != outputPoints.size() * 3U) {
        throw std::runtime_error("invalid exact Cherchi coordinate export inputs");
    }
    const bigrational exactMultiplier(multiplier);
    for (std::size_t vertex = 0; vertex < outputPoints.size(); ++vertex) {
        const genericPoint* point = outputPoints[vertex];
        if (point == nullptr || (!point->isExplicit3D() && !point->isLPI()))
            throw std::runtime_error("exact binary64 export requires explicit Cherchi points");
        bigrational x, y, z;
        if (!point->getExactXYZCoordinates(x, y, z))
            throw std::runtime_error("Cherchi exact output point evaluation failed");
        const std::array<double, 3> exactRounded = {
            nearestBinary64(x / exactMultiplier),
            nearestBinary64(y / exactMultiplier),
            nearestBinary64(z / exactMultiplier),
        };
        // This function is exclusive to --dump-exact-rational. The selected
        // arrangement walk already authenticated point/triangle identities.
        // Approximate implicit constructions can lose many ulps through
        // cancellation; their proximity must not decide whether the exact
        // export is actually exact-rounded. Ordinary/provenance-only output
        // keeps the original approximate coordinates. Topology/containment of
        // this derived binary64 representation still requires a separate audit.
        for (std::size_t axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(exactRounded[axis]))
                throw std::runtime_error("exact Cherchi coordinate is not representable as finite binary64");
            outputCoords[3U * vertex + axis] = exactRounded[axis];
        }
    }
}

std::array<uint, 3> orientedMinimumFirst(std::array<uint, 3> face) {
    const auto minimum = std::min_element(face.begin(), face.end());
    std::rotate(face.begin(), minimum, face.end());
    return face;
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

void writeExactCoordinateSidecar(
    const char* filename,
    const std::vector<const genericPoint*>& outputPoints,
    const std::vector<uint>& canonicalOrder,
    const double multiplier,
    const std::vector<double>& outputCoords) {
    if (!std::isfinite(multiplier) || multiplier == 0.0) {
        throw std::runtime_error("invalid Cherchi output scale multiplier");
    }
    const bigrational exactMultiplier(multiplier);
    if (outputCoords.size() % 3U != 0U ||
        canonicalOrder.size() != outputCoords.size() / 3U)
        throw std::runtime_error("output coordinate array has invalid size");
    if (outputPoints.size() != outputCoords.size() / 3U)
        throw std::runtime_error("Cherchi output point/coordinate count mismatch");
    const int fd = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) throw std::runtime_error(std::string("cannot create exact sidecar: ") + std::strerror(errno));
    std::FILE* output = fdopen(fd, "w");
    if (output == nullptr) {
        close(fd);
        throw std::runtime_error("cannot stream exact sidecar");
    }
    std::fprintf(output, "vertex_index\tstatus\texact_x\texact_y\texact_z\t"
                         "memory_x\tmemory_y\tmemory_z\tobj_x\tobj_y\tobj_z\n");
    const std::size_t vertexCount = canonicalOrder.size();
    for (std::size_t vertex = 0; vertex < vertexCount; ++vertex) {
        const uint sourceVertex = canonicalOrder[vertex];
        if (sourceVertex >= outputPoints.size())
            throw std::runtime_error("canonical exact-sidecar vertex index is invalid");
        const std::array<std::string, 3> key = {
            canonicalCoordinate(outputCoords[3U * sourceVertex]),
            canonicalCoordinate(outputCoords[3U * sourceVertex + 1U]),
            canonicalCoordinate(outputCoords[3U * sourceVertex + 2U]),
        };
        const std::array<std::string, 3> memory = {
            roundTripCoordinate(outputCoords[3U * sourceVertex]),
            roundTripCoordinate(outputCoords[3U * sourceVertex + 1U]),
            roundTripCoordinate(outputCoords[3U * sourceVertex + 2U]),
        };
        const genericPoint* point = outputPoints[sourceVertex];
        if (point == nullptr) throw std::runtime_error("null Cherchi output point");
        if (!point->isExplicit3D() && !point->isLPI()) {
            const char* kind = point->isTPI() ? "SKIPPED_TPI" :
                               point->isLNC() ? "SKIPPED_LNC" :
                               "SKIPPED_UNSUPPORTED_POINT_TYPE";
            std::fprintf(output, "%zu\t%s\t-\t-\t-\t%s\t%s\t%s\t%s\t%s\t%s\n",
                         vertex, kind, memory[0].c_str(), memory[1].c_str(), memory[2].c_str(),
                         key[0].c_str(), key[1].c_str(), key[2].c_str());
            continue;
        }
        bigrational exactX, exactY, exactZ;
        if (!point->getExactXYZCoordinates(exactX, exactY, exactZ))
            throw std::runtime_error("Cherchi exact output point evaluation failed");
        const std::array<std::string, 3> exact = {
            exactRationalString(exactX / exactMultiplier),
            exactRationalString(exactY / exactMultiplier),
            exactRationalString(exactZ / exactMultiplier),
        };
        const std::array<double, 3> exactRounded = {
            nearestBinary64(exactX / exactMultiplier),
            nearestBinary64(exactY / exactMultiplier),
            nearestBinary64(exactZ / exactMultiplier),
        };
        const bool coordinateMatchesExact =
            exactRounded[0] == outputCoords[3U * sourceVertex] &&
            exactRounded[1] == outputCoords[3U * sourceVertex + 1U] &&
            exactRounded[2] == outputCoords[3U * sourceVertex + 2U];
        std::fprintf(output, "%zu\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                     vertex, coordinateMatchesExact
                         ? "MATCHED" : "APPROXIMATION_DISAGREES",
                     exact[0].c_str(), exact[1].c_str(), exact[2].c_str(),
                     memory[0].c_str(), memory[1].c_str(), memory[2].c_str(),
                     key[0].c_str(), key[1].c_str(), key[2].c_str());
    }
    const bool writeFailed = std::ferror(output) != 0;
    if (std::fclose(output) != 0 || writeFailed)
        throw std::runtime_error("failed to write exact sidecar");
}

std::string pointTypeToken(const genericPoint& point) {
    if (point.isExplicit3D()) return "EXPLICIT3D";
    if (point.isLPI()) return "LPI";
    if (point.isTPI()) return "TPI";
    if (point.isLNC()) return "LNC";
    return "UNSUPPORTED";
}

void writeConstructionSidecar(
    const char* filename,
    const std::vector<const genericPoint*>& outputPoints,
    const std::vector<uint>& canonicalOrder,
    const double multiplier) {
    if (!std::isfinite(multiplier) || multiplier == 0.0 ||
        outputPoints.size() != canonicalOrder.size()) {
        throw std::runtime_error("invalid construction provenance inputs");
    }
    const bigrational exactMultiplier(multiplier);
    const int fd = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) throw std::runtime_error(std::string("cannot create construction sidecar: ") + std::strerror(errno));
    std::FILE* output = fdopen(fd, "w");
    if (output == nullptr) {
        close(fd);
        throw std::runtime_error("cannot stream construction sidecar");
    }
    std::fprintf(output, "vertex_index\tpoint_type\toperand_count\tparameter_t\t"
                         "operand0_x\toperand0_y\toperand0_z\t"
                         "operand1_x\toperand1_y\toperand1_z\t"
                         "operand2_x\toperand2_y\toperand2_z\t"
                         "operand3_x\toperand3_y\toperand3_z\t"
                         "operand4_x\toperand4_y\toperand4_z\t"
                         "operand5_x\toperand5_y\toperand5_z\t"
                         "operand6_x\toperand6_y\toperand6_z\t"
                         "operand7_x\toperand7_y\toperand7_z\t"
                         "operand8_x\toperand8_y\toperand8_z\n");
    for (std::size_t vertex = 0; vertex < canonicalOrder.size(); ++vertex) {
        const uint pointIndex = canonicalOrder[vertex];
        if (pointIndex >= outputPoints.size() || outputPoints[pointIndex] == nullptr)
            throw std::runtime_error("construction provenance vertex index is invalid");
        const genericPoint& point = *outputPoints[pointIndex];
        std::vector<const explicitPoint3D*> operands;
        double parameter = 0.0;
        if (point.isExplicit3D()) {
            operands.push_back(&point.toExplicit3D());
        } else if (point.isLPI()) {
            const auto& p = point.toLPI();
            operands = {&p.P(), &p.Q(), &p.R(), &p.S(), &p.T()};
        } else if (point.isTPI()) {
            const auto& p = point.toTPI();
            operands = {&p.V1(), &p.V2(), &p.V3(), &p.W1(), &p.W2(), &p.W3(),
                        &p.U1(), &p.U2(), &p.U3()};
        } else if (point.isLNC()) {
            const auto& p = point.toLNC();
            operands = {&p.P(), &p.Q()};
            parameter = p.T();
        } else {
            throw std::runtime_error("unsupported Cherchi point construction type");
        }
        std::fprintf(output, "%zu\t%s\t%zu\t%.17g", vertex,
                     pointTypeToken(point).c_str(), operands.size(), parameter);
        for (const explicitPoint3D* operand : operands) {
            bigrational x, y, z;
            if (!operand->getExactXYZCoordinates(x, y, z))
                throw std::runtime_error("failed to evaluate Cherchi construction operand");
            const std::array<std::string, 3> exact = {
                exactRationalString(x / exactMultiplier),
                exactRationalString(y / exactMultiplier),
                exactRationalString(z / exactMultiplier),
            };
            std::fprintf(output, "\t%s\t%s\t%s", exact[0].c_str(),
                         exact[1].c_str(), exact[2].c_str());
        }
        for (std::size_t operand = operands.size(); operand < 9; ++operand)
            std::fprintf(output, "\t-\t-\t-");
        std::fprintf(output, "\n");
    }
    const bool writeFailed = std::ferror(output) != 0;
    if (std::fclose(output) != 0 || writeFailed)
        throw std::runtime_error("failed to write construction sidecar");
}

void writeInputObjectLabelSidecar(
    const char* filename,
    const std::vector<std::bitset<NBIT>>& faceLabels) {
    const int fd = open(filename, O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) throw std::runtime_error(std::string("cannot create face-label sidecar: ") + std::strerror(errno));
    std::FILE* output = fdopen(fd, "w");
    if (output == nullptr) {
        close(fd);
        throw std::runtime_error("cannot stream face-label sidecar");
    }
    std::fprintf(output, "face_index\tinput_object_mask\n");
    for (std::size_t face = 0; face < faceLabels.size(); ++face) {
        std::fprintf(output, "%zu\t0x%llx\n", face,
                     static_cast<unsigned long long>(faceLabels[face].to_ullong()));
    }
    const bool writeFailed = std::ferror(output) != 0;
    if (std::fclose(output) != 0 || writeFailed)
        throw std::runtime_error("failed to write face-label sidecar");
}

} // namespace

int main(int argc, char** argv) {
    const bool dumpExact = argc == 6 && std::string(argv[5]) == "--dump-exact-rational";
    const bool dumpActualProvenance = argc == 6 &&
        std::string(argv[5]) == "--dump-actual-provenance";
    const bool dumpProvenance = dumpExact || dumpActualProvenance;
    if (argc != 5 && !dumpProvenance) {
        std::cerr << "usage: cherchi_probe intersection|subtraction inputA.obj inputB.obj output.obj [--dump-exact-rational|--dump-actual-provenance]\n";
        return 2;
    }
    const std::string operation(argv[1]);
    const BoolOp op = operation == "intersection" ? INTERSECTION :
                      operation == "subtraction" ? SUBTRACTION : NONE;
    if (op == NONE) {
        std::cerr << "unsupported Boolean operation: " << operation << '\n';
        return 2;
    }

    const tbb::global_control single_thread(tbb::global_control::max_allowed_parallelism, 1);
    std::vector<double> in_coords, out_coords;
    std::vector<uint> in_tris, in_labels, out_tris;
    std::vector<std::bitset<NBIT>> out_labels;
    loadMultipleFiles({argv[2], argv[3]}, in_coords, in_tris, in_labels);
    point_arena exactArena;
    std::vector<genericPoint*> exactArrangementPoints;
    std::vector<const genericPoint*> exactOutputPoints;
    if (dumpProvenance) {
        runBooleanPipelineRetainingArrangementPoints(
            in_coords, in_tris, in_labels, op, exactArena,
            exactArrangementPoints, exactOutputPoints, out_coords, out_tris, out_labels);
        if (exactArrangementPoints.size() < 5)
            throw std::runtime_error("arrangement point list lacks scale sentinel");
        const double multiplier = exactArrangementPoints.back()->toExplicit3D().X();
        alignExactPointsWithApproximateCoordinates(exactOutputPoints, multiplier,
                                                   out_coords);
        if (dumpExact)
            replaceWithNearestExactCoordinates(exactOutputPoints, multiplier, out_coords);
    } else {
        booleanPipeline(in_coords, in_tris, in_labels, op, out_coords, out_tris, out_labels);
    }
    if (out_coords.empty() || out_tris.empty()) {
        std::cerr << "Boolean result is empty\n";
        return 3;
    }
    const int fd = open(argv[4], O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) {
        std::cerr << "cannot create output exclusively " << argv[4] << ": " << std::strerror(errno) << '\n';
        return 4;
    }
    std::FILE* output = fdopen(fd, "w");
    if (!output) {
        std::cerr << "cannot stream output " << argv[4] << ": " << std::strerror(errno) << '\n';
        close(fd);
        return 4;
    }
    const std::size_t vertexCount = out_coords.size() / 3;
    std::vector<std::array<std::string, 3>> coordinateKeys(vertexCount);
    for (std::size_t vertex = 0; vertex < vertexCount; ++vertex) {
        coordinateKeys[vertex] = {
            canonicalCoordinate(out_coords[3 * vertex]),
            canonicalCoordinate(out_coords[3 * vertex + 1]),
            canonicalCoordinate(out_coords[3 * vertex + 2]),
        };
    }
    std::vector<uint> order(vertexCount);
    std::iota(order.begin(), order.end(), 0U);
    std::sort(order.begin(), order.end(), [&](const uint lhs, const uint rhs) {
        return coordinateKeys[lhs] < coordinateKeys[rhs];
    });
    std::vector<uint> oldToNew(vertexCount);
    for (std::size_t newIndex = 0; newIndex < vertexCount; ++newIndex) {
        oldToNew[order[newIndex]] = static_cast<uint>(newIndex);
        const auto& key = coordinateKeys[order[newIndex]];
        std::fprintf(output, "v %s %s %s\n", key[0].c_str(), key[1].c_str(), key[2].c_str());
    }
    struct CanonicalFace {
        std::array<uint, 3> vertices;
        std::bitset<NBIT> inputObjects;
    };
    if (dumpProvenance && out_labels.size() != out_tris.size() / 3U)
        throw std::runtime_error("Cherchi output face labels do not align with output triangles");
    std::vector<CanonicalFace> faces;
    faces.reserve(out_tris.size() / 3);
    for (std::size_t triangle = 0; triangle < out_tris.size() / 3U; ++triangle) {
        const std::size_t index = 3U * triangle;
        faces.push_back({orientedMinimumFirst({
                             oldToNew[out_tris[index]], oldToNew[out_tris[index + 1]],
                             oldToNew[out_tris[index + 2]],
                         }),
                         dumpProvenance ? out_labels[triangle] : std::bitset<NBIT>{}});
    }
    std::sort(faces.begin(), faces.end(), [](const CanonicalFace& lhs, const CanonicalFace& rhs) {
        if (lhs.vertices != rhs.vertices) return lhs.vertices < rhs.vertices;
        return lhs.inputObjects.to_ullong() < rhs.inputObjects.to_ullong();
    });
    for (const auto& face : faces)
        std::fprintf(output, "f %u %u %u\n", face.vertices[0] + 1,
                     face.vertices[1] + 1, face.vertices[2] + 1);
    const bool write_failed = std::ferror(output) != 0;
    if (std::fclose(output) != 0 || write_failed) {
        std::cerr << "failed to write output\n";
        return 5;
    }
    if (dumpProvenance) {
        const std::string sidecar = std::string(argv[4]) + ".exact.tsv";
        if (exactArrangementPoints.size() < 5)
            throw std::runtime_error("arrangement point list lacks sentinels");
        const double multiplier = exactArrangementPoints.back()->toExplicit3D().X();
        writeExactCoordinateSidecar(sidecar.c_str(), exactOutputPoints, order,
                                    multiplier, out_coords);
        const std::string constructionSidecar = std::string(argv[4]) + ".construction.tsv";
        writeConstructionSidecar(constructionSidecar.c_str(), exactOutputPoints,
                                 order, multiplier);
        const std::string faceSidecar = std::string(argv[4]) + ".faces.tsv";
        std::vector<std::bitset<NBIT>> sortedFaceLabels;
        sortedFaceLabels.reserve(faces.size());
        for (const auto& face : faces) sortedFaceLabels.push_back(face.inputObjects);
        writeInputObjectLabelSidecar(faceSidecar.c_str(), sortedFaceLabels);
    }
    std::cout << "op=" << operation << " input_vertices=" << in_coords.size() / 3
              << " input_triangles=" << in_tris.size() / 3
              << " output_vertices=" << out_coords.size() / 3
              << " output_triangles=" << out_tris.size() / 3 << '\n';
    return 0;
}
