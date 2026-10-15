#include "spin/QuadricMomentEvaluator.h"
#include "spin/QuadricDiscreteMomentDifferential.h"
#include "spin/PolyhedralMomentDifferential.h"
#include "spin/QuadricSourceMotion.h"

#include "manufacturing/ManufacturingGeometry.h"
#include "manufacturing/SlicerValidation.h"
#include "physics/MassProperties.h"
#include "verification/ExactMaterialRegion.h"
#include "verification/ExactBoundaryAudit.h"
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
#include "verification/ExactLibiglBoolean.h"
#include "verification/ExactRationalCoordinates.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <igl/copyleft/cgal/assign_scalar.h>
#endif

#include <igl/MeshBooleanType.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/convex_hull.h>
#include <igl/winding_number.h>
#include <igl/bfs_orient.h>
#include <CGAL/Gmpq.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <chrono>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <tuple>
#include <vector>

#ifndef MAKE_IT_STAND_LIBIGL_SOURCE_ID
#define MAKE_IT_STAND_LIBIGL_SOURCE_ID "unversioned"
#endif


namespace
{

#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
using RationalPoint = std::array<CGAL::Gmpq,3>;
using RationalCoordinates = std::vector<std::array<std::string,3>>;

std::string rationalText(const CGAL::Gmpq& value)
{
    std::ostringstream out; out.imbue(std::locale::classic()); out << value;
    return out.str();
}

RationalCoordinates rationalTargetVertices(const spin::TriangleMesh& mesh,
    const spin::SpinTargetFrame& frame)
{
    RationalCoordinates out(mesh.vertices.rows());
    for (int v=0;v<mesh.vertices.rows();++v) for(int a=0;a<3;++a)
    {
        CGAL::Gmpq value(0);
        for(int b=0;b<3;++b)
            value += CGAL::Gmpq(frame.axes(b,a)) *
                (CGAL::Gmpq(mesh.vertices(v,b))-CGAL::Gmpq(frame.origin(b)));
        out[v][a]=rationalText(value);
    }
    return out;
}

void assignRationalVertices(const RationalCoordinates& coordinates, Eigen::MatrixXd* vertices)
{
    vertices->resize(coordinates.size(),3);
    for(std::size_t v=0;v<coordinates.size();++v)for(int a=0;a<3;++a)
        igl::copyleft::cgal::assign_scalar(
            CGAL::Epeck::FT(verification::parseExactRationalCoordinate(coordinates[v][a])),(*vertices)(v,a));
}

std::vector<std::string> rationalSurfaceMoments(const RationalCoordinates& coordinates,
    const Eigen::MatrixXi& faces, double density)
{
    std::vector<RationalPoint> points(coordinates.size());
    for(std::size_t v=0;v<points.size();++v)for(int a=0;a<3;++a)
        points[v][a]=verification::parseExactRationalCoordinate(coordinates[v][a]);
    std::array<CGAL::Gmpq,10> moments{};
    constexpr std::array<std::array<int,2>,6> products{{{{0,0}},{{1,1}},{{2,2}},{{1,2}},{{0,2}},{{0,1}}}};
    for(int f=0;f<faces.rows();++f)
    {
        const auto& a=points.at(faces(f,0)); const auto& b=points.at(faces(f,1)); const auto& c=points.at(faces(f,2));
        const auto volume=(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6;
        moments[0]+=volume;
        RationalPoint sum;
        for(int i=0;i<3;++i){sum[i]=a[i]+b[i]+c[i];moments[1+i]+=volume*sum[i]/4;}
        for(int i=0;i<6;++i)
        {
            const int u=products[i][0],v=products[i][1];
            moments[4+i]+=volume*(sum[u]*sum[v]+a[u]*a[v]+b[u]*b[v]+c[u]*c[v])/20;
        }
    }
    std::vector<std::string> out; out.reserve(10);
    for(auto& value:moments)out.push_back(rationalText(value*CGAL::Gmpq(density)));
    return out;
}

MassMoments assignedRationalMoments(const std::vector<std::string>& exact)
{
    MassMoments out;
    std::array<double,10> m;
    for(int i=0;i<10;++i)igl::copyleft::cgal::assign_scalar(
        CGAL::Epeck::FT(verification::parseExactRationalCoordinate(exact.at(i))),m[i]);
    out.mass=m[0];out.firstMoment=Eigen::Vector3d(m[1],m[2],m[3]);
    out.secondMoment<<m[4],m[9],m[8],m[9],m[5],m[7],m[8],m[7],m[6];
    return out;
}
#endif

// A sufficient, deliberately incomplete certificate for represented doubles.
// Exact rational interval arithmetic bounds the interior, not just vertices.
bool certifyNegativeOnDesignBox(
    const Eigen::MatrixXd& vertices,
    const spin::QuadricCoefficients10& quadric,
    std::string* certificate
)
{
    using Q = CGAL::Gmpq;
    if (vertices.rows() == 0 || vertices.cols() != 3 || !vertices.allFinite())
        return false;
    std::array<Q, 3> lo, hi;
    for (int i = 0; i < 3; ++i)
    {
        lo[i] = Q(vertices.col(i).minCoeff());
        hi[i] = Q(vertices.col(i).maxCoeff());
    }
    const auto& p = quadric.values();
    Q upper(p(0));
    const auto addUpper = [&](const Q& coefficient, const Q& lower, const Q& higher) {
        upper += coefficient * (coefficient >= Q(0) ? higher : lower);
    };
    for (int i = 0; i < 3; ++i)
    {
        addUpper(Q(p(1+i)), lo[i], hi[i]);
        const Q a = lo[i]*lo[i], b = hi[i]*hi[i];
        const Q lower = lo[i] <= Q(0) && hi[i] >= Q(0) ? Q(0) : std::min(a,b);
        addUpper(Q(p(4+i)), lower, std::max(a,b));
    }
    const std::array<std::array<int,3>,3> cross = {{{1,2,7},{0,2,8},{0,1,9}}};
    for (const auto& term : cross)
    {
        const int i=term[0], j=term[1];
        const std::array<Q,4> products = {{lo[i]*lo[j],lo[i]*hi[j],hi[i]*lo[j],hi[i]*hi[j]}};
        addUpper(Q(2)*Q(p(term[2])), *std::min_element(products.begin(),products.end()),
            *std::max_element(products.begin(),products.end()));
    }
    if (!(upper < Q(0))) return false;
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    identity << "ExactGmpqAabbNegative-v1;upper=" << upper << ';';
    for (int i=0;i<3;++i) identity << lo[i] << ',' << hi[i] << ';';
    *certificate = identity.str();
    return true;
}

// A truncated unbounded quadric is a valid proxy only when its artificial
// closure is provably outside the entire design domain. The domain is inside
// the convex hull of its boundary vertices; each cap triangle must be strictly
// separated from that hull. Exact
// rational projection intervals make this a conservative proof; failure to
// find a separating direction remains unknown and is rejected.
void requireTruncationCapsOutsideDesignDomain(
    const Eigen::MatrixXd& designVertices,
    const spin::QuadricTriangulation& triangulation,
    const std::vector<std::array<std::string,3>>* exactDesign=nullptr
)
{
    if (triangulation.capFaceCount == 0) return;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(!triangulation.exactVertices.empty()) {
        if(!exactDesign || exactDesign->empty() ||
           triangulation.faceIsCap.size()!=triangulation.closedBoundary.faces.rows())
            throw std::runtime_error("rational truncation cap requires authority domain and provenance");
        RationalPoint minimum,maximum;
        for(int axis=0;axis<3;++axis) {
            minimum[axis]=maximum[axis]=verification::parseExactRationalCoordinate(exactDesign->at(0)[axis]);
            for(const auto& point:*exactDesign) {
                const auto x=verification::parseExactRationalCoordinate(point[axis]);
                minimum[axis]=std::min(minimum[axis],x);maximum[axis]=std::max(maximum[axis],x);
            }
        }
        for(int f=0;f<triangulation.closedBoundary.faces.rows();++f)if(triangulation.faceIsCap(f)) {
            bool separated=false;
            for(int axis=0;axis<3 && !separated;++axis) {
                bool above=true,below=true;
                for(int k=0;k<3;++k) {
                    const auto x=verification::parseExactRationalCoordinate(
                        triangulation.exactVertices.at(triangulation.closedBoundary.faces(f,k))[axis]);
                    above&=x>maximum[axis];below&=x<minimum[axis];
                }
                separated=above || below;
            }
            if(!separated)throw std::runtime_error("rational truncation cap is not certified outside authority design domain");
        }
        return;
    }
#endif
    if (designVertices.cols() != 3 || designVertices.rows() == 0 ||
        !designVertices.allFinite() ||
        triangulation.faceIsCap.size() != triangulation.closedBoundary.faces.rows())
        throw std::runtime_error("truncation cap is not certified outside design domain: incomplete cap provenance");

    using Q = CGAL::Gmpq;
    const auto& vertices = triangulation.closedBoundary.vertices;
    const auto& faces = triangulation.closedBoundary.faces;

    // The historical five-direction probe is deliberately retained as a
    // cheap first pass.  If it cannot certify a cap, use the complete
    // separating-axis set for the convex hull of the design boundary and the
    // cap triangle.  This is still a conservative certificate for the actual
    // design domain because the domain is contained in that hull, while it
    // avoids rejecting valid oblique caps merely because their separator is
    // not axis-aligned or parallel to the cap plane.
    Eigen::MatrixXd hullVertices;
    Eigen::MatrixXi hullFaces;
    bool hullReady = false;
    const auto exactHullSeparated = [&](const std::array<Eigen::Vector3d,3>& capPoints) {
        if (!hullReady)
        {
            try
            {
                igl::copyleft::cgal::convex_hull(designVertices, hullVertices, hullFaces);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error("truncation cap is not certified outside design domain: convex hull construction failed");
            }
            if (hullVertices.cols() != 3 || hullVertices.rows() < 4 ||
                hullFaces.cols() != 3 || hullFaces.rows() == 0 ||
                !hullVertices.allFinite() || (hullFaces.array() < 0).any() ||
                (hullFaces.array() >= hullVertices.rows()).any())
                throw std::runtime_error("truncation cap is not certified outside design domain: degenerate design convex hull");
            hullReady = true;
        }

        std::vector<Eigen::Vector3d> axes;
        axes.reserve(static_cast<std::size_t>(hullFaces.rows() + 4 * hullFaces.rows()));
        const auto addAxis = [&](const Eigen::Vector3d& axis) {
            if (axis.allFinite() && axis.squaredNorm() > 1e-24)
                axes.push_back(axis);
        };
        std::set<std::pair<int,int>> hullEdges;
        for (Eigen::Index f = 0; f < hullFaces.rows(); ++f)
        {
            const Eigen::Vector3d a = hullVertices.row(hullFaces(f,0)).transpose();
            const Eigen::Vector3d b = hullVertices.row(hullFaces(f,1)).transpose();
            const Eigen::Vector3d c = hullVertices.row(hullFaces(f,2)).transpose();
            addAxis((b-a).cross(c-a));
            for (int k = 0; k < 3; ++k)
            {
                int first = hullFaces(f,k), second = hullFaces(f,(k+1)%3);
                if (first > second) std::swap(first,second);
                hullEdges.emplace(first,second);
            }
        }
        addAxis((capPoints[1]-capPoints[0]).cross(capPoints[2]-capPoints[0]));
        for (int k = 0; k < 3; ++k)
            for (const auto& edge : hullEdges)
            {
                const Eigen::Vector3d hullEdge =
                    hullVertices.row(edge.second).transpose() -
                    hullVertices.row(edge.first).transpose();
                const Eigen::Vector3d capEdge =
                    capPoints[(k+1)%3] - capPoints[k];
                addAxis(hullEdge.cross(capEdge));
            }

        const auto projectionBounds = [&](const Eigen::MatrixXd& points,
                                          const Eigen::Vector3d& axis,
                                          Q* minimum, Q* maximum) {
            bool first = true;
            for (Eigen::Index row = 0; row < points.rows(); ++row)
            {
                const Q value = Q(axis.x()) * Q(points(row,0)) +
                    Q(axis.y()) * Q(points(row,1)) + Q(axis.z()) * Q(points(row,2));
                if (first) { *minimum = *maximum = value; first = false; }
                else { *minimum = std::min(*minimum,value); *maximum = std::max(*maximum,value); }
            }
        };
        for (const auto& axis : axes)
        {
            Q hullMinimum, hullMaximum, capMinimum, capMaximum;
            projectionBounds(hullVertices,axis,&hullMinimum,&hullMaximum);
            Eigen::MatrixXd capMatrix(3,3);
            for (int i = 0; i < 3; ++i) capMatrix.row(i) = capPoints[i].transpose();
            projectionBounds(capMatrix,axis,&capMinimum,&capMaximum);
            if (capMaximum < hullMinimum || hullMaximum < capMinimum)
                return true;
        }
        return false;
    };
    Eigen::Vector3d designCenter = Eigen::Vector3d::Zero();
    for (Eigen::Index row = 0; row < designVertices.rows(); ++row)
        designCenter += designVertices.row(row).transpose();
    designCenter /= static_cast<double>(designVertices.rows());
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        if (triangulation.faceIsCap(face) == 0) continue;
        const int a = faces(face, 0), b = faces(face, 1), c = faces(face, 2);
        if (a < 0 || b < 0 || c < 0 || a >= vertices.rows() ||
            b >= vertices.rows() || c >= vertices.rows())
            throw std::runtime_error("truncation cap is not certified outside design domain: invalid cap face");
        const Eigen::Vector3d pa = vertices.row(a).transpose();
        const Eigen::Vector3d pb = vertices.row(b).transpose();
        const Eigen::Vector3d pc = vertices.row(c).transpose();
        const Eigen::Vector3d normal = (pb - pa).cross(pc - pa);
        if (!normal.allFinite() || normal.squaredNorm() == 0.0)
            throw std::runtime_error("truncation cap is not certified outside design domain");

        const std::array<Eigen::Vector3d, 5> directions = {{
            Eigen::Vector3d::UnitX(), Eigen::Vector3d::UnitY(), Eigen::Vector3d::UnitZ(),
            normal, (pa + pb + pc) / 3.0 - designCenter
        }};
        bool separated = false;
        for (const auto& direction : directions)
        {
            if (!direction.allFinite() || direction.squaredNorm() == 0.0) continue;
            Q designMinimum, designMaximum, capMinimum, capMaximum;
            bool first = true;
            const auto projection = [&](double x, double y, double z) {
                return Q(direction.x()) * Q(x) + Q(direction.y()) * Q(y) + Q(direction.z()) * Q(z);
            };
            for (Eigen::Index row = 0; row < designVertices.rows(); ++row)
            {
                const Q value = projection(designVertices(row, 0), designVertices(row, 1), designVertices(row, 2));
                if (first)
                {
                    designMinimum = designMaximum = value;
                    first = false;
                }
                else
                {
                    designMinimum = std::min(designMinimum, value);
                    designMaximum = std::max(designMaximum, value);
                }
            }
            const std::array<Eigen::Vector3d, 3> capPoints = {{pa, pb, pc}};
            first = true;
            for (const auto& point : capPoints)
            {
                const Q value = projection(point.x(), point.y(), point.z());
                if (first)
                {
                    capMinimum = capMaximum = value;
                    first = false;
                }
                else
                {
                    capMinimum = std::min(capMinimum, value);
                    capMaximum = std::max(capMaximum, value);
                }
            }
            if (capMinimum > designMaximum || capMaximum < designMinimum)
            {
                separated = true;
                break;
            }
        }
        if (!separated)
        {
            const std::array<Eigen::Vector3d, 3> capPoints = {{pa, pb, pc}};
            if (!exactHullSeparated(capPoints))
                throw std::runtime_error("truncation cap is not certified outside design domain");
        }
    }
}


spin::TriangleMesh toTarget(
    const spin::TriangleMesh& model,
    const spin::SpinTargetFrame& frame
)
{
    spin::TriangleMesh target = model;
    for (Eigen::Index row = 0; row < model.vertices.rows(); ++row)
    {
        target.vertices.row(row) =
            frame.axes.transpose() *
            (model.vertices.row(row).transpose() - frame.origin);
    }
    return target;
}


spin::TriangleMesh toModel(
    const spin::TriangleMesh& target,
    const spin::SpinTargetFrame& frame
)
{
    spin::TriangleMesh model = target;
    for (Eigen::Index row = 0; row < target.vertices.rows(); ++row)
    {
        model.vertices.row(row) =
            frame.axes * target.vertices.row(row).transpose() + frame.origin;
    }
    return model;
}


std::string meshHash(const spin::TriangleMesh& mesh, const char* version)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << version << ';' << mesh.vertices.rows() << ';'
           << mesh.faces.rows() << ';' << std::hexfloat;
    for (Eigen::Index row = 0; row < mesh.vertices.rows(); ++row)
    {
        stream << mesh.vertices(row, 0) << ';'
               << mesh.vertices(row, 1) << ';'
               << mesh.vertices(row, 2) << ';';
    }
    for (Eigen::Index row = 0; row < mesh.faces.rows(); ++row)
    {
        stream << mesh.faces(row, 0) << ';'
               << mesh.faces(row, 1) << ';'
               << mesh.faces(row, 2) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


void canonicalizeBooleanMesh(
    spin::TriangleMesh* mesh,
    Eigen::VectorXi* birthFaces,
    std::vector<std::array<std::string, 3>>* exactVertices = nullptr,
    bool rationalAuthority = false
)
{
    const bool hasExactVertices = exactVertices && !exactVertices->empty();
    if (rationalAuthority && !hasExactVertices)
        throw std::runtime_error("rational authority canonicalization requires actual coordinates");
    if (hasExactVertices &&
        exactVertices->size() != static_cast<std::size_t>(mesh->vertices.rows()))
    {
        throw std::runtime_error("invalid exact Boolean coordinate tape size");
    }
    for (Eigen::Index row = 0; row < mesh->vertices.rows(); ++row)
    {
        for (Eigen::Index column = 0; column < mesh->vertices.cols(); ++column)
        {
            if (mesh->vertices(row, column) == 0.0)
            {
                mesh->vertices(row, column) = 0.0;
            }
        }
    }
    std::vector<int> vertexOrder(
        static_cast<std::size_t>(mesh->vertices.rows())
    );
    for (int index = 0; index < mesh->vertices.rows(); ++index)
    {
        vertexOrder[static_cast<std::size_t>(index)] = index;
    }
    std::sort(
        vertexOrder.begin(),
        vertexOrder.end(),
        [&](int left, int right) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if (rationalAuthority)
            {
                for(int a=0;a<3;++a)
                {
                    const auto l=verification::parseExactRationalCoordinate((*exactVertices)[left][a]);
                    const auto r=verification::parseExactRationalCoordinate((*exactVertices)[right][a]);
                    if(l<r)return true; if(r<l)return false;
                }
                return left<right;
            }
#endif
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                if (mesh->vertices(left, coordinate) <
                    mesh->vertices(right, coordinate))
                {
                    return true;
                }
                if (mesh->vertices(right, coordinate) <
                    mesh->vertices(left, coordinate))
                {
                    return false;
                }
            }
            return left < right;
        }
    );
    std::vector<Eigen::Vector3d> uniqueVertices;
    uniqueVertices.reserve(vertexOrder.size());
    std::vector<int> uniqueInputs;
    uniqueInputs.reserve(vertexOrder.size());
    std::vector<int> remap(vertexOrder.size());
    for (const int input : vertexOrder)
    {
        const Eigen::Vector3d point = mesh->vertices.row(input);
        bool newVertex=uniqueVertices.empty() || point != uniqueVertices.back();
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        if(rationalAuthority && !uniqueVertices.empty())
        {
            newVertex=false;
            for(int a=0;a<3;++a)newVertex |=
                verification::parseExactRationalCoordinate((*exactVertices)[input][a]) !=
                verification::parseExactRationalCoordinate((*exactVertices)[uniqueInputs.back()][a]);
        }
#endif
        if (
            newVertex
        )
        {
            uniqueVertices.push_back(point);
            uniqueInputs.push_back(input);
        }
        else if (hasExactVertices)
        {
            // A floating-coordinate weld must not silently merge distinct
            // rational Boolean vertices. Only exact duplicates share a row.
            const auto& previous = (*exactVertices)[uniqueInputs.back()];
            const auto& current = (*exactVertices)[input];
            for (int axis = 0; axis < 3; ++axis)
            {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                if (verification::parseExactRationalCoordinate(previous[axis]) !=
                    verification::parseExactRationalCoordinate(current[axis]))
#else
                if (previous[axis] != current[axis])
#endif
                {
                    throw std::runtime_error(
                        "distinct exact Boolean vertices collide under floating assignment");
                }
            }
        }
        remap[static_cast<std::size_t>(input)] =
            static_cast<int>(uniqueVertices.size()) - 1;
    }

    struct CanonicalFace
    {
        std::array<int, 3> indices;
        int birth = 0;
    };
    std::vector<CanonicalFace> faces;
    faces.reserve(static_cast<std::size_t>(mesh->faces.rows()));
    for (Eigen::Index row = 0; row < mesh->faces.rows(); ++row)
    {
        std::array<int, 3> indices = {{
            remap[static_cast<std::size_t>(mesh->faces(row, 0))],
            remap[static_cast<std::size_t>(mesh->faces(row, 1))],
            remap[static_cast<std::size_t>(mesh->faces(row, 2))]
        }};
        const Eigen::Vector3d& a =
            uniqueVertices[static_cast<std::size_t>(indices[0])];
        const Eigen::Vector3d& b =
            uniqueVertices[static_cast<std::size_t>(indices[1])];
        const Eigen::Vector3d& c =
            uniqueVertices[static_cast<std::size_t>(indices[2])];
        if (
            indices[0] == indices[1] || indices[1] == indices[2] ||
            indices[2] == indices[0]
        )
        {
            if (hasExactVertices)
                throw std::runtime_error("Boolean coordinate normalization collapsed a face");
            continue;
        }
        if (!rationalAuthority && (b - a).cross(c - a).squaredNorm() == 0.0)
        {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if (hasExactVertices)
            {
                // Rounded edge subtraction/cross products may erase a real
                // thin face. Test the represented binary64 points exactly;
                // retaining the original indices also retains the SAME
                // rational face and birth facet. No perturbation or repair.
                std::array<CGAL::Gmpq,3> u,v;
                for(int axis=0;axis<3;++axis) {
                    u[axis]=CGAL::Gmpq(b(axis))-CGAL::Gmpq(a(axis));
                    v[axis]=CGAL::Gmpq(c(axis))-CGAL::Gmpq(a(axis));
                }
                bool nonzero=false;
                for(int axis=0;axis<3;++axis) {
                    const int j=(axis+1)%3,k=(axis+2)%3;
                    nonzero |= u[j]*v[k] != u[k]*v[j];
                }
                if (!nonzero)
                    throw std::runtime_error("binary64 Boolean boundary contains an exactly degenerate face");
            }
            else
#endif
                continue; // Preserve the legacy non-tape profile's behavior.
        }
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        if(rationalAuthority)
        {
            std::array<RationalPoint,3> p;
            for(int k=0;k<3;++k)for(int s=0;s<3;++s)
                p[k][s]=verification::parseExactRationalCoordinate((*exactVertices)[uniqueInputs[indices[k]]][s]);
            bool nonzero=false;
            for(int s=0;s<3;++s)
            {
                const int u=(s+1)%3,v=(s+2)%3;
                nonzero |= (p[1][u]-p[0][u])*(p[2][v]-p[0][v]) != (p[1][v]-p[0][v])*(p[2][u]-p[0][u]);
            }
            if(!nonzero)throw std::runtime_error("rational Boolean authority contains a degenerate face");
        }
#endif
        const int minimum = static_cast<int>(
            std::min_element(indices.begin(), indices.end()) - indices.begin()
        );
        std::array<int, 3> rotated = {{
            indices[minimum],
            indices[(minimum + 1) % 3],
            indices[(minimum + 2) % 3]
        }};
        faces.push_back({rotated, (*birthFaces)(row)});
    }
    std::sort(
        faces.begin(),
        faces.end(),
        [](const CanonicalFace& left, const CanonicalFace& right) {
            return std::tie(
                left.indices[0],
                left.indices[1],
                left.indices[2],
                left.birth
            ) < std::tie(
                right.indices[0],
                right.indices[1],
                right.indices[2],
                right.birth
            );
        }
    );

    std::vector<unsigned char> referenced(uniqueVertices.size(), 0U);
    for (const CanonicalFace& face : faces)
    {
        for (const int index : face.indices)
        {
            referenced[static_cast<std::size_t>(index)] = 1U;
        }
    }
    std::vector<int> compactRemap(uniqueVertices.size(), -1);
    int compactVertexCount = 0;
    for (std::size_t index = 0; index < uniqueVertices.size(); ++index)
    {
        if (referenced[index] != 0U)
        {
            compactRemap[index] = compactVertexCount++;
        }
    }
    Eigen::MatrixXd sortedVertices(compactVertexCount, 3);
    std::vector<std::array<std::string, 3>> sortedExactVertices;
    if (hasExactVertices)
    {
        sortedExactVertices.reserve(static_cast<std::size_t>(compactVertexCount));
    }
    int destination = 0;
    for (std::size_t index = 0; index < uniqueVertices.size(); ++index)
    {
        if (referenced[index] != 0U)
        {
            sortedVertices.row(destination++) = uniqueVertices[index];
            if (hasExactVertices)
            {
                sortedExactVertices.push_back((*exactVertices)[uniqueInputs[index]]);
            }
        }
    }
    mesh->vertices = std::move(sortedVertices);
    if (hasExactVertices)
    {
        *exactVertices = std::move(sortedExactVertices);
    }
    mesh->faces.resize(static_cast<Eigen::Index>(faces.size()), 3);
    birthFaces->resize(static_cast<Eigen::Index>(faces.size()));
    for (int row = 0; row < static_cast<int>(faces.size()); ++row)
    {
        mesh->faces.row(row) <<
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[0])],
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[1])],
            compactRemap[static_cast<std::size_t>(
                faces[static_cast<std::size_t>(row)].indices[2])];
        (*birthFaces)(row) = faces[static_cast<std::size_t>(row)].birth;
    }
}


Eigen::Matrix<double, 10, 1> rawMonomials(const Eigen::Vector3d& point)
{
    Eigen::Matrix<double, 10, 1> values;
    values <<
        1.0,
        point.x(), point.y(), point.z(),
        point.x() * point.x(),
        point.y() * point.y(),
        point.z() * point.z(),
        point.y() * point.z(),
        point.x() * point.z(),
        point.x() * point.y();
    return values;
}


Eigen::Matrix<double, 10, 1> coefficientMonomials(
    const Eigen::Vector3d& point
)
{
    Eigen::Matrix<double, 10, 1> values = rawMonomials(point);
    values.tail<3>() *= 2.0;
    return values;
}


struct QuadraturePoint
{
    std::array<double, 3> barycentric;
    double weight;
};


const std::array<QuadraturePoint, 7>& degreeFiveRule()
{
    static const std::array<QuadraturePoint, 7> rule = {{
        {{{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}}, 0.225},
        {{{0.059715871789770, 0.470142064105115, 0.470142064105115}},
            0.132394152788506},
        {{{0.470142064105115, 0.059715871789770, 0.470142064105115}},
            0.132394152788506},
        {{{0.470142064105115, 0.470142064105115, 0.059715871789770}},
            0.132394152788506},
        {{{0.797426985353087, 0.101286507323456, 0.101286507323456}},
            0.125939180544827},
        {{{0.101286507323456, 0.797426985353087, 0.101286507323456}},
            0.125939180544827},
        {{{0.101286507323456, 0.101286507323456, 0.797426985353087}},
            0.125939180544827}
    }};
    return rule;
}


void accumulateInterfaceDerivative(
    const spin::TriangleMesh& material,
    const Eigen::VectorXi& birthFaces,
    int designFaceCount,
    int quadricInterfaceFaceCount,
    const spin::QuadricCoefficients10& quadric,
    double density,
    double minimumGradient,
    Eigen::Matrix<double, 10, 10>* jacobian,
    int* interfaceTriangleCount,
    bool* derivativeDefined
)
{
    *jacobian = Eigen::Matrix<double, 10, 10>::Zero();
    *interfaceTriangleCount = 0;
    *derivativeDefined = true;
    for (Eigen::Index face = 0; face < material.faces.rows(); ++face)
    {
        const int birth = birthFaces(face);
        if (
            birth < designFaceCount ||
            birth >= designFaceCount + quadricInterfaceFaceCount
        )
        {
            continue;
        }
        ++(*interfaceTriangleCount);
        const Eigen::Vector3d a =
            material.vertices.row(material.faces(face, 0));
        const Eigen::Vector3d b =
            material.vertices.row(material.faces(face, 1));
        const Eigen::Vector3d c =
            material.vertices.row(material.faces(face, 2));
        const double area = 0.5 * (b - a).cross(c - a).norm();
        if (!std::isfinite(area) || area <= 0.0)
        {
            *derivativeDefined = false;
            continue;
        }
        for (const QuadraturePoint& sample : degreeFiveRule())
        {
            const Eigen::Vector3d point =
                sample.barycentric[0] * a +
                sample.barycentric[1] * b +
                sample.barycentric[2] * c;
            const double gradientNorm = quadric.gradient(point).norm();
            if (!std::isfinite(gradientNorm) || gradientNorm <= minimumGradient)
            {
                *derivativeDefined = false;
                continue;
            }
            *jacobian += density * area * sample.weight *
                rawMonomials(point) *
                coefficientMonomials(point).transpose() /
                gradientNorm;
        }
    }
    if (*interfaceTriangleCount == 0 || !jacobian->allFinite())
    {
        *derivativeDefined = false;
    }
    const Eigen::Matrix<double, 10, 10> tangentProjection =
        Eigen::Matrix<double, 10, 10>::Identity() -
        quadric.values() * quadric.values().transpose();
    *jacobian *= tangentProjection;
}


std::string topologyHash(
    const spin::TriangleMesh& material,
    const Eigen::VectorXi& birthFaces
)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricMomentTopology-v1;" << material.vertices.rows() << ';'
           << material.faces.rows() << ';';
    for (Eigen::Index face = 0; face < material.faces.rows(); ++face)
    {
        stream << material.faces(face, 0) << ';'
               << material.faces(face, 1) << ';'
               << material.faces(face, 2) << ';'
               << birthFaces(face) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


// Recompute the entire interface vote when the inexpensive double predicate
// cannot decide it. Subtracting two nearly parallel represented edges can
// round away a genuine nonzero area; exact arithmetic uses these same input
// coordinates, not a repaired mesh or a relaxed ambiguity threshold.
bool exactQuadricInterfaceInteriorIsPositive(
    const spin::QuadricCoefficients10& quadric,
    const spin::QuadricTriangulation& triangulation,
    const spin::QuadricCoefficientAuthority10* authority=nullptr)
{
    using Q = CGAL::Gmpq;
    const auto& vertices = triangulation.closedBoundary.vertices;
    const auto& faces = triangulation.closedBoundary.faces;
    const auto quadratic = quadric.quadraticMatrix();
    const auto linear = quadric.linearVector();
    const std::array<std::array<int,3>,3> indices{{{{4,9,8}},{{9,5,7}},{{8,7,6}}}};
    const Q threshold(1e-10);
    int positive = 0, negative = 0;
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        if (triangulation.faceIsCap(face) != 0) continue;
        std::array<Q, 3> a, b, c, u, v, center;
        for (int corner = 0; corner < 3; ++corner)
            if (faces(face, corner) < 0 || faces(face, corner) >= vertices.rows())
                throw std::invalid_argument("quadric-side witness found an invalid interface index");
        for (int axis = 0; axis < 3; ++axis)
        {
            a[axis] = Q(vertices(faces(face, 0), axis));
            b[axis] = Q(vertices(faces(face, 1), axis));
            c[axis] = Q(vertices(faces(face, 2), axis));
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
            if(!triangulation.exactVertices.empty()) {
                a[axis]=verification::parseExactRationalCoordinate(triangulation.exactVertices.at(faces(face,0))[axis]);
                b[axis]=verification::parseExactRationalCoordinate(triangulation.exactVertices.at(faces(face,1))[axis]);
                c[axis]=verification::parseExactRationalCoordinate(triangulation.exactVertices.at(faces(face,2))[axis]);
            }
#endif
            u[axis] = b[axis] - a[axis];
            v[axis] = c[axis] - a[axis];
            center[axis] = (a[axis] + b[axis] + c[axis]) / Q(3);
        }
        const std::array<Q, 3> area{
            u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0]};
        Q areaSquared(0), gradientSquared(0), dot(0);
        for (int axis = 0; axis < 3; ++axis)
        {
            Q gradient=authority ? Q(authority->geometry()[axis+1]) : Q(linear(axis));
            for (int other = 0; other < 3; ++other)
                gradient += Q(2) * (authority ? Q(authority->geometry()[indices[axis][other]]) : Q(quadratic(axis, other))) * center[other];
            areaSquared += area[axis] * area[axis];
            gradientSquared += gradient * gradient;
            dot += area[axis] * gradient;
        }
        if (areaSquared == Q(0) || gradientSquared == Q(0))
            throw std::runtime_error("quadric-side witness is degenerate on an interface face");
        if (dot*dot <= threshold*threshold*areaSquared*gradientSquared)
            throw std::runtime_error("quadric-side witness is ambiguous on an interface face");
        if (dot > Q(0)) ++positive;
        else ++negative;
    }
    if (positive + negative == 0)
        throw std::runtime_error("quadric-side witness has no non-cap interface faces");
    if (positive != 0 && negative != 0)
    {
        std::ostringstream diagnostic;
        diagnostic << "outward quadric interface faces disagree on the interior sign"
            << ": positive_exterior_faces=" << positive
            << " negative_exterior_faces=" << negative;
        throw std::runtime_error(diagnostic.str());
    }
    return negative != 0;
}

} // namespace


namespace spin
{

namespace {
thread_local ScopedQuadricGeometryObserver::Observer geometryObserver;
void observeGeometry(const char* stage, const QuadricMomentEvaluation& evaluation)
{
    if (geometryObserver) geometryObserver(stage, evaluation);
}
}

ScopedQuadricGeometryObserver::ScopedQuadricGeometryObserver(Observer observer)
    : previous_(std::move(geometryObserver))
{
    geometryObserver = std::move(observer);
}

ScopedQuadricGeometryObserver::~ScopedQuadricGeometryObserver()
{
    geometryObserver = std::move(previous_);
}

bool quadricSelectorHasStrictInactiveCentroidMargin(
    const QuadricCoefficients10& quadric,
    const QuadricTriangulation& triangulation,
    bool selectedPositiveRegion)
{
    using Q = CGAL::Gmpq;
    using Point = std::array<Q,3>;
    const auto& vertices = triangulation.closedBoundary.vertices;
    const auto& faces = triangulation.closedBoundary.faces;
    if (vertices.cols()!=3 || !vertices.allFinite() || faces.cols()!=3 ||
        faces.rows()==0 || triangulation.faceIsCap.size()!=faces.rows())
        return false;
    Q volumeSix(0);
    Point firstTwentyFour{Q(0),Q(0),Q(0)};
    for (Eigen::Index face=0;face<faces.rows();++face)
    {
        std::array<Point,3> p;
        for (int corner=0;corner<3;++corner)
        {
            const int index=faces(face,corner);
            if(index<0 || index>=vertices.rows()) return false;
            for(int axis=0;axis<3;++axis) {
                p[corner][axis]=Q(vertices(index,axis));
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                if(!triangulation.exactVertices.empty())
                    p[corner][axis]=verification::parseExactRationalCoordinate(triangulation.exactVertices.at(index)[axis]);
#endif
            }
        }
        const Q determinant=p[0][0]*(p[1][1]*p[2][2]-p[1][2]*p[2][1])+
            p[0][1]*(p[1][2]*p[2][0]-p[1][0]*p[2][2])+
            p[0][2]*(p[1][0]*p[2][1]-p[1][1]*p[2][0]);
        volumeSix+=determinant;
        for(int axis=0;axis<3;++axis)
            firstTwentyFour[axis]+=determinant*(p[0][axis]+p[1][axis]+p[2][axis]);
    }
    if(volumeSix<=Q(6)*Q(1e-15)) return false;
    Point center;
    for(int axis=0;axis<3;++axis)
        center[axis]=firstTwentyFour[axis]/(Q(4)*volumeSix);
    const auto quadratic=quadric.quadraticMatrix();
    const auto linear=quadric.linearVector();
    Q value(quadric.constant());
    for(int axis=0;axis<3;++axis)
    {
        value+=Q(linear(axis))*center[axis];
        for(int other=0;other<3;++other)
            value+=center[axis]*Q(quadratic(axis,other))*center[other];
    }
    // The actual threshold is 1e-10*max(1,scale), hence this STRICT
    // lower-bound test proves inactivity regardless of winding or scale.
    // q(center)=0 is inside this open band, not its availability boundary.
    const Q threshold(1e-10);
    if(value<=-threshold || value>=threshold) return false;
    try
    {
        // Check EVERY non-cap face exactly, including the unchanged squared
        // 1e-10 alignment margin; no majority vote or sampled witness.
        return exactQuadricInterfaceInteriorIsPositive(quadric,triangulation)==
            selectedPositiveRegion;
    }
    catch(const std::runtime_error&)
    {
        return false;
    }
}

bool quadricSolidInteriorIsPositiveRegion(
    const QuadricCoefficients10& quadric,
    const QuadricTriangulation& triangulation)
{
    if(!triangulation.exactVertices.empty())
        return exactQuadricInterfaceInteriorIsPositive(quadric,triangulation);
    const auto& vertices = triangulation.closedBoundary.vertices;
    const auto& faces = triangulation.closedBoundary.faces;
    if (vertices.cols() != 3 || !vertices.allFinite() || faces.cols() != 3 ||
        faces.rows() == 0 || triangulation.faceIsCap.size() != faces.rows())
        throw std::invalid_argument(
            "quadric-side witness requires finite source geometry and face-cap provenance");

    // A projective two-sheet source can have separated closed components.
    // Its aggregate centroid can lie between them, and remote ray chords
    // need not align with grad(q) at the chord centroid. Require a strict
    // interior witness in EVERY component, with unanimous implicit sign.
    // This does not reorient faces or select a side by a face majority.
    if (triangulation.paperProjectiveHyperboloid)
    {
        Eigen::MatrixXi unusedOrientedFaces;
        Eigen::VectorXi labels;
        igl::bfs_orient(faces, unusedOrientedFaces, labels);
        const int count = labels.size() == faces.rows() ? labels.maxCoeff() + 1 : 0;
        bool allWitnessed = count > 1;
        std::optional<bool> commonSign;
        for (int component = 0; allWitnessed && component < count; ++component)
        {
            const Eigen::Index size = (labels.array() == component).count();
            Eigen::MatrixXi componentFaces(size, 3);
            Eigen::Index row = 0;
            for (Eigen::Index face = 0; face < faces.rows(); ++face)
                if (labels(face) == component) componentFaces.row(row++) = faces.row(face);
            const MassMoments componentMoments = computeSurfaceMassMoments(vertices, componentFaces, 1.0);
            Eigen::Vector3d point;
            if (!componentMoments.tryCenterOfMass(point, 1e-15)) { allWitnessed = false; break; }
            const double winding = igl::winding_number(vertices, componentFaces, point);
            const double value = quadric.evaluate(point);
            const Eigen::Vector3d absolute = point.cwiseAbs();
            const double scale = std::abs(quadric.constant()) +
                quadric.linearVector().cwiseAbs().dot(absolute) +
                absolute.dot(quadric.quadraticMatrix().cwiseAbs() * absolute);
            if (!std::isfinite(winding) || winding <= 0.5 + 1e-6 ||
                !std::isfinite(value) || std::abs(value) <= 1e-10 * std::max(1.0, scale))
            { allWitnessed = false; break; }
            const bool sign = value > 0.0;
            if (commonSign && *commonSign != sign)
                throw std::runtime_error("quadric component interior witnesses disagree on the interior sign");
            commonSign = sign;
        }
        if (allWitnessed && commonSign) return *commonSign;
    }

    // Preserve the fast historical rule only when the volume centroid is
    // actually certified inside the outward source. A non-convex solid's
    // centroid is not, by itself, an interior point.
    const MassMoments moments = computeSurfaceMassMoments(vertices, faces, 1.0);
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    if (moments.tryCenterOfMass(centroid, 1e-15))
    {
        const double winding = igl::winding_number(vertices, faces, centroid);
        const double value = quadric.evaluate(centroid);
        const Eigen::Vector3d absolute = centroid.cwiseAbs();
        const Eigen::Matrix3d absoluteQuadratic = quadric.quadraticMatrix().cwiseAbs();
        const double scale = std::abs(quadric.constant()) +
            quadric.linearVector().cwiseAbs().dot(absolute) +
            absolute.dot(absoluteQuadratic * absolute);
        const double threshold = 1e-10 * std::max(1.0, scale);
        if (std::isfinite(winding) && std::abs(winding) > 0.5 + 1e-6 &&
            std::isfinite(value) && std::abs(value) > threshold)
            return value > 0.0;
    }

    int positiveExteriorVotes = 0;
    int negativeExteriorVotes = 0;
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        if (triangulation.faceIsCap(face) != 0) continue;
        const int ia = faces(face, 0), ib = faces(face, 1), ic = faces(face, 2);
        if (ia < 0 || ib < 0 || ic < 0 || ia >= vertices.rows() ||
            ib >= vertices.rows() || ic >= vertices.rows())
            throw std::invalid_argument("quadric-side witness found an invalid interface index");
        const Eigen::Vector3d a = vertices.row(ia).transpose();
        const Eigen::Vector3d b = vertices.row(ib).transpose();
        const Eigen::Vector3d c = vertices.row(ic).transpose();
        const Eigen::Vector3d areaVector = (b - a).cross(c - a);
        const Eigen::Vector3d gradient = quadric.gradient((a + b + c) / 3.0);
        const double denominator = areaVector.norm() * gradient.norm();
        if (!std::isfinite(denominator) || denominator <= 0.0)
            return exactQuadricInterfaceInteriorIsPositive(quadric, triangulation);
        const double alignment = areaVector.dot(gradient) / denominator;
        if (!std::isfinite(alignment) || std::abs(alignment) <= 1e-10)
            return exactQuadricInterfaceInteriorIsPositive(quadric, triangulation);
        if (alignment > 0.0) ++positiveExteriorVotes;
        else ++negativeExteriorVotes;
    }
    if (positiveExteriorVotes + negativeExteriorVotes == 0)
        throw std::runtime_error("quadric-side witness has no non-cap interface faces");
    if (positiveExteriorVotes != 0 && negativeExteriorVotes != 0)
        return exactQuadricInterfaceInteriorIsPositive(quadric, triangulation);

    // The source interface is oriented outward. Positive alignment means q
    // increases when leaving the solid, so the interior lies on q<0.
    return negativeExteriorVotes != 0;
}

std::string quadricBooleanBackendIdentity()
{
    return std::string("libigl::copyleft::cgal::mesh_boolean;libigl_source=") +
        MAKE_IT_STAND_LIBIGL_SOURCE_ID;
}

const char* quadricMomentDerivativeMethodToken(QuadricMomentDerivativeMethod method) noexcept
{
    switch (method)
    {
    case QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape:
        return "legacy-continuous-implicit-shape-v1";
    case QuadricMomentDerivativeMethod::DiscreteFullRankQuadric:
        return "discrete-fullrank-quadric-boolean-transport-v2";
    case QuadricMomentDerivativeMethod::DiscreteReducedStratum:
        return "discrete-rank-deficient-stratum-boolean-transport-v1";
    }
    return "unknown";
}

void validateQuadricMomentEvaluatorSettings(const QuadricMomentEvaluatorSettings& settings)
{
#if !MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(settings.triangulation.rationalSourceCoordinates)
        throw std::invalid_argument("rational source coordinates require exact geometry checks");
    if (settings.exactBinary64SelfIntersection)
        throw std::invalid_argument("exact binary64 self-intersection policy requires exact geometry checks");
    if (settings.exactBinary64SourceBoundary)
        throw std::invalid_argument("exact binary64 source boundary policy requires exact geometry checks");
    if (settings.exactBinary64MaterialBoundary)
        throw std::invalid_argument("exact binary64 material boundary policy requires exact geometry checks");
#endif
    if (settings.emptyRegionPolicy != QuadricEmptyRegionPolicy::RejectEmpty &&
        settings.emptyRegionPolicy != QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1)
        throw std::invalid_argument("unknown quadric empty-region policy");
    if (settings.derivativeMethod != QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape &&
        settings.derivativeMethod != QuadricMomentDerivativeMethod::DiscreteFullRankQuadric &&
        settings.derivativeMethod != QuadricMomentDerivativeMethod::DiscreteReducedStratum)
        throw std::invalid_argument("unknown quadric derivative method");
    if (settings.capRadiusPolicy != QuadricCapRadiusPolicy::RejectUncertified &&
        settings.capRadiusPolicy != QuadricCapRadiusPolicy::GrowToCertifiedV1)
        throw std::invalid_argument("unknown quadric cap radius policy");
    if (!std::isfinite(settings.minimumInterfaceGradient) || settings.minimumInterfaceGradient <= 0.0)
        throw std::invalid_argument("minimum interface gradient must be finite and positive");
    const auto& t = settings.triangulation;
    if(t.compactRationalSourceCoordinates && !t.rationalSourceCoordinates)
        throw std::invalid_argument("compact rational source requires rational source coordinates");
    if(t.rationalSourceCoordinates &&
       settings.derivativeMethod!=QuadricMomentDerivativeMethod::DiscreteFullRankQuadric)
        throw std::invalid_argument("rational source coordinates require discrete full-rank derivatives");
    if (t.icosphereSubdivisions < 0 || t.icosphereSubdivisions > 7 ||
        (t.maximumIcosphereSubdivisions != -1 &&
         (t.maximumIcosphereSubdivisions < t.icosphereSubdivisions ||
          t.maximumIcosphereSubdivisions > 7)) ||
        t.paraboloidGridResolution < 1 || t.cylinderAngularResolution < 8 ||
        t.hyperbolicCylinderGridResolution < 2 ||
        t.parabolicCylinderGridResolution < 1 ||
        t.maximumRadiusDoublings < 0 || !std::isfinite(t.truncationRadius) || t.truncationRadius < 0.0)
        throw std::invalid_argument("invalid quadric triangulation settings");
}


std::string computeQuadricEvaluatorRequestHash(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& settings)
{
    validateQuadricMomentEvaluatorSettings(settings);
    const auto& t = settings.triangulation;
    const std::string base = computeSpinDesignRequestHash(request, settings.exactBinary64SourceBoundary);
    if (settings.derivativeMethod == QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape &&
        settings.emptyRegionPolicy == QuadricEmptyRegionPolicy::RejectEmpty &&
        settings.capRadiusPolicy == QuadricCapRadiusPolicy::RejectUncertified &&
        !settings.exactBinary64SelfIntersection &&
        !settings.exactBinary64SourceBoundary &&
        !settings.exactBinary64MaterialBoundary &&
        settings.minimumInterfaceGradient == 1e-12 && t.icosphereSubdivisions == 3 &&
        t.maximumIcosphereSubdivisions == -1 &&
        t.paraboloidGridResolution == 24 && t.cylinderAngularResolution == 96 &&
        t.hyperbolicCylinderGridResolution == 24 &&
        t.parabolicCylinderGridResolution == 24 &&
        t.truncationRadius == 0.0 && t.maximumRadiusDoublings == 8 &&
        !t.retainFixedNormalProvenance && !t.paperProjectiveHyperboloid &&
        !t.algebraicFiniteMap && !t.fixedNormalParaboloid && !t.rationalSourceCoordinates)
        return base;
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    identity << "QuadricEvaluatorRequest-v2;" << base << ';'
        << static_cast<int>(settings.derivativeMethod) << ';'
        << static_cast<int>(settings.emptyRegionPolicy) << ';'
        << static_cast<int>(settings.capRadiusPolicy) << ';'
        << t.icosphereSubdivisions << ';'
        << t.maximumIcosphereSubdivisions << ';'
        << t.paraboloidGridResolution << ';'
        << t.cylinderAngularResolution << ';'
        << t.hyperbolicCylinderGridResolution << ';'
        << t.parabolicCylinderGridResolution << ';'
        << t.maximumRadiusDoublings << ';'
        << t.retainFixedNormalProvenance << ';' << std::hexfloat
        << t.truncationRadius << ';' << settings.minimumInterfaceGradient << ';';
    // Preserve every existing legacy settings identity byte-for-byte.
    if (settings.emptyRegionPolicy == QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1)
        identity << "empty-region=certify-strictly-separated-v1;";
    if (settings.capRadiusPolicy == QuadricCapRadiusPolicy::GrowToCertifiedV1)
        identity << "cap-radius=grow-to-certified-v1;";
    if (t.paperProjectiveHyperboloid)
        identity << "hyperboloid-source=paper-projective-triangles-v1;";
    if (t.algebraicFiniteMap)
        identity << "finite-map=algebraic-v1;";
    if (t.fixedNormalParaboloid)
        identity << "paraboloid-source=finite-normal-paraboloid-v1;";
    if (t.rationalSourceCoordinates)
        identity << "source-coordinates=rational-enclosed-v1;absolute-width=2^-384;";
    if (t.compactRationalSourceCoordinates)
        identity << "representative=enclosed-dyadic-1024-v1;";
    if (settings.exactBinary64SelfIntersection)
        identity << "material-self-intersection=original-binary64-exact-indexed-simplex-v1;";
    if (settings.exactBinary64SourceBoundary)
        identity << "canonical-source-boundary=original-binary64-exact-outward-positive-v1;";
    if (settings.exactBinary64MaterialBoundary)
        identity << "solver-material-boundary=original-binary64-exact-signed-shell-v2-connected-embedding-or-nef;";
    if (settings.exactBinary64SourceBoundary && settings.exactBinary64MaterialBoundary &&
        settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteFullRankQuadric)
        identity << "motion-coordinates=same-boolean-rational-before-assignment-v1;";
    return manufacturing::computeContentSha256(identity.str());
}

manufacturing::BoundaryAudit auditQuadricMaterialBoundary(
    const TriangleMesh& mesh, bool exactBinary64SelfIntersection)
{
    return auditQuadricMaterialBoundary(mesh,exactBinary64SelfIntersection,false);
}

static manufacturing::BoundaryAudit auditQuadricMaterialBoundaryImpl(
    const TriangleMesh& mesh, bool exactBinary64SelfIntersection, bool exactBinary64MaterialBoundary)
{
    auto result=manufacturing::auditTriangleBoundary(mesh.vertices,mesh.faces,
        manufacturing::BoundaryOrientation::AnyNonzero);
    if (exactBinary64MaterialBoundary)
    {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
        // Fail closed even when the old heuristic happens to accept a mesh.
        // Exact nonzero areas are geometric facts, NOT a conditioning or
        // printability certificate. Never alter coordinates or triangulation.
        result.nondegenerateFaces=false;
        if (!result.matrixShapeValid || !result.finiteVertices || !result.faceIndicesValid)
            return result;
        // An embedded connected closed orientable boundary of positive exact
        // volume bounds one finite regular solid (Jordan-Brouwer). The exact
        // indexed-simplex certificate already proves every required premise;
        // building a second Nef representation adds no region information.
        // Multiple boundary components still require signed-shell nesting.
        const auto single=verification::auditExactTriangleBoundary(mesh.vertices,mesh.faces,
            verification::ExactBoundaryOrientation::AnyNonzero);
        const auto& connected=single.facts();
        if(connected.boundaryValid && connected.componentVolumeSigns.size()==1 &&
           connected.componentVolumeSigns.front()>0) {
            result.nondegenerateFaces=true;result.degenerateFaceCount=0;
            result.edgeManifold=connected.edgeManifold;result.vertexManifold=connected.vertexManifold;
            result.closed=connected.closed;result.consistentlyOriented=connected.consistentlyOriented;
            result.componentOrientationsValid=connected.componentOrientationsValid;
            result.selfIntersections.completed=connected.intersectionCompleted;
            result.selfIntersections.intersectionFree=connected.intersectionFree;
            result.selfIntersections.intersectingPairCount=connected.intersectingPairs.size();
            result.selfIntersections.narrowPhasePairCount=connected.exactPairTests;
            result.selfIntersections.firstFace=-1;result.selfIntersections.secondFace=-1;
            return result;
        }
        const auto region=verification::makeExactMaterialRegion(mesh.vertices,mesh.faces);
        const auto& facts=region.boundaryFacts();
        result.nondegenerateFaces=region.facts().completed && region.facts().valid &&
            facts.boundaryValid && facts.degenerateFaces.empty();
        result.degenerateFaceCount=static_cast<int>(facts.degenerateFaces.size());
        result.edgeManifold=facts.edgeManifold;
        result.vertexManifold=facts.vertexManifold;
        result.closed=facts.closed;
        result.consistentlyOriented=facts.consistentlyOriented;
        result.componentOrientationsValid=facts.componentOrientationsValid;
        result.selfIntersections.completed=facts.intersectionCompleted;
        result.selfIntersections.intersectionFree=facts.intersectionFree;
        result.selfIntersections.intersectingPairCount=facts.intersectingPairs.size();
        result.selfIntersections.narrowPhasePairCount=facts.exactPairTests;
        result.selfIntersections.firstFace=facts.intersectingPairs.empty()?-1:facts.intersectingPairs.front()[0];
        result.selfIntersections.secondFace=facts.intersectingPairs.empty()?-1:facts.intersectingPairs.front()[1];
        return result;
#else
        throw std::runtime_error("exact binary64 material boundary policy unavailable");
#endif
    }
    if (!exactBinary64SelfIntersection || !result.matrixShapeValid ||
        !result.finiteVertices || !result.faceIndicesValid || !result.nondegenerateFaces)
        return result;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    // Only the self-intersection predicate changes. Never translate/trim the
    // queried simplices or override existing area, volume or topology gates.
    const auto exact=verification::auditExactTriangleBoundary(mesh.vertices,mesh.faces,
        verification::ExactBoundaryOrientation::AnyNonzero);
    const auto& facts=exact.facts();
    result.selfIntersections.completed=facts.intersectionCompleted;
    result.selfIntersections.intersectionFree=facts.intersectionFree;
    result.selfIntersections.intersectingPairCount=facts.intersectingPairs.size();
    result.selfIntersections.narrowPhasePairCount=facts.exactPairTests;
    result.selfIntersections.firstFace=facts.intersectingPairs.empty()?-1:facts.intersectingPairs.front()[0];
    result.selfIntersections.secondFace=facts.intersectingPairs.empty()?-1:facts.intersectingPairs.front()[1];
#else
    throw std::runtime_error("exact binary64 self-intersection policy unavailable");
#endif
    return result;
}

manufacturing::BoundaryAudit auditQuadricMaterialBoundary(
    const TriangleMesh& mesh, bool exactBinary64SelfIntersection, bool exactBinary64MaterialBoundary)
{
    struct ProofEntry {
        std::string indexedBoundary;
        manufacturing::BoundaryAudit proof;
    };
    static thread_local std::array<ProofEntry,2> completedProofs;
    static thread_local std::size_t replacement=0;
    std::ostringstream key;key.imbue(std::locale::classic());
    key<<"QuadricIndexedBoundaryProof-v1;"<<exactBinary64SelfIntersection<<';'
       <<exactBinary64MaterialBoundary<<';'<<mesh.vertices.rows()<<';'<<mesh.vertices.cols()<<';'
       <<mesh.faces.rows()<<';'<<mesh.faces.cols()<<';'<<std::hexfloat;
    for(Eigen::Index row=0;row<mesh.vertices.rows();++row)
        for(Eigen::Index axis=0;axis<mesh.vertices.cols();++axis)key<<mesh.vertices(row,axis)<<';';
    for(Eigen::Index row=0;row<mesh.faces.rows();++row)
        for(Eigen::Index corner=0;corner<mesh.faces.cols();++corner)key<<mesh.faces(row,corner)<<';';
    const std::string indexedBoundary=key.str();
    for(const auto& entry:completedProofs)
        if(entry.indexedBoundary==indexedBoundary)return entry.proof;
    const auto proof=auditQuadricMaterialBoundaryImpl(mesh,exactBinary64SelfIntersection,exactBinary64MaterialBoundary);
    if(proof.strictBoundaryValid() && proof.signedVolume>0.) {
        completedProofs[replacement]={indexedBoundary,proof};
        replacement=(replacement+1)%completedProofs.size();
    }
    return proof;
}

static QuadricMomentEvaluation evaluateQuadricMaterialMomentsImpl(
    const SpinDesignRequest& request,
    const QuadricCoefficients10& quadricTarget,
    const QuadricMomentEvaluatorSettings& settings,
    const QuadricCoefficientAuthority10* authority,
    QuadricAuthorityMomentAttachment* authorityDerivative=nullptr
)
{
    const auto evaluationStart = std::chrono::steady_clock::now();
    validateSpinDesignRequest(request, settings.exactBinary64SourceBoundary);
    validateQuadricMomentEvaluatorSettings(settings);
    const auto stage=[&](const char* name){if(settings.diagnosticStageObserver)settings.diagnosticStageObserver(name);};

    TriangleMesh designTarget =
        toTarget(request.designDomain, request.targetFrame);
    QuadricMomentEvaluation result;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(settings.triangulation.rationalSourceCoordinates) {
        result.designExactVerticesTarget=rationalTargetVertices(request.designDomain,request.targetFrame);
        assignRationalVertices(result.designExactVerticesTarget,&designTarget.vertices);
    }
#endif
    result.derivativeMethod = settings.derivativeMethod;
    result.booleanBackendIdentity = quadricBooleanBackendIdentity();
    QuadricTriangulationSettings triangulationSettings = settings.triangulation;
    if (settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteFullRankQuadric ||
        settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteReducedStratum)
        triangulationSettings.retainFixedNormalProvenance = true;
    int capGrowths = 0;
    stage("source-construction-begin");
    for (;;)
    {
        result.quadricTriangulation = authority
            ? triangulateQuadric(*authority,designTarget.vertices,triangulationSettings)
            : triangulateQuadric(quadricTarget,designTarget.vertices,triangulationSettings);
        try
        {
            requireTruncationCapsOutsideDesignDomain(
                designTarget.vertices, result.quadricTriangulation,
                settings.triangulation.rationalSourceCoordinates?&result.designExactVerticesTarget:nullptr);
            break;
        }
        catch (const std::runtime_error& error)
        {
            if (settings.capRadiusPolicy != QuadricCapRadiusPolicy::GrowToCertifiedV1 ||
                result.quadricTriangulation.capFaceCount == 0 ||
                capGrowths >= settings.triangulation.maximumRadiusDoublings ||
                !std::isfinite(result.quadricTriangulation.truncationRadius) ||
                result.quadricTriangulation.truncationRadius <= 0.0)
                throw;
            triangulationSettings.truncationRadius =
                result.quadricTriangulation.truncationRadius * 2.0;
            triangulationSettings.maximumRadiusDoublings = 0;
            ++capGrowths;
            (void)error;
        }
    }

    observeGeometry("source-before-side", result);
    result.quadricSolidIsPositiveRegion = authority
        ? exactQuadricInterfaceInteriorIsPositive(quadricTarget,result.quadricTriangulation,authority)
        : quadricSolidInteriorIsPositiveRegion(quadricTarget, result.quadricTriangulation);
    stage("source-construction-end");

    Eigen::VectorXi birthFaces;
    const igl::MeshBooleanType operation =
        result.quadricSolidIsPositiveRegion
        ? igl::MESH_BOOLEAN_TYPE_INTERSECT
        : igl::MESH_BOOLEAN_TYPE_MINUS;
    bool booleanSucceeded;
    stage("boolean-begin");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(settings.triangulation.rationalSourceCoordinates) {
        auto exact=verification::computeExactLibiglBooleanConstructionFromRationals(
            result.designExactVerticesTarget,designTarget.faces,
            result.quadricTriangulation.exactVertices,result.quadricTriangulation.closedBoundary.faces,
            result.quadricSolidIsPositiveRegion?verification::ExactMaterialBooleanOperation::Intersection:
                verification::ExactMaterialBooleanOperation::Difference);
        booleanSucceeded=exact.succeeded;
        result.materialBoundaryTarget.faces=std::move(exact.faces);
        birthFaces=std::move(exact.birthFaces);
        result.materialExactVerticesTarget=std::move(exact.vertices);
        assignRationalVertices(result.materialExactVerticesTarget,&result.materialBoundaryTarget.vertices);
    } else if(settings.exactBinary64SourceBoundary && settings.exactBinary64MaterialBoundary &&
       settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteFullRankQuadric) {
        // ONE Boolean creates both the actual binary64 boundary and its motion
        // tape. Never replay a separate construction or alter triangulation.
        auto exact=verification::computeExactLibiglBooleanConstruction(
            designTarget.vertices,designTarget.faces,
            result.quadricTriangulation.closedBoundary.vertices,
            result.quadricTriangulation.closedBoundary.faces,
            result.quadricSolidIsPositiveRegion?verification::ExactMaterialBooleanOperation::Intersection:
                verification::ExactMaterialBooleanOperation::Difference);
        booleanSucceeded=exact.succeeded;
        result.materialBoundaryTarget.faces=std::move(exact.faces);
        birthFaces=std::move(exact.birthFaces);
        result.materialExactVerticesTarget=std::move(exact.vertices);
        result.materialBoundaryTarget.vertices.resize(result.materialExactVerticesTarget.size(),3);
        for(std::size_t v=0;v<result.materialExactVerticesTarget.size();++v)for(int axis=0;axis<3;++axis) {
            const auto rational=verification::parseExactRationalCoordinate(result.materialExactVerticesTarget[v][axis]);
            igl::copyleft::cgal::assign_scalar(CGAL::Epeck::FT(rational),result.materialBoundaryTarget.vertices(v,axis));
        }
    } else
#endif
    booleanSucceeded = igl::copyleft::cgal::mesh_boolean(
        designTarget.vertices,
        designTarget.faces,
        result.quadricTriangulation.closedBoundary.vertices,
        result.quadricTriangulation.closedBoundary.faces,
        operation,
        result.materialBoundaryTarget.vertices,
        result.materialBoundaryTarget.faces,
        birthFaces
    );
    stage("boolean-end");
    if (geometryObserver) {
        result.materialBirthFaces = birthFaces;
        observeGeometry("boolean-before-normalization", result);
    }
    if (!booleanSucceeded)
    {
        // In the pinned libigl mesh_boolean implementation, false is returned
        // when winding-number propagation rejects the stitched arrangement
        // produced by remesh_self_intersections. Keep this separate from an
        // empty Boolean result: it is a geometry-kernel failure, not evidence
        // that the requested material region is empty.
        std::ostringstream diagnostic;
        diagnostic << "BOOLEAN_ARRANGEMENT_PWN_FAILURE: libigl rejected the "
            << "stitched Boolean arrangement during winding-number propagation; operation="
            << (operation == igl::MESH_BOOLEAN_TYPE_INTERSECT ? "INTERSECT" : "MINUS")
            << " design_vertices=" << designTarget.vertices.rows()
            << " design_faces=" << designTarget.faces.rows()
            << " source_vertices=" << result.quadricTriangulation.closedBoundary.vertices.rows()
            << " source_faces=" << result.quadricTriangulation.closedBoundary.faces.rows()
            << " output_vertices=" << result.materialBoundaryTarget.vertices.rows()
            << " output_faces=" << result.materialBoundaryTarget.faces.rows()
            << " birth_faces=" << birthFaces.size();
        throw std::runtime_error(diagnostic.str());
    }
    if (
        result.materialBoundaryTarget.faces.rows() == 0 ||
        birthFaces.size() != result.materialBoundaryTarget.faces.rows()
    )
    {
        if (settings.emptyRegionPolicy == QuadricEmptyRegionPolicy::CertifyStrictlySeparatedV1)
        {
            const bool emptyBoolean = booleanSucceeded &&
                result.materialBoundaryTarget.vertices.rows() == 0 &&
                result.materialBoundaryTarget.faces.rows() == 0 && birthFaces.size() == 0;
            if (emptyBoolean && certifyNegativeOnDesignBox(
                    designTarget.vertices, quadricTarget, &result.emptyRegionCertificate))
            {
                const Eigen::Matrix3d modelToTarget = request.targetFrame.axes.transpose();
                result.totalMomentsTarget = request.fixedShellMoments.transformed(
                    modelToTarget, -modelToTarget * request.targetFrame.origin);
                if (!result.totalMomentsTarget.isFinite() || result.totalMomentsTarget.mass() <= 0.0)
                    throw std::runtime_error("certified empty variable region requires positive finite total fixed-shell mass");
                result.materialRegion = QuadricMaterialRegion::CertifiedEmpty;
                result.materialBoundaryTarget.vertices.resize(0,3);
                result.materialBoundaryTarget.faces.resize(0,3);
                result.materialBoundary = result.materialBoundaryTarget;
                result.derivativeDefined = settings.derivativeMethod ==
                    QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape;
                result.derivativeDiagnostic = result.derivativeDefined
                    ? "strict analytic beta separation; locally constant zero analytic variable moments"
                    : "empty analytic region certified; discrete boundary separation unverified";
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
                if (!result.derivativeDefined && result.quadricTriangulation.type == QuadricSurfaceType::Ellipsoid)
                {
                    try
                    {
                        const auto a = verification::makeExactMaterialRegion(designTarget.vertices,designTarget.faces);
                        const auto& boundary = result.quadricTriangulation.closedBoundary;
                        const auto b = verification::makeExactMaterialRegion(boundary.vertices,boundary.faces);
                        const auto separation = verification::checkExactBoundaryClearance(a,b,0.0);
                        result.derivativeDefined = separation.completed &&
                            separation.applicable.value_or(false) && separation.contact.has_value() &&
                            !*separation.contact && separation.squaredDistanceExact.has_value();
                        result.emptyRegionCertificate += ";discrete-clearance=" + separation.clearanceDigest;
                        if (separation.squaredDistanceExact)
                            result.emptyRegionCertificate += ";squared=" + *separation.squaredDistanceExact;
                        result.derivativeDiagnostic = result.derivativeDefined
                            ? "strict analytic and exact discrete boundary separation; locally constant zero discrete moments"
                            : "empty region certified; discrete boundary contact or incomplete clearance";
                    }
                    catch (const std::exception&)
                    {
                        // A failed clearance proof never grants a derivative.
                        result.derivativeDefined = false;
                        result.derivativeDiagnostic = "empty region certified; discrete clearance proof failed";
                    }
                }
#endif
                result.emptyRegionCertificate += result.derivativeDefined ? ";derivative=defined" : ";derivative=undefined";
                const std::string binding = computeQuadricEvaluatorRequestHash(request,settings) + ';' +
                    quadricTarget.hash() + ';' + result.emptyRegionCertificate;
                result.geometryHash = manufacturing::computeContentSha256("CertifiedEmptyVariableGeometry-v1;" + binding);
                result.topologyHash = manufacturing::computeContentSha256("CertifiedEmptyVariableTopology-v1;" +
                    meshHash(designTarget,"DesignTarget-v1"));
                result.provenanceTopologyHash = result.topologyHash;
                result.derivativeScientificHash = manufacturing::computeContentSha256(
                    "CertifiedEmptyMomentDifferential-v2;" + result.booleanBackendIdentity + ';' +
                    binding + ';' + result.geometryHash);
                result.residualEvaluationSeconds = std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-evaluationStart).count();
                return result;
            }
            std::ostringstream diagnostic;
            diagnostic << "UNRESOLVED_GEOMETRY_FAILURE: boolean_success="
                << booleanSucceeded << " vertices="
                << result.materialBoundaryTarget.vertices.rows() << " faces="
                << result.materialBoundaryTarget.faces.rows() << " birth_faces="
                << birthFaces.size();
            throw std::runtime_error(diagnostic.str());
        }
        throw std::runtime_error(
            "quadric/design-domain Boolean produced no auditable material"
        );
    }
    canonicalizeBooleanMesh(&result.materialBoundaryTarget, &birthFaces,
        &result.materialExactVerticesTarget,settings.triangulation.rationalSourceCoordinates);
    result.materialBirthFaces=birthFaces;
    observeGeometry("boolean-after-normalization", result);
    // Keep the immutable proof from this actual Boolean alive through motion
    // transport. Reuse is permitted only after an exact indexed-input match.
    std::optional<verification::ExactMaterialRegion> rationalMaterialAuthority;
    stage("material-authority-audit-begin");
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(settings.triangulation.rationalSourceCoordinates) {
        rationalMaterialAuthority.emplace(verification::makeExactMaterialRegionFromRationals(
            result.materialExactVerticesTarget,result.materialBoundaryTarget.faces));
        if(!rationalMaterialAuthority->facts().completed || !rationalMaterialAuthority->facts().valid)
            throw std::runtime_error("rational material authority is not a valid signed-shell boundary");
    } else
#endif
    {
    const manufacturing::BoundaryAudit materialAudit =
        auditQuadricMaterialBoundary(result.materialBoundaryTarget,
            settings.exactBinary64SelfIntersection,settings.exactBinary64MaterialBoundary);
    if (!materialAudit.strictBoundaryValid() || materialAudit.signedVolume <= 0.0)
    {
        double minimumSquaredDoubleArea =
            std::numeric_limits<double>::infinity();
        Eigen::Index minimumAreaFace = -1;
        for (Eigen::Index face = 0;
             face < result.materialBoundaryTarget.faces.rows(); ++face)
        {
            const Eigen::Vector3d a = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 0));
            const Eigen::Vector3d b = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 1));
            const Eigen::Vector3d c = result.materialBoundaryTarget.vertices.row(
                result.materialBoundaryTarget.faces(face, 2));
            const double squaredDoubleArea =
                (b - a).cross(c - a).squaredNorm();
            if (squaredDoubleArea < minimumSquaredDoubleArea)
            {
                minimumSquaredDoubleArea = squaredDoubleArea;
                minimumAreaFace = face;
            }
        }
        std::ostringstream diagnostic;
        diagnostic <<
            "quadric material Boolean is not a positive closed manifold: "
            "volume=" << materialAudit.signedVolume <<
            " matrix_shape=" << materialAudit.matrixShapeValid <<
            " finite=" << materialAudit.finiteVertices <<
            " indices=" << materialAudit.faceIndicesValid <<
            " nondegenerate=" << materialAudit.nondegenerateFaces <<
            " degenerate_faces=" << materialAudit.degenerateFaceCount <<
            " min_squared_double_area=" << minimumSquaredDoubleArea <<
            " unreferenced_vertices=" <<
                materialAudit.unreferencedVertexCount <<
            " closed=" << materialAudit.closed <<
            " edge_manifold=" << materialAudit.edgeManifold <<
            " vertex_manifold=" << materialAudit.vertexManifold <<
            " oriented=" << materialAudit.consistentlyOriented <<
            " component_orientation=" <<
                materialAudit.componentOrientationsValid <<
            " intersection_free=" <<
                materialAudit.selfIntersections.intersectionFree <<
            " intersection_completed=" <<
                materialAudit.selfIntersections.completed <<
            " intersection_pairs=" <<
                materialAudit.selfIntersections.intersectingPairCount <<
            " first_faces=" << materialAudit.selfIntersections.firstFace <<
            ',' << materialAudit.selfIntersections.secondFace;
        if (minimumAreaFace >= 0)
        {
            const Eigen::Vector3i indices =
                result.materialBoundaryTarget.faces.row(minimumAreaFace);
            const Eigen::Vector3d a = result.materialBoundaryTarget.vertices.row(indices[0]);
            const Eigen::Vector3d b = result.materialBoundaryTarget.vertices.row(indices[1]);
            const Eigen::Vector3d c = result.materialBoundaryTarget.vertices.row(indices[2]);
            diagnostic << std::setprecision(17)
                << " minimum_area_face=" << minimumAreaFace
                << " area_witness_vertices=" << indices[0] << ',' << indices[1]
                << ',' << indices[2]
                << " area_witness_birth_face=" << birthFaces[minimumAreaFace]
                << " area_witness_points="
                << a.x() << ',' << a.y() << ',' << a.z() << ';'
                << b.x() << ',' << b.y() << ',' << b.z() << ';'
                << c.x() << ',' << c.y() << ',' << c.z();
        }
        throw std::runtime_error(diagnostic.str());
    }

    }
    stage("material-authority-audit-end");
    stage("moments-begin");
    MassMoments variable;
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(settings.triangulation.rationalSourceCoordinates) {
        result.variableMomentsExactTarget=rationalSurfaceMoments(result.materialExactVerticesTarget,
            result.materialBoundaryTarget.faces,request.materialDensity);
        variable=assignedRationalMoments(result.variableMomentsExactTarget);
    } else
#endif
    variable = computeSurfaceMassMoments(
        result.materialBoundaryTarget.vertices,
        result.materialBoundaryTarget.faces,
        request.materialDensity
    );
    if (!variable.isFinite() || variable.mass <= 0.0)
    {
        throw std::runtime_error(
            "quadric material has invalid integrated mass moments"
        );
    }
    result.variableMomentsTarget =
        physics::RawMassMoments10::fromMassMoments(variable);
    const Eigen::Matrix3d modelToTarget = request.targetFrame.axes.transpose();
    const Eigen::Vector3d targetTranslation =
        -modelToTarget * request.targetFrame.origin;
    const physics::RawMassMoments10 fixedTarget =
        request.fixedShellMoments.transformed(
            modelToTarget,
            targetTranslation
        );
    result.totalMomentsTarget = result.variableMomentsTarget + fixedTarget;
    stage("moments-end");

    const auto residualEnd = std::chrono::steady_clock::now();
    if(authority) {
        // Never apply represented-coefficient source velocities to the new
        // authority mapping. The actual value/material remains available for
        // its independent checks; Newton must reject this absent derivative.
        result.derivativeDefined=false;
        result.derivativeDiagnostic="coefficient authority source motion is not yet certified; Jacobian unavailable";
        if(authorityDerivative) {
            authorityDerivative->coefficientPolicy=authority->policy;
            authorityDerivative->precisionBits=authority->precisionBits;
            authorityDerivative->coefficientAuthority=authority->geometry();
            authorityDerivative->sourceConstructionIdentity=result.quadricTriangulation.rationalConstructionIdentity;
            auto sourceSettings=triangulationSettings;
            sourceSettings.truncationRadius=result.quadricTriangulation.truncationRadius;
            sourceSettings.maximumRadiusDoublings=0;
            stage("authority-source-motion-begin");
            const auto motion=evaluateQuadricSourceMotion(*authority,designTarget.vertices,
                sourceSettings,result.quadricTriangulation);
            stage("authority-source-motion-end");
            if(!motion.defined) {
                authorityDerivative->differential.diagnostic=motion.diagnostic;
                result.derivativeDiagnostic=motion.diagnostic;
            } else {
                authorityDerivative->sourceMotionScientificHash=motion.scientificHash;
                authorityDerivative->sourceDerivativeLower=motion.lower;
                authorityDerivative->sourceDerivativeUpper=motion.upper;
                stage("authority-material-transport-begin");
                authorityDerivative->differential=differentiateRationalPolyhedralMaterialMoments(
                    designTarget,result.quadricTriangulation.closedBoundary,motion.coefficientJacobian,
                    result.materialBoundaryTarget,birthFaces,request.materialDensity,
                    result.materialExactVerticesTarget,result.quadricTriangulation.exactVertices,
                    result.designExactVerticesTarget,settings.diagnosticStageObserver,
                    rationalMaterialAuthority?&*rationalMaterialAuthority:nullptr);
                stage("authority-material-transport-end");
                const auto& summary=authorityDerivative->differential.summary;
                result.derivativeDefined=summary.defined;
                result.momentJacobian=summary.jacobian;
                result.derivativeDiagnostic=authorityDerivative->differential.diagnostic;
                result.vertexProvenance=summary.vertexProvenance;
                result.constructedIntersectionVertexCount=summary.constructedIntersectionVertexCount;
                result.maximumConstraintResidual=summary.maximumConstraintResidual;
                result.minimumConstraintSingularRatio=summary.minimumConstraintSingularRatio;
                authorityDerivative->defined=summary.defined;
                if(summary.defined) {
                    std::vector<std::string> keys=result.vertexProvenance,faces;
                    std::sort(keys.begin(),keys.end());
                    for(int f=0;f<result.materialBoundaryTarget.faces.rows();++f) {
                        std::array<std::string,3> corners;
                        for(int k=0;k<3;++k)corners[k]=result.vertexProvenance[result.materialBoundaryTarget.faces(f,k)];
                        const int first=static_cast<int>(std::min_element(corners.begin(),corners.end())-corners.begin());
                        faces.push_back(corners[first]+";"+corners[(first+1)%3]+";"+
                            corners[(first+2)%3]+";"+std::to_string(birthFaces(f)));
                    }
                    std::sort(faces.begin(),faces.end());
                    std::ostringstream topology;topology<<"AuthoritySourceFeatureTopology-v1;";
                    for(const auto& key:keys)topology<<key<<'\n';
                    for(const auto& face:faces)topology<<face<<'\n';
                    result.provenanceTopologyHash=manufacturing::computeContentSha256(topology.str());
                }
                for(int f=0;f<birthFaces.size();++f)
                    if(birthFaces(f)>=designTarget.faces.rows())++result.interfaceTriangleCount;
            }
        }
    }
    else if (settings.derivativeMethod == QuadricMomentDerivativeMethod::LegacyContinuousImplicitShape)
    {
        accumulateInterfaceDerivative(
            result.materialBoundaryTarget,
            birthFaces,
            static_cast<int>(designTarget.faces.rows()),
            result.quadricTriangulation.interfaceFaceCount,
            quadricTarget,
            request.materialDensity,
            settings.minimumInterfaceGradient,
            &result.momentJacobian,
            &result.interfaceTriangleCount,
            &result.derivativeDefined
        );
        result.derivativeDiagnostic = result.derivativeDefined
            ? "legacy continuous implicit-shape quadrature; not the discrete geometry derivative"
            : "legacy interface derivative undefined";
    }
    else if (settings.derivativeMethod == QuadricMomentDerivativeMethod::DiscreteFullRankQuadric)
    {
        QuadricDiscreteMomentDifferential differential;
        if(!result.materialExactVerticesTarget.empty()) {
            auto sourceSettings=triangulationSettings;
            sourceSettings.truncationRadius=result.quadricTriangulation.truncationRadius;
            sourceSettings.maximumRadiusDoublings=0;
            stage("source-motion-begin");
            const auto motion=settings.triangulation.rationalSourceCoordinates
                ? evaluateQuadricSourceMotion(quadricTarget,designTarget.vertices,sourceSettings,result.quadricTriangulation)
                : evaluateQuadricSourceMotion(quadricTarget,designTarget.vertices,sourceSettings);
            stage("source-motion-end");
            stage("material-transport-begin");
            if(!motion.defined) differential.diagnostic=motion.diagnostic;
            else if(motion.triangulation.hash!=result.quadricTriangulation.hash)
                differential.diagnostic="invalid exact motion source differs from actual Boolean operand";
            else if(settings.triangulation.rationalSourceCoordinates)
                differential=differentiatePolyhedralMaterialMoments(
                    designTarget,result.quadricTriangulation.closedBoundary,motion.velocity,
                    result.materialBoundaryTarget,birthFaces,request.materialDensity,
                    settings.exactBinary64SelfIntersection,settings.exactBinary64SourceBoundary,
                    settings.exactBinary64MaterialBoundary,result.materialExactVerticesTarget,
                    result.quadricTriangulation.exactVertices,result.designExactVerticesTarget,
                    settings.diagnosticStageObserver,
                    rationalMaterialAuthority ? &*rationalMaterialAuthority : nullptr);
            else differential=differentiatePolyhedralMaterialMoments(
                designTarget,result.quadricTriangulation.closedBoundary,motion.velocity,
                result.materialBoundaryTarget,birthFaces,request.materialDensity,
                settings.exactBinary64SelfIntersection,settings.exactBinary64SourceBoundary,
                settings.exactBinary64MaterialBoundary,result.materialExactVerticesTarget,
                {},{},settings.diagnosticStageObserver);
            stage("material-transport-end");
        } else differential = differentiateDiscreteQuadricMoments(
            designTarget, result.quadricTriangulation, quadricTarget,
            result.materialBoundaryTarget, birthFaces, request.materialDensity,triangulationSettings);
        result.momentJacobian = differential.jacobian;
        result.derivativeDefined = differential.defined;
        result.derivativeDiagnostic = differential.diagnostic;
        result.vertexProvenance = differential.vertexProvenance;
        result.constructedIntersectionVertexCount = differential.constructedIntersectionVertexCount;
        result.maximumConstraintResidual = differential.maximumConstraintResidual;
        result.minimumConstraintSingularRatio = differential.minimumConstraintSingularRatio;
        for (int f = 0; f < birthFaces.size(); ++f)
            if (birthFaces(f) >= designTarget.faces.rows()) ++result.interfaceTriangleCount;
        if (differential.defined)
        {
            std::vector<std::string> vertices = result.vertexProvenance;
            std::sort(vertices.begin(), vertices.end());
            std::vector<std::string> faces;
            for (int f = 0; f < result.materialBoundaryTarget.faces.rows(); ++f)
            {
                std::array<std::string,3> keys;
                for (int k = 0; k < 3; ++k)
                    keys[k] = result.vertexProvenance[result.materialBoundaryTarget.faces(f,k)];
                const int first = static_cast<int>(std::min_element(keys.begin(),keys.end())-keys.begin());
                faces.push_back(keys[first]+";"+keys[(first+1)%3]+";"+keys[(first+2)%3]+";"+std::to_string(birthFaces(f)));
            }
            std::sort(faces.begin(),faces.end());
            std::ostringstream identity;
            identity << "QuadricSourceFeatureTopology-v1;";
            for (const auto& v : vertices) identity << v << '\n';
            for (const auto& f : faces) identity << f << '\n';
            result.provenanceTopologyHash = manufacturing::computeContentSha256(identity.str());
        }
    }
    else
    {
        result.stratumDifferential = differentiateDiscreteQuadricStratumMoments(
            designTarget, result.quadricTriangulation, quadricTarget,
            result.materialBoundaryTarget, birthFaces, request.materialDensity,
            triangulationSettings);
        result.derivativeDefined = false;
        result.derivativeDiagnostic = result.stratumDifferential.diagnostic;
        result.vertexProvenance = result.stratumDifferential.vertexProvenance;
        result.constructedIntersectionVertexCount =
            result.stratumDifferential.constructedIntersectionVertexCount;
        result.maximumConstraintResidual =
            result.stratumDifferential.maximumConstraintResidual;
        result.minimumConstraintSingularRatio =
            result.stratumDifferential.minimumConstraintSingularRatio;
        if (result.stratumDifferential.defined)
        {
            std::vector<std::string> vertices = result.vertexProvenance;
            std::sort(vertices.begin(), vertices.end());
            std::vector<std::string> faces;
            for (int f = 0; f < result.materialBoundaryTarget.faces.rows(); ++f)
            {
                std::array<std::string, 3> keys;
                for (int k = 0; k < 3; ++k)
                    keys[k] = result.vertexProvenance[
                        result.materialBoundaryTarget.faces(f, k)];
                const int first = static_cast<int>(
                    std::min_element(keys.begin(), keys.end()) - keys.begin());
                faces.push_back(keys[first] + ";" + keys[(first + 1) % 3] + ";" +
                    keys[(first + 2) % 3] + ";" + std::to_string(birthFaces(f)));
            }
            std::sort(faces.begin(), faces.end());
            std::ostringstream identity;
            identity << "QuadricSourceFeatureTopology-v1;";
            for (const auto& vertex : vertices) identity << vertex << '\n';
            for (const auto& face : faces) identity << face << '\n';
            result.provenanceTopologyHash =
                manufacturing::computeContentSha256(identity.str());
        }
        for (int f = 0; f < birthFaces.size(); ++f)
            if (birthFaces(f) >= designTarget.faces.rows())
                ++result.interfaceTriangleCount;
    }
    const auto derivativeEnd = std::chrono::steady_clock::now();
    result.residualEvaluationSeconds =
        std::chrono::duration<double>(residualEnd - evaluationStart).count();
    result.derivativeEvaluationSeconds =
        std::chrono::duration<double>(derivativeEnd - residualEnd).count();
    result.materialBoundary =
        toModel(result.materialBoundaryTarget, request.targetFrame);
    result.geometryHash = meshHash(
        result.materialBoundaryTarget,
        "QuadricMomentGeometry-v1"
    );
    if(settings.triangulation.rationalSourceCoordinates) {
        std::ostringstream identity;
        identity << "RationalMaterialGeometry-v1;";
        for(const auto& point:result.materialExactVerticesTarget)
            for(const auto& coordinate:point)identity << coordinate << ';';
        for(int f=0;f<result.materialBoundaryTarget.faces.rows();++f)
            for(int k=0;k<3;++k)identity << result.materialBoundaryTarget.faces(f,k) << ';';
        for(const auto& moment:result.variableMomentsExactTarget)identity << moment << ';';
        result.geometryHash=manufacturing::computeContentSha256(identity.str());
    }
    result.topologyHash = topologyHash(
        result.materialBoundaryTarget,
        birthFaces
    );
    std::ostringstream scientific;
    scientific.imbue(std::locale::classic());
    scientific << "QuadricMomentDifferentialScientific-v3;"
        << quadricMomentDerivativeMethodToken(result.derivativeMethod) << ';'
        << result.derivativeDefined << ';' << result.derivativeDiagnostic << ';'
        << result.booleanBackendIdentity << ';'
        << computeQuadricEvaluatorRequestHash(request, settings) << ';'
        << result.quadricTriangulation.hash << ';' << result.geometryHash << ';'
        << result.provenanceTopologyHash << ';' << std::hexfloat
        << settings.minimumInterfaceGradient << ';';
    for (int j = 0; j < 10; ++j) scientific << quadricTarget.values()(j) << ';';
    if(authority) {
        scientific<<authority->policy<<';'<<authority->precisionBits<<';'
                  <<result.quadricTriangulation.rationalConstructionIdentity<<';';
        for(const auto& coefficient:authority->geometry())scientific<<coefficient<<';';
    }
    if (result.stratumDifferential.defined)
    {
        scientific << "reduced-basis="
                   << result.stratumDifferential.coefficientBasis.rows() << 'x'
                   << result.stratumDifferential.coefficientBasis.cols() << ';';
        for (Eigen::Index row = 0;
             row < result.stratumDifferential.coefficientBasis.rows(); ++row)
            for (Eigen::Index column = 0;
                 column < result.stratumDifferential.coefficientBasis.cols(); ++column)
                scientific << result.stratumDifferential.coefficientBasis(row, column) << ';';
        for (Eigen::Index row = 0;
             row < result.stratumDifferential.jacobian.rows(); ++row)
            for (Eigen::Index column = 0;
                 column < result.stratumDifferential.jacobian.cols(); ++column)
                scientific << result.stratumDifferential.jacobian(row, column) << ';';
    }
    result.derivativeScientificHash = manufacturing::computeContentSha256(scientific.str());
    if(authorityDerivative) {
        std::ostringstream identity;identity<<"AuthorityActualBooleanMomentDifferential-v1;"
            <<result.derivativeScientificHash<<';'<<authorityDerivative->sourceMotionScientificHash<<';'
            <<authorityDerivative->defined<<';'
            <<authorityDerivative->differential.maximumConstraintResidualExact<<';';
        for(const auto& row:authorityDerivative->differential.coefficientJacobian)
            for(const auto& entry:row)identity<<entry<<';';
        for(const auto& matrix:authorityDerivative->differential.materialVertexJacobian)
            for(const auto& row:matrix)for(const auto& entry:row)identity<<entry<<';';
        authorityDerivative->scientificHash=manufacturing::computeContentSha256(identity.str());
    }
    return result;
}

QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,const QuadricCoefficients10& q,
    const QuadricMomentEvaluatorSettings& settings) {
    return evaluateQuadricMaterialMomentsImpl(request,q,settings,nullptr);
}

QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,const QuadricCoefficientAuthority10& authority,
    const QuadricMomentEvaluatorSettings& settings) {
    if(!settings.triangulation.rationalSourceCoordinates ||
       settings.derivativeMethod!=QuadricMomentDerivativeMethod::DiscreteFullRankQuadric ||
       settings.emptyRegionPolicy!=QuadricEmptyRegionPolicy::RejectEmpty)
        throw std::invalid_argument("coefficient authority value path requires explicit rational full-rank nonempty profile");
    QuadricCoefficients10::Vector display;
    const auto values=authority.geometryDisplay();
    for(int i=0;i<10;++i)display(i)=values[i];
    return evaluateQuadricMaterialMomentsImpl(request,QuadricCoefficients10::fromNormalized(display),settings,&authority);
}

QuadricMomentEvaluation evaluateQuadricMaterialMoments(
    const SpinDesignRequest& request,const QuadricCoefficientAuthority10& authority,
    const QuadricMomentEvaluatorSettings& settings,QuadricAuthorityMomentAttachment& attachment) {
    attachment=QuadricAuthorityMomentAttachment{};
    if(!settings.triangulation.rationalSourceCoordinates ||
       settings.derivativeMethod!=QuadricMomentDerivativeMethod::DiscreteFullRankQuadric ||
       settings.emptyRegionPolicy!=QuadricEmptyRegionPolicy::RejectEmpty)
        throw std::invalid_argument("coefficient authority derivative requires explicit rational full-rank nonempty profile");
    QuadricCoefficients10::Vector display;
    const auto values=authority.geometryDisplay();
    for(int i=0;i<10;++i)display(i)=values[i];
    return evaluateQuadricMaterialMomentsImpl(request,QuadricCoefficients10::fromNormalized(display),settings,&authority,&attachment);
}

QuadricAuthorityMaterialContext makeQuadricAuthorityMaterialContext(
    const SpinDesignRequest& input,double density) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(!std::isfinite(density) || density<=0)
        throw std::invalid_argument("fixed material density must be positive and finite");
    QuadricAuthorityMaterialContext out;
    out.request_=input;out.fixedDensity_=density;
    out.request_.fixedShellMoments=physics::RawMassMoments10{};
    validateSpinDesignRequest(out.request_);
    const SpinTargetFrame modelFrame;
    const auto outerModel=rationalTargetVertices(input.outerSurface,modelFrame);
    const auto designModel=rationalTargetVertices(input.designDomain,modelFrame);
    const auto outerTarget=rationalTargetVertices(input.outerSurface,input.targetFrame);
    const auto designTarget=rationalTargetVertices(input.designDomain,input.targetFrame);
    out.outerRegion_=std::make_shared<verification::ExactMaterialRegion>(
        verification::makeExactMaterialRegionFromRationals(outerTarget,input.outerSurface.faces));
    out.designRegion_=std::make_shared<verification::ExactMaterialRegion>(
        verification::makeExactMaterialRegionFromRationals(designTarget,input.designDomain.faces));
    for(const auto* region:{out.outerRegion_.get(),out.designRegion_.get()})
        if(!region->facts().completed || !region->facts().valid || region->facts().explicitEmpty ||
           !region->facts().checkerReversedComponents.empty())
            throw std::invalid_argument("fixed material requires valid original oriented geometry");
    const auto subset=verification::checkExactMaterialSubset(*out.designRegion_,*out.outerRegion_);
    if(!subset.completed || !subset.aSubsetB || !*subset.aSubsetB)
        throw std::invalid_argument("design domain containment not established");
    const auto om=rationalSurfaceMoments(outerModel,input.outerSurface.faces,density);
    const auto dm=rationalSurfaceMoments(designModel,input.designDomain.faces,density);
    const auto ot=rationalSurfaceMoments(outerTarget,input.outerSurface.faces,density);
    const auto dt=rationalSurfaceMoments(designTarget,input.designDomain.faces,density);
    std::vector<std::string> modelMoments(10);
    for(int i=0;i<10;++i) {
        modelMoments[i]=rationalText(CGAL::Gmpq(om[i])-CGAL::Gmpq(dm[i]));
        out.fixedMomentsTarget_[i]=rationalText(CGAL::Gmpq(ot[i])-CGAL::Gmpq(dt[i]));
    }
    if(CGAL::Gmpq(modelMoments[0])<0)throw std::invalid_argument("negative fixed material mass");
    out.request_.fixedShellMoments=physics::RawMassMoments10::fromMassMoments(assignedRationalMoments(modelMoments));
    validateSpinDesignRequest(out.request_);
    CGAL::Gmpq a[3][3],inverse[3][3],b[3];
    for(int i=0;i<3;++i) {
        b[i]=0;
        for(int j=0;j<3;++j) {
            a[i][j]=CGAL::Gmpq(input.targetFrame.axes(j,i));
            b[i]-=a[i][j]*CGAL::Gmpq(input.targetFrame.origin(j));
            out.modelToTarget_[i][j]=rationalText(a[i][j]);
        }
        out.modelToTarget_[i][3]=rationalText(b[i]);
    }
    const auto det=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-
        a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    if(det==0)throw std::invalid_argument("singular exact target frame");
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
        const int r0=(j+1)%3,r1=(j+2)%3,c0=(i+1)%3,c1=(i+2)%3;
        inverse[i][j]=(a[r0][c0]*a[r1][c1]-a[r0][c1]*a[r1][c0])/det;
        out.targetToModel_[i][j]=rationalText(inverse[i][j]);
    }
    for(int i=0;i<3;++i) {
        CGAL::Gmpq translation=0;
        for(int j=0;j<3;++j)translation-=inverse[i][j]*b[j];
        out.targetToModel_[i][3]=rationalText(translation);
    }
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"GeometryOwnedFixedMaterial-v1;"<<computeSpinDesignRequestHash(out.request_)<<';'
        <<rationalText(CGAL::Gmpq(density))<<';'<<out.outerRegion_->facts().regionDigest<<';'
        <<out.designRegion_->facts().regionDigest<<';';
    for(const auto& v:out.fixedMomentsTarget_)identity<<v<<';';
    for(const auto* affine:{&out.modelToTarget_,&out.targetToModel_})
        for(const auto& row:*affine)for(const auto& v:row)identity<<v<<';';
    out.identity_=manufacturing::computeContentSha256(identity.str());
    return out;
#else
    throw std::runtime_error("geometry-owned authority material requires exact geometry backend");
#endif
}

QuadricAuthorityTotalMomentEvaluation evaluateQuadricMaterialMoments(
    const QuadricAuthorityMaterialContext& context,const QuadricCoefficientAuthority10& authority,
    const QuadricMomentEvaluatorSettings& settings) {
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    QuadricAuthorityTotalMomentEvaluation out;
    out.material=evaluateQuadricMaterialMoments(context.request(),authority,settings,out.derivative);
    if(!out.derivative.defined || out.material.variableMomentsExactTarget.size()!=10)
        throw std::runtime_error("actual authority material differential unavailable");
    out.contextIdentity=context.identity();
    std::ostringstream identity;identity.imbue(std::locale::classic());
    identity<<"GeometryOwnedTotalMomentEvaluation-v1;"<<out.contextIdentity<<';'
        <<out.derivative.scientificHash<<';';
    for(int i=0;i<10;++i) {
        out.totalMomentsTarget[i]=rationalText(CGAL::Gmpq(context.fixedMomentsTarget()[i])+
            CGAL::Gmpq(out.material.variableMomentsExactTarget[i]));
        identity<<out.totalMomentsTarget[i]<<';';
    }
    out.scientificHash=manufacturing::computeContentSha256(identity.str());
    return out;
#else
    throw std::runtime_error("geometry-owned authority evaluation requires exact geometry backend");
#endif
}

} // namespace spin
