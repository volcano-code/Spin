#include "spin/PlanarDiskIntersection.h"
#include "physics/MassProperties.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <list>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace spin
{
namespace
{

using Wide = long double;
using Moments = std::array<Wide, 6>;
constexpr Wide kPi = 3.141592653589793238462643383279502884L;
constexpr Wide kTau = 2.0L * kPi;

Wide finiteDouble(Wide value, const char* label)
{
    const double converted = static_cast<double>(value);
    if (!std::isfinite(converted))
        throw std::overflow_error(label);
    return static_cast<Wide>(converted);
}

Moments integrateLine(Wide x0, Wide y0, Wide x1, Wide y1)
{
    const Wide dx = x1 - x0;
    const Wide dy = y1 - y0;
    const Wide x2 = x0*x0;
    const Wide y2 = y0*y0;
    const Wide dx2 = dx*dx;
    const Wide dy2 = dy*dy;

    return {
        dy * (x0 + dx/2.0L),
        dy/2.0L * (x2 + x0*dx + dx2/3.0L),
        dy * (x0*y0 + (x0*dy + y0*dx)/2.0L + dx*dy/3.0L),
        dy/3.0L * (x2*x0 + 1.5L*x2*dx + x0*dx2 + dx2*dx/4.0L),
        dy * (x0*y2 + (2.0L*x0*y0*dy + dx*y2)/2.0L +
              (x0*dy2 + 2.0L*dx*y0*dy)/3.0L + dx*dy2/4.0L),
        dy/2.0L * (x2*y0 + (x2*dy + 2.0L*x0*dx*y0)/2.0L +
                   (2.0L*x0*dx*dy + dx2*y0)/3.0L + dx2*dy/4.0L)
    };
}

Moments integrateArc(Wide radius, Wide begin, Wide end)
{
    const auto fArea = [](Wide t) { return t/2.0L + std::sin(2.0L*t)/4.0L; };
    const auto fFirstX = [](Wide t) {
        const Wide s = std::sin(t);
        return s - s*s*s/3.0L;
    };
    const auto fFirstY = [](Wide t) { return -std::pow(std::cos(t), 3)/3.0L; };
    const auto fSecondXX = [](Wide t) {
        return 3.0L*t/8.0L + std::sin(2.0L*t)/4.0L + std::sin(4.0L*t)/32.0L;
    };
    const auto fSecondYY = [](Wide t) { return t/8.0L - std::sin(4.0L*t)/32.0L; };
    const auto fSecondXY = [](Wide t) { return -std::pow(std::cos(t), 4)/4.0L; };
    const Wide r2 = radius*radius;
    const Wide r3 = r2*radius;
    const Wide r4 = r2*r2;
    return {
        r2*(fArea(end) - fArea(begin)),
        r3/2.0L*(fFirstX(end) - fFirstX(begin)),
        r3*(fFirstY(end) - fFirstY(begin)),
        r4/3.0L*(fSecondXX(end) - fSecondXX(begin)),
        r4*(fSecondYY(end) - fSecondYY(begin)),
        r4/2.0L*(fSecondXY(end) - fSecondXY(begin))
    };
}

bool insideRegion(Wide x, Wide y,
                  const std::vector<PlanarBoundaryLoop>& loops,
                  Wide centerX, Wide centerY)
{
    int winding = 0;
    for (const auto& loop : loops)
    {
        std::size_t count = loop.size();
        if (count > 1 && loop.front().x() == loop.back().x() &&
                         loop.front().y() == loop.back().y())
            --count;
        for (std::size_t i = 0; i < count; ++i)
        {
            const Wide ax = static_cast<Wide>(loop[i].x()) - centerX;
            const Wide ay = static_cast<Wide>(loop[i].y()) - centerY;
            const Wide bx = static_cast<Wide>(loop[(i+1)%count].x()) - centerX;
            const Wide by = static_cast<Wide>(loop[(i+1)%count].y()) - centerY;
            const Wide cross = (bx-ax)*(y-ay) - (x-ax)*(by-ay);
            if (ay <= y)
            {
                if (by > y && cross > 0.0L) ++winding;
            }
            else if (by <= y && cross < 0.0L)
            {
                --winding;
            }
        }
    }
    return winding != 0;
}

void add(Moments& total, const Moments& part)
{
    for (std::size_t i = 0; i < total.size(); ++i) total[i] += part[i];
}

void addCircleIntersections(Wide x0, Wide y0, Wide dx, Wide dy, Wide radius,
                            std::vector<Wide>& angles)
{
    const Wide a = dx*dx + dy*dy;
    const Wide b = 2.0L*(x0*dx + y0*dy);
    const Wide c = x0*x0 + y0*y0 - radius*radius;
    const Wide discriminant = b*b - 4.0L*a*c;
    if (discriminant < 0.0L) return;
    const Wide root = std::sqrt(std::max(0.0L, discriminant));
    const std::array<Wide,2> roots = {
        (-b-root)/(2.0L*a), (-b+root)/(2.0L*a)
    };
    for (Wide t : roots)
    {
        if (t < 0.0L || t > 1.0L) continue;
        Wide angle = std::atan2(y0+t*dy, x0+t*dx);
        if (angle < 0.0L) angle += kTau;
        angles.push_back(angle);
    }
}

void integrateClippedEdge(Wide x0, Wide y0, Wide x1, Wide y1, Wide radius,
                          Moments& total)
{
    const Wide dx = x1-x0;
    const Wide dy = y1-y0;
    const Wide a = dx*dx + dy*dy;
    const Wide b = 2.0L*(x0*dx + y0*dy);
    const Wide c = x0*x0 + y0*y0 - radius*radius;
    const Wide discriminant = b*b - 4.0L*a*c;
    Wide begin = 0.0L, end = 1.0L;
    if (discriminant < 0.0L)
    {
        if (c > 0.0L) return;
    }
    else
    {
        const Wide root = std::sqrt(std::max(0.0L, discriminant));
        const Wide t0 = (-b-root)/(2.0L*a);
        const Wide t1 = (-b+root)/(2.0L*a);
        begin = std::max(begin, t0);
        end = std::min(end, t1);
        if (!(end > begin)) return;
    }
    add(total, integrateLine(x0+begin*dx, y0+begin*dy,
                             x0+end*dx, y0+end*dy));
}

PlanarDiskIntersectionTopology classifyTopology(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk)
{
    if (!disk.center.allFinite() || !std::isfinite(disk.radius) || disk.radius < 0.0)
        throw std::invalid_argument("disk center must be finite and radius nonnegative");
    (void)integratePlanarBoundaryLoops(loops);

    PlanarDiskIntersectionTopology result;
    if (disk.radius == 0.0)
    {
        result.signature = "ZERO_RADIUS";
        return result;
    }

    const Wide cx = static_cast<Wide>(disk.center.x());
    const Wide cy = static_cast<Wide>(disk.center.y());
    const Wide radius = static_cast<Wide>(disk.radius);
    const Wide machine = std::numeric_limits<Wide>::epsilon();
    std::ostringstream signature;
    for (std::size_t loopIndex = 0; loopIndex < loops.size(); ++loopIndex)
    {
        const auto& loop = loops[loopIndex];
        std::size_t count = loop.size();
        if (count > 1 && loop.front().x() == loop.back().x() &&
                         loop.front().y() == loop.back().y())
            --count;
        std::size_t strictCrossings = 0;
        for (std::size_t i = 0; i < count; ++i)
        {
            const auto& p = loop[i];
            const auto& q = loop[(i+1)%count];
            const Wide x0 = static_cast<Wide>(p.x())-cx;
            const Wide y0 = static_cast<Wide>(p.y())-cy;
            const Wide dx = static_cast<Wide>(q.x())-static_cast<Wide>(p.x());
            const Wide dy = static_cast<Wide>(q.y())-static_cast<Wide>(p.y());
            const Wide a = dx*dx + dy*dy;
            const Wide b = 2.0L*(x0*dx + y0*dy);
            const Wide c = x0*x0 + y0*y0 - radius*radius;
            const Wide x1 = x0 + dx;
            const Wide y1 = y0 + dy;
            const Wide endValue = x1*x1 + y1*y1 - radius*radius;
            const Wide valueScale = std::max({1.0L, std::abs(x0*x0+y0*y0),
                                               std::abs(x1*x1+y1*y1), radius*radius});
            const Wide valueTolerance = 128.0L*machine*valueScale;
            if (std::abs(c) <= valueTolerance || std::abs(endValue) <= valueTolerance ||
                !(a > 0.0L))
            {
                result.nearTransition = true;
                continue;
            }
            const Wide discriminant = b*b - 4.0L*a*c;
            const Wide discriminantScale = std::max({1.0L, std::abs(b*b),
                                                       std::abs(4.0L*a*c)});
            const Wide discriminantTolerance = 128.0L*machine*discriminantScale;
            if (std::abs(discriminant) <= discriminantTolerance)
            {
                result.nearTransition = true;
                continue;
            }
            if (discriminant < 0.0L) continue;

            const Wide root = std::sqrt(discriminant);
            const std::array<Wide, 2> roots = {
                (-b-root)/(2.0L*a), (-b+root)/(2.0L*a)
            };
            const Wide parameterTolerance = 128.0L*machine;
            for (const Wide t : roots)
            {
                if (std::abs(t) <= parameterTolerance ||
                    std::abs(t-1.0L) <= parameterTolerance)
                    result.nearTransition = true;
                else if (t > 0.0L && t < 1.0L)
                    ++strictCrossings;
            }
        }
        if (loopIndex != 0) signature << ';';
        signature << strictCrossings;
    }
    result.signature = signature.str();
    return result;
}

} // namespace

PlanarDiskIntersectionTopology classifyPlanarDiskIntersectionTopology(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk)
{
    return classifyTopology(loops, disk);
}

PlanarRawMoments6 integratePlanarDiskIntersection(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk)
{
    if (!disk.center.allFinite() || !std::isfinite(disk.radius) || disk.radius <= 0.0)
        throw std::invalid_argument("disk center must be finite and radius positive");
    (void)integratePlanarBoundaryLoops(loops); // Validate loop structure and net orientation.

    const Wide cx = static_cast<Wide>(disk.center.x());
    const Wide cy = static_cast<Wide>(disk.center.y());
    const Wide radius = static_cast<Wide>(disk.radius);
    Moments local{};
    std::vector<Wide> angles;

    for (const auto& loop : loops)
    {
        std::size_t count = loop.size();
        if (count > 1 && loop.front().x() == loop.back().x() &&
                         loop.front().y() == loop.back().y())
            --count;
        for (std::size_t i = 0; i < count; ++i)
        {
            const auto& p = loop[i];
            const auto& q = loop[(i+1)%count];
            const Wide x0 = static_cast<Wide>(p.x())-cx;
            const Wide y0 = static_cast<Wide>(p.y())-cy;
            const Wide x1 = static_cast<Wide>(q.x())-cx;
            const Wide y1 = static_cast<Wide>(q.y())-cy;
            addCircleIntersections(x0, y0, x1-x0, y1-y0, radius, angles);
            integrateClippedEdge(x0, y0, x1, y1, radius, local);
        }
    }

    std::sort(angles.begin(), angles.end());
    const Wide dedupTolerance = 64.0L*std::numeric_limits<Wide>::epsilon()*kTau;
    angles.erase(std::unique(angles.begin(), angles.end(), [=](Wide a, Wide b) {
        return std::abs(a-b) <= dedupTolerance;
    }), angles.end());

    if (angles.empty())
    {
        if (insideRegion(radius, 0.0L, loops, cx, cy))
            add(local, integrateArc(radius, 0.0L, kTau));
    }
    else
    {
        for (std::size_t i = 0; i < angles.size(); ++i)
        {
            const Wide begin = angles[i];
            const Wide end = (i+1 < angles.size()) ? angles[i+1] : angles[0]+kTau;
            const Wide middle = (begin+end)/2.0L;
            if (insideRegion(radius*std::cos(middle), radius*std::sin(middle),
                             loops, cx, cy))
                add(local, integrateArc(radius, begin, end));
        }
    }

    if (local[0] < 0.0L)
        throw std::runtime_error("polygon/disk intersection produced negative oriented area");
    if (local[0] == 0.0L)
        return {};

    const Wide area = local[0];
    const Wide firstX = local[1] + cx*area;
    const Wide firstY = local[2] + cy*area;
    const Wide secondXX = local[3] + 2.0L*cx*local[1] + cx*cx*area;
    const Wide secondYY = local[4] + 2.0L*cy*local[2] + cy*cy*area;
    const Wide secondXY = local[5] + cx*local[2] + cy*local[1] + cx*cy*area;
    const double convertedArea = static_cast<double>(finiteDouble(area, "intersection area overflowed"));
    if (!(convertedArea > 0.0))
        throw std::overflow_error("positive intersection area underflowed double precision");
    return {convertedArea,
            static_cast<double>(finiteDouble(firstX, "intersection first moment overflowed")),
            static_cast<double>(finiteDouble(firstY, "intersection first moment overflowed")),
            static_cast<double>(finiteDouble(secondXX, "intersection second moment overflowed")),
            static_cast<double>(finiteDouble(secondYY, "intersection second moment overflowed")),
            static_cast<double>(finiteDouble(secondXY, "intersection mixed moment overflowed"))};
}

PlanarRegionDifferenceBoundary planarRegionDifferenceBoundary(
    const std::vector<PlanarBoundaryLoop>& loops,
    const PlanarCircleDisk& disk)
{
    if (!disk.center.allFinite() || !std::isfinite(disk.radius) || disk.radius <= 0.0)
        throw std::invalid_argument("disk center must be finite and radius positive");
    (void)integratePlanarBoundaryLoops(loops);

    PlanarRegionDifferenceBoundary result;
    const PlanarDiskIntersectionTopology topology =
        classifyTopology(loops, disk);
    if (topology.nearTransition)
    {
        result.diagnostic = "polygon/disk boundary is tangent or near a topology transition";
        return result;
    }

    double coordinateScale = std::max(1.0, disk.center.cwiseAbs().maxCoeff() + disk.radius);
    for (const auto& loop : loops)
        for (const auto& point : loop)
            coordinateScale = std::max(coordinateScale, point.cwiseAbs().maxCoeff());
    const double pointTolerance = 256.0 * std::numeric_limits<double>::epsilon() * coordinateScale;
    const Wide cx = static_cast<Wide>(disk.center.x());
    const Wide cy = static_cast<Wide>(disk.center.y());
    const Wide radius = static_cast<Wide>(disk.radius);
    std::vector<PlanarBoundaryPiece> pieces;
    std::vector<Wide> angles;

    const auto circlePoint = [&](Wide angle) {
        return Eigen::Vector2d(
            disk.center.x() + disk.radius * static_cast<double>(std::cos(angle)),
            disk.center.y() + disk.radius * static_cast<double>(std::sin(angle)));
    };
    const auto addRoots = [&](Wide x0, Wide y0, Wide dx, Wide dy,
                              std::vector<Wide>& roots) {
        const Wide a = dx*dx + dy*dy;
        if (!(a > 0.0L)) return;
        const Wide b = 2.0L*(x0*dx + y0*dy);
        const Wide c = x0*x0 + y0*y0 - radius*radius;
        const Wide discriminant = b*b - 4.0L*a*c;
        if (discriminant <= 0.0L) return;
        const Wide root = std::sqrt(discriminant);
        const Wide values[2] = {(-b-root)/(2.0L*a), (-b+root)/(2.0L*a)};
        for (const Wide t : values)
        {
            if (t > 0.0L && t < 1.0L)
            {
                roots.push_back(t);
                Wide angle = std::atan2(y0+t*dy, x0+t*dx);
                if (angle < 0.0L) angle += kTau;
                angles.push_back(angle);
            }
        }
    };

    for (const auto& loop : loops)
    {
        std::size_t count = loop.size();
        if (count > 1 && loop.front().x() == loop.back().x() &&
                         loop.front().y() == loop.back().y())
            --count;
        for (std::size_t i = 0; i < count; ++i)
        {
            const Eigen::Vector2d& p = loop[i];
            const Eigen::Vector2d& q = loop[(i+1)%count];
            const Wide x0 = static_cast<Wide>(p.x()) - cx;
            const Wide y0 = static_cast<Wide>(p.y()) - cy;
            const Wide dx = static_cast<Wide>(q.x()) - static_cast<Wide>(p.x());
            const Wide dy = static_cast<Wide>(q.y()) - static_cast<Wide>(p.y());
            std::vector<Wide> roots{0.0L, 1.0L};
            addRoots(x0, y0, dx, dy, roots);
            std::sort(roots.begin(), roots.end());
            const Wide tTolerance = 64.0L * std::numeric_limits<Wide>::epsilon();
            roots.erase(std::unique(roots.begin(), roots.end(), [=](Wide a, Wide b) {
                return std::abs(a-b) <= tTolerance;
            }), roots.end());
            for (std::size_t interval = 0; interval + 1 < roots.size(); ++interval)
            {
                const Wide begin = roots[interval];
                const Wide end = roots[interval+1];
                const Wide middle = (begin+end)/2.0L;
                const Wide mx = x0 + middle*dx;
                const Wide my = y0 + middle*dy;
                if (mx*mx + my*my <= radius*radius) continue;
                const Eigen::Vector2d start(
                    static_cast<double>(static_cast<Wide>(p.x()) + begin*dx),
                    static_cast<double>(static_cast<Wide>(p.y()) + begin*dy));
                const Eigen::Vector2d endPoint(
                    static_cast<double>(static_cast<Wide>(p.x()) + end*dx),
                    static_cast<double>(static_cast<Wide>(p.y()) + end*dy));
                PlanarBoundaryPiece piece;
                piece.type = PlanarBoundaryPieceType::LineSegment;
                piece.start = start;
                piece.end = endPoint;
                pieces.push_back(std::move(piece));
            }
        }
    }

    std::sort(angles.begin(), angles.end());
    const Wide angleTolerance = 128.0L * std::numeric_limits<double>::epsilon();
    angles.erase(std::unique(angles.begin(), angles.end(), [=](Wide a, Wide b) {
        return std::abs(a-b) <= angleTolerance;
    }), angles.end());
    const std::size_t intervalCount = angles.empty() ? 1 : angles.size();
    for (std::size_t i = 0; i < intervalCount; ++i)
    {
        const Wide begin = angles.empty() ? 0.0L : angles[i];
        const Wide end = angles.empty() ? kTau :
            (i+1 < angles.size() ? angles[i+1] : angles[0]+kTau);
        const Wide middle = (begin+end)/2.0L;
        if (!insideRegion(radius*std::cos(middle), radius*std::sin(middle),
                          loops, cx, cy))
            continue;
        PlanarBoundaryPiece piece;
        piece.type = PlanarBoundaryPieceType::CircularArc;
        piece.center = disk.center;
        piece.radius = disk.radius;
        // The polygon-minus-disk cavity is clockwise so the remaining
        // material stays on the left of the directed boundary.
        piece.startAngle = static_cast<double>(end);
        piece.sweep = -static_cast<double>(end-begin);
        piece.start = circlePoint(end);
        piece.end = circlePoint(begin);
        pieces.push_back(std::move(piece));
    }

    if (pieces.empty())
    {
        result.defined = true;
        return result;
    }

    std::vector<Eigen::Vector2d> nodes;
    std::vector<std::size_t> startNode(pieces.size());
    std::vector<std::size_t> endNode(pieces.size());
    auto findNode = [&](const Eigen::Vector2d& point) {
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            if ((nodes[i]-point).norm() <= pointTolerance) return i;
        }
        nodes.push_back(point);
        return nodes.size()-1;
    };
    std::vector<std::vector<std::size_t>> outgoing;
    std::vector<std::size_t> incoming;
    for (std::size_t i = 0; i < pieces.size(); ++i)
    {
        startNode[i] = findNode(pieces[i].start);
        endNode[i] = findNode(pieces[i].end);
    }
    outgoing.resize(nodes.size());
    incoming.resize(nodes.size(), 0);
    for (std::size_t i = 0; i < pieces.size(); ++i)
    {
        outgoing[startNode[i]].push_back(i);
        ++incoming[endNode[i]];
    }
    for (std::size_t node = 0; node < nodes.size(); ++node)
    {
        if (outgoing[node].size() != 1 || incoming[node] != 1)
        {
            result.diagnostic = "clipped boundary endpoint arrangement is not a set of closed cycles";
            return result;
        }
    }

    std::vector<bool> used(pieces.size(), false);
    for (std::size_t first = 0; first < pieces.size(); ++first)
    {
        if (used[first]) continue;
        PlanarBoundaryCycle cycle;
        std::size_t current = first;
        const std::size_t cycleStartNode = startNode[first];
        for (std::size_t guard = 0; guard <= pieces.size(); ++guard)
        {
            if (used[current])
            {
                result.diagnostic = "clipped boundary cycles intersect or repeat a piece";
                return result;
            }
            used[current] = true;
            cycle.push_back(pieces[current]);
            const std::size_t nextNode = endNode[current];
            if (nextNode == cycleStartNode) break;
            current = outgoing[nextNode].front();
            if (guard == pieces.size())
            {
                result.diagnostic = "clipped boundary cycle did not close";
                return result;
            }
        }
        result.cycles.push_back(std::move(cycle));
    }
    result.defined = true;
    return result;
}

PlanarRawMoments6 integratePlanarBoundaryCycles(
    const std::vector<PlanarBoundaryCycle>& cycles)
{
    Moments total{};
    for (const auto& cycle : cycles)
    {
        if (cycle.empty())
            throw std::invalid_argument("planar boundary cycle must not be empty");
        for (const PlanarBoundaryPiece& piece : cycle)
        {
            if (!piece.start.allFinite() || !piece.end.allFinite())
                throw std::invalid_argument("planar boundary piece coordinates must be finite");
            if (piece.type == PlanarBoundaryPieceType::LineSegment)
            {
                add(total, integrateLine(piece.start.x(), piece.start.y(),
                                         piece.end.x(), piece.end.y()));
                continue;
            }
            if (!piece.center.allFinite() || !std::isfinite(piece.radius) ||
                piece.radius <= 0.0 || !std::isfinite(piece.startAngle) ||
                !std::isfinite(piece.sweep) || piece.sweep == 0.0 ||
                std::abs(piece.sweep) > 2.0*static_cast<double>(kPi) + 1e-12)
                throw std::invalid_argument("circular boundary arc parameters are invalid");
            const Wide begin = static_cast<Wide>(piece.startAngle);
            const Wide end = begin + static_cast<Wide>(piece.sweep);
            Moments local = piece.sweep > 0.0
                ? integrateArc(piece.radius, begin, end)
                : integrateArc(piece.radius, end, begin);
            if (piece.sweep < 0.0)
                for (Wide& value : local) value = -value;
            const Wide cx = static_cast<Wide>(piece.center.x());
            const Wide cy = static_cast<Wide>(piece.center.y());
            // Green primitives are not translation invariant on an OPEN
            // arc: e.g. integral(x dy) gains cx*(y_end-y_begin). Close the
            // arc with its analytic chord, translate the resulting signed
            // region moments, then remove that chord in global coordinates.
            // Whole-circle arcs need no special case; their chord vanishes.
            const Wide x0 = piece.radius*std::cos(begin);
            const Wide y0 = piece.radius*std::sin(begin);
            const Wide x1 = piece.radius*std::cos(end);
            const Wide y1 = piece.radius*std::sin(end);
            add(local, integrateLine(x1, y1, x0, y0));
            const Wide area = local[0];
            local[1] += cx*area;
            local[2] += cy*area;
            local[3] += 2.0L*cx*local[1] - cx*cx*area;
            local[4] += 2.0L*cy*local[2] - cy*cy*area;
            local[5] += cx*local[2] + cy*local[1] - cx*cy*area;
            const Moments globalChord = integrateLine(x1+cx, y1+cy, x0+cx, y0+cy);
            for (std::size_t i = 0; i < local.size(); ++i)
                local[i] -= globalChord[i];
            add(total, local);
        }
    }
    if (total[0] < 0.0L)
        throw std::runtime_error("oriented planar boundary cycles have negative material area");
    return {
        static_cast<double>(finiteDouble(total[0], "boundary area overflowed")),
        static_cast<double>(finiteDouble(total[1], "boundary first moment overflowed")),
        static_cast<double>(finiteDouble(total[2], "boundary first moment overflowed")),
        static_cast<double>(finiteDouble(total[3], "boundary second moment overflowed")),
        static_cast<double>(finiteDouble(total[4], "boundary second moment overflowed")),
        static_cast<double>(finiteDouble(total[5], "boundary mixed moment overflowed"))};
}

bool planarRegionsHavePositiveAreaIntersection(
    const std::vector<PlanarBoundaryLoop>& first,
    const std::vector<PlanarBoundaryLoop>& second)
{
    using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
    using Polygon = CGAL::Polygon_2<Kernel>;
    using PolygonSet = CGAL::Polygon_set_2<Kernel>;
    using PolygonWithHoles = typename PolygonSet::Polygon_with_holes_2;
    using Exact = Kernel::FT;

    const auto makeRegion = [](const std::vector<PlanarBoundaryLoop>& loops) {
        using LocalKernel = CGAL::Exact_predicates_exact_constructions_kernel;
        using LocalPolygon = CGAL::Polygon_2<LocalKernel>;
        using LocalSet = CGAL::Polygon_set_2<LocalKernel>;
        using LocalExact = LocalKernel::FT;
        std::vector<std::pair<LocalExact, LocalPolygon>> ordered;
        ordered.reserve(loops.size());
        for (const auto& loop : loops)
        {
            if (loop.size() < 3)
                throw std::invalid_argument("planar overlap region loop needs at least three points");
            LocalPolygon polygon;
            for (const auto& point : loop)
            {
                if (!point.allFinite())
                    throw std::invalid_argument("planar overlap region coordinates must be finite");
                polygon.push_back(LocalKernel::Point_2(point.x(), point.y()));
            }
            if (!polygon.is_simple())
                throw std::invalid_argument("planar overlap region loop is not simple");
            const LocalExact area = polygon.area();
            if (area == LocalExact(0))
                throw std::invalid_argument("planar overlap region loop has zero area");
            ordered.emplace_back(area < LocalExact(0) ? -area : area, std::move(polygon));
        }
        std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
            return a.first > b.first;
        });
        LocalSet region;
        for (const auto& item : ordered)
        {
            if (item.second.orientation() == CGAL::COUNTERCLOCKWISE)
                region.join(item.second);
            else
                region.difference(item.second);
        }
        return region;
    };

    PolygonSet intersection = makeRegion(first);
    const PolygonSet other = makeRegion(second);
    intersection.intersection(other);
    std::list<PolygonWithHoles> components;
    intersection.polygons_with_holes(std::back_inserter(components));
    const auto absolute = [](const Exact& value) { return value < Exact(0) ? -value : value; };
    for (const auto& component : components)
    {
        if (component.is_unbounded())
            throw std::runtime_error("bounded planar overlap produced an unbounded component");
        Exact area = absolute(component.outer_boundary().area());
        for (auto hole = component.holes_begin(); hole != component.holes_end(); ++hole)
            area -= absolute(hole->area());
        if (area > Exact(0))
            return true;
    }
    return false;
}

physics::RawMassMoments10 integratePlanarBoundaryCyclesExtruded(
    const std::vector<PlanarBoundaryCycle>& cycles,
    const double z0,
    const double z1,
    const double volumeDensity)
{
    if (!std::isfinite(z0) || !std::isfinite(z1) || !(z1 > z0))
        throw std::invalid_argument("extrusion interval must be finite and increasing");
    if (!std::isfinite(volumeDensity) || !(volumeDensity > 0.0))
        throw std::invalid_argument("extrusion density must be finite and positive");

    const PlanarRawMoments6 planar = integratePlanarBoundaryCycles(cycles);
    const double height = z1 - z0;
    const double zFirstIntegral = (z1*z1 - z0*z0) / 2.0;
    const double zSecondIntegral = (z1*z1*z1 - z0*z0*z0) / 3.0;

    MassMoments moments;
    moments.mass = volumeDensity * height * planar.area;
    moments.firstMoment <<
        volumeDensity * height * planar.firstX,
        volumeDensity * height * planar.firstY,
        volumeDensity * zFirstIntegral * planar.area;
    moments.secondMoment <<
        volumeDensity * height * planar.secondXX,
        volumeDensity * height * planar.secondXY,
        volumeDensity * zFirstIntegral * planar.firstX,
        volumeDensity * height * planar.secondXY,
        volumeDensity * height * planar.secondYY,
        volumeDensity * zFirstIntegral * planar.firstY,
        volumeDensity * zFirstIntegral * planar.firstX,
        volumeDensity * zFirstIntegral * planar.firstY,
        volumeDensity * zSecondIntegral * planar.area;
    if (!std::isfinite(moments.mass) || !moments.firstMoment.allFinite() ||
        !moments.secondMoment.allFinite())
        throw std::overflow_error("extruded boundary moments are non-finite");
    return physics::RawMassMoments10::fromMassMoments(moments);
}

} // namespace spin
