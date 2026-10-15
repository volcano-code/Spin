#include "spin/PlanarMassMoments.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace spin
{
namespace
{

using WideMoments = std::array<long double, 6>;


double finiteDouble(long double value, const char* label)
{
    const double result = static_cast<double>(value);
    if (!std::isfinite(result))
        throw std::overflow_error(label);
    return result;
}


bool samePoint(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
{
    return a.x() == b.x() && a.y() == b.y();
}


WideMoments translateToWorld(const WideMoments& local,
                             long double originX,
                             long double originY)
{
    const long double area = local[0];
    const long double firstX = local[1];
    const long double firstY = local[2];
    return {
        area,
        firstX + originX * area,
        firstY + originY * area,
        local[3] + 2.0L * originX * firstX + originX * originX * area,
        local[4] + 2.0L * originY * firstY + originY * originY * area,
        local[5] + originX * firstY + originY * firstX + originX * originY * area
    };
}


} // namespace


PlanarMassMoments6 PlanarRawMoments6::weighted(double arealDensity) const
{
    if (!std::isfinite(arealDensity) || arealDensity <= 0.0)
        throw std::invalid_argument("planar areal density must be finite and positive");
    const double mass = finiteDouble(
        static_cast<long double>(area) * arealDensity,
        "density-weighted planar area overflowed double precision");
    if (!(mass > 0.0))
        throw std::overflow_error("density-weighted planar area underflowed double precision");
    return {
        mass,
        finiteDouble(static_cast<long double>(firstX) * arealDensity,
                     "density-weighted planar first moment overflowed double precision"),
        finiteDouble(static_cast<long double>(firstY) * arealDensity,
                     "density-weighted planar first moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondXX) * arealDensity,
                     "density-weighted planar second moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondYY) * arealDensity,
                     "density-weighted planar second moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondXY) * arealDensity,
                     "density-weighted planar mixed moment overflowed double precision")
    };
}


PlanarMassMoments6& PlanarMassMoments6::operator+=(const PlanarMassMoments6& other)
{
    const PlanarMassMoments6 sum{
        finiteDouble(static_cast<long double>(mass) + other.mass,
                     "accumulated planar mass overflowed double precision"),
        finiteDouble(static_cast<long double>(firstX) + other.firstX,
                     "accumulated planar first moment overflowed double precision"),
        finiteDouble(static_cast<long double>(firstY) + other.firstY,
                     "accumulated planar first moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondXX) + other.secondXX,
                     "accumulated planar second moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondYY) + other.secondYY,
                     "accumulated planar second moment overflowed double precision"),
        finiteDouble(static_cast<long double>(secondXY) + other.secondXY,
                     "accumulated planar mixed moment overflowed double precision")
    };
    *this = sum;
    return *this;
}


Eigen::Matrix<double, 6, 1> PlanarMassMoments6::vector() const
{
    Eigen::Matrix<double, 6, 1> result;
    result << mass, firstX, firstY, secondXX, secondYY, secondXY;
    return result;
}


PlanarRawMoments6 integratePlanarBoundaryLoops(
    const std::vector<PlanarBoundaryLoop>& loops)
{
    if (loops.empty())
        throw std::invalid_argument("at least one planar boundary loop is required");

    const auto firstNonempty = std::find_if(loops.begin(), loops.end(),
        [](const PlanarBoundaryLoop& loop) { return !loop.empty(); });
    if (firstNonempty == loops.end())
        throw std::invalid_argument("planar boundary loops cannot all be empty");

    const long double originX = static_cast<long double>((*firstNonempty)[0].x());
    const long double originY = static_cast<long double>((*firstNonempty)[0].y());
    WideMoments local{};

    for (const PlanarBoundaryLoop& loop : loops)
    {
        std::size_t count = loop.size();
        if (count >= 2 && samePoint(loop.front(), loop.back()))
            --count;
        if (count < 3)
            throw std::invalid_argument("each planar boundary loop needs at least three vertices");

        WideMoments loopMoments{};
        for (std::size_t i = 0; i < count; ++i)
        {
            const Eigen::Vector2d& p = loop[i];
            const Eigen::Vector2d& q = loop[(i + 1) % count];
            if (!p.allFinite())
                throw std::invalid_argument("planar boundary coordinates must be finite");
            if (samePoint(p, q))
                throw std::invalid_argument("consecutive planar boundary vertices must differ");

            const long double x0 = static_cast<long double>(p.x()) - originX;
            const long double y0 = static_cast<long double>(p.y()) - originY;
            const long double x1 = static_cast<long double>(q.x()) - originX;
            const long double y1 = static_cast<long double>(q.y()) - originY;
            const long double cross = x0 * y1 - x1 * y0;

            loopMoments[0] += cross / 2.0L;
            loopMoments[1] += (x0 + x1) * cross / 6.0L;
            loopMoments[2] += (y0 + y1) * cross / 6.0L;
            loopMoments[3] += (x0*x0 + x0*x1 + x1*x1) * cross / 12.0L;
            loopMoments[4] += (y0*y0 + y0*y1 + y1*y1) * cross / 12.0L;
            loopMoments[5] += (2.0L*x0*y0 + x0*y1 + x1*y0 + 2.0L*x1*y1) * cross / 24.0L;
        }
        if (loopMoments[0] == 0.0L)
            throw std::invalid_argument("planar boundary loop must enclose nonzero signed area");
        for (std::size_t k = 0; k < local.size(); ++k)
            local[k] += loopMoments[k];
    }

    if (!(local[0] > 0.0L))
        throw std::invalid_argument("oriented planar region must have positive total area");

    const WideMoments world = translateToWorld(local, originX, originY);
    const double area = finiteDouble(world[0], "planar area overflowed double precision");
    if (!(area > 0.0))
        throw std::overflow_error("positive planar area underflowed double precision");
    return {
        area,
        finiteDouble(world[1], "planar first moment overflowed double precision"),
        finiteDouble(world[2], "planar first moment overflowed double precision"),
        finiteDouble(world[3], "planar second moment overflowed double precision"),
        finiteDouble(world[4], "planar second moment overflowed double precision"),
        finiteDouble(world[5], "planar mixed moment overflowed double precision")
    };
}


} // namespace spin
