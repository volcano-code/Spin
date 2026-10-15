// Exact surface integrals and Reynolds coefficients of the emitted material.
// Geometry, topology and source motion remain the caller's registered inputs.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gmpxx.h>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using Point = std::array<mpq_class, 3>;
using Face = std::array<std::size_t, 3>;
constexpr std::array<std::array<int, 2>, 6> pairs{{{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}}};

mpq_class parseRational(const std::string &text) {
    const auto slash = text.find('/');
    if (slash != std::string::npos && text.find('/', slash + 1) != std::string::npos)
        throw std::invalid_argument("INVALID_RATIONAL_TOKEN");
    mpz_class numerator, denominator = 1;
    const auto first = text.substr(0, slash);
    if (first.empty() || mpz_set_str(numerator.get_mpz_t(), first.c_str(), 0))
        throw std::invalid_argument("INVALID_RATIONAL_NUMERATOR");
    if (slash != std::string::npos) {
        const auto second = text.substr(slash + 1);
        if (second.empty() || mpz_set_str(denominator.get_mpz_t(), second.c_str(), 0) || denominator == 0)
            throw std::invalid_argument("INVALID_RATIONAL_DENOMINATOR");
    }
    mpq_class value(numerator, denominator);
    value.canonicalize();
    return value;
}

Point areaVector(const Point &a, const Point &b, const Point &c) {
    Point u, v;
    for (int k = 0; k < 3; ++k) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; }
    return {u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]};
}

struct Boundary {
    std::vector<Point> points, areas;
    std::vector<Face> faces;
};

Boundary boundary(const std::vector<std::vector<std::string>> &rows, const py::array &input) {
    Boundary result;
    if (rows.size() < 4) throw std::invalid_argument("INVALID_VERTEX_COUNT");
    for (const auto &row : rows) {
        if (row.size() != 3) throw std::invalid_argument("INVALID_VERTEX_DIMENSION");
        result.points.push_back({parseRational(row[0]), parseRational(row[1]), parseRational(row[2])});
    }
    const std::string kind = py::str(input.dtype().attr("kind"));
    if (kind != "i" && kind != "u") throw std::invalid_argument("NONINTEGER_FACES");
    const auto f = py::array_t<long long, py::array::c_style | py::array::forcecast>::ensure(input);
    if (!f || f.ndim() != 2 || f.shape(1) != 3 || f.shape(0) < 4)
        throw std::invalid_argument("INVALID_FACE_DIMENSION");
    const auto triangles = f.unchecked<2>();
    std::map<std::pair<std::size_t,std::size_t>, std::pair<int,int>> edges;
    for (py::ssize_t i = 0; i < f.shape(0); ++i) {
        Face face;
        for (int j = 0; j < 3; ++j) {
            const auto index = triangles(i,j);
            if (index < 0 || static_cast<std::size_t>(index) >= rows.size())
                throw std::invalid_argument("INVALID_FACE_INDEX");
            face[j] = static_cast<std::size_t>(index);
        }
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            throw std::invalid_argument("REPEATED_FACE_VERTEX");
        Point area = areaVector(result.points[face[0]], result.points[face[1]], result.points[face[2]]);
        if (area[0] == 0 && area[1] == 0 && area[2] == 0)
            throw std::invalid_argument("EXACT_ZERO_AREA");
        result.faces.push_back(face);
        result.areas.push_back(std::move(area));
        for (int j = 0; j < 3; ++j) {
            auto u = face[j], v = face[(j+1)%3];
            auto &edge = edges[{std::min(u,v), std::max(u,v)}];
            ++edge.first; edge.second += u < v ? 1 : -1;
        }
    }
    for (const auto &entry : edges)
        if (entry.second.first != 2 || entry.second.second != 0)
            throw std::invalid_argument("OPEN_NONMANIFOLD_OR_MISORIENTED_BOUNDARY");
    return result;
}

struct RationalSum {
    std::vector<std::optional<mpq_class>> partials;
    void add(mpq_class value) {
        std::size_t level = 0;
        while (level < partials.size() && partials[level]) {
            value += *partials[level]; partials[level].reset(); ++level;
        }
        if (level == partials.size()) partials.emplace_back(std::move(value));
        else partials[level] = std::move(value);
    }
    mpq_class total() const {
        mpq_class value = 0;
        for (const auto &partial : partials) if (partial) value += *partial;
        return value;
    }
};

struct RoundToDouble {
    double operator()(const mpq_class &value, bool preserveNonzero = false) {
        if (value == 0) return 0.;
        const bool negative = value < 0;
        mpz_class numerator = value.get_num(), denominator = value.get_den();
        if (negative) numerator = -numerator;
        long exponent = static_cast<long>(mpz_sizeinbase(numerator.get_mpz_t(),2)) -
                        static_cast<long>(mpz_sizeinbase(denominator.get_mpz_t(),2));
        if (exponent > 1024) throw std::invalid_argument("RATIONAL_COEFFICIENT_NOT_REPRESENTABLE");
        if (exponent < -1075) {
            if (preserveNonzero) throw std::invalid_argument("RATIONAL_COEFFICIENT_NOT_REPRESENTABLE");
            return negative ? -0. : 0.;
        }
        if (exponent >= 0) {
            mpz_class threshold = denominator << exponent;
            if (numerator < threshold) --exponent;
        } else {
            mpz_class scaled = numerator << -exponent;
            if (scaled < denominator) --exponent;
        }
        const long spacing = exponent >= -1022 ? exponent - 52 : -1074;
        if (spacing < 0) numerator <<= -spacing;
        else denominator <<= spacing;
        mpz_class significand, remainder;
        mpz_tdiv_qr(significand.get_mpz_t(), remainder.get_mpz_t(), numerator.get_mpz_t(), denominator.get_mpz_t());
        remainder *= 2;
        if (remainder > denominator || (remainder == denominator && mpz_odd_p(significand.get_mpz_t())))
            ++significand;
        double rounded = std::ldexp(significand.get_d(), static_cast<int>(spacing));
        if (negative) rounded = -rounded;
        if (!std::isfinite(rounded) || (preserveNonzero && value != 0 && rounded == 0))
            throw std::invalid_argument("RATIONAL_COEFFICIENT_NOT_REPRESENTABLE");
        return rounded;
    }
};

mpq_class densityValue(double density) {
    if (!std::isfinite(density) || density <= 0) throw std::invalid_argument("INVALID_DENSITY");
    return mpq_class(density);
}

py::array_t<double> integrateBoundary(const Boundary &mesh, const mpq_class &rho) {
    std::array<RationalSum,10> moments;
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        const auto &face = mesh.faces[i]; const auto &area = mesh.areas[i];
        Point sums{{0,0,0}};
        std::array<mpq_class,6> second{};
        for (const auto index : face) {
            const auto &v = mesh.points[index];
            for (int p = 0; p < 3; ++p) {
                sums[p] += v[p];
            }
            for (int m = 0; m < 6; ++m) second[m] += v[pairs[m][0]]*v[pairs[m][1]];
        }
        // The exact area vector is orthogonal to both triangle edges, so
        // area.dot(x) is the same K at all three vertices. Substituting that
        // identity in the original divergence face integrals cancels their
        // repeated cubic terms exactly: K/6, K*sum(x)/24 and
        // K*(sum(x)sum(y)+sum(x*y))/120. No geometric approximation occurs.
        mpq_class constant = 0;
        for (int c = 0; c < 3; ++c) constant += area[c]*mesh.points[face[0]][c];
        moments[0].add(constant/6);
        for (int p = 0; p < 3; ++p) moments[p+1].add(constant*sums[p]/24);
        for (int m = 0; m < 6; ++m) {
            const int p = pairs[m][0], q = pairs[m][1];
            moments[m+4].add(constant*(sums[p]*sums[q]+second[m])/120);
        }
    }
    py::array_t<double> output(10); auto out = output.mutable_unchecked<1>();
    RoundToDouble rounded;
    for (int m = 0; m < 10; ++m) out(m) = rounded(moments[m].total()*rho);
    if (out(0) <= 0) throw std::invalid_argument("NONPOSITIVE_MASS");
    return output;
}

py::array_t<double> integrateMoments(const std::vector<std::vector<std::string>> &vertices,
                                    const py::array &faces, double density) {
    const auto rho = densityValue(density);
    return integrateBoundary(boundary(vertices, faces), rho);
}

py::array_t<double> selectedTransport(const Boundary &mesh, const mpq_class &rho,
                                     const std::vector<std::size_t> &selected) {
    py::array_t<double> output(std::vector<py::ssize_t>{static_cast<py::ssize_t>(selected.size()),3,10,3});
    auto out = output.mutable_unchecked<4>(); RoundToDouble rounded;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto &face = mesh.faces[selected[i]]; const auto &area = mesh.areas[selected[i]];
        Point sums{{0,0,0}}; std::array<mpq_class,6> second{};
        for (auto index : face) {
            const auto &v = mesh.points[index];
            for (int c = 0; c < 3; ++c) sums[c] += v[c];
            for (int m = 0; m < 6; ++m) second[m] += v[pairs[m][0]]*v[pairs[m][1]];
        }
        for (int j = 0; j < 3; ++j) {
            const auto &v = mesh.points[face[j]];
            std::array<mpq_class,10> expectations;
            expectations[0] = mpq_class(1,3);
            for (int c = 0; c < 3; ++c) expectations[c+1] = (sums[c]+v[c])/12;
            for (int m = 0; m < 6; ++m) {
                const int p = pairs[m][0], q = pairs[m][1];
                expectations[m+4] = (sums[p]*sums[q]+second[m]+v[p]*sums[q]+v[q]*sums[p]+2*v[p]*v[q])/60;
            }
            for (int m = 0; m < 10; ++m)
                for (int c = 0; c < 3; ++c) out(i,j,m,c) = rounded(rho*area[c]*expectations[m]/2, true);
        }
    }
    return output;
}

py::array_t<double> transportCoefficients(const std::vector<std::vector<std::string>> &vertices,
                                         const py::array &faces, double density) {
    const auto rho = densityValue(density);
    const auto mesh = boundary(vertices, faces);
    std::vector<std::size_t> selected;
    selected.reserve(mesh.faces.size());
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) selected.push_back(i);
    return selectedTransport(mesh, rho, selected);
}

std::vector<std::size_t> selectedIndices(const Boundary &mesh,
                                       const std::vector<long long> &indices) {
    std::vector<std::size_t> selected;
    selected.reserve(indices.size());
    long long previous = -1;
    for (const auto index : indices) {
        if (index < 0 || static_cast<std::size_t>(index) >= mesh.faces.size() || index <= previous)
            throw std::invalid_argument("INVALID_ORDERED_TRANSPORT_FACE_SELECTION");
        previous = index;
        selected.push_back(static_cast<std::size_t>(index));
    }
    return selected;
}

py::array_t<double> transportCoefficientsSelected(
        const std::vector<std::vector<std::string>> &vertices,
        const py::array &faces, double density, const std::vector<long long> &indices) {
    const auto rho = densityValue(density);
    // Always validate every face, including stationary faces and an empty selection.
    const auto mesh = boundary(vertices, faces);
    return selectedTransport(mesh, rho, selectedIndices(mesh, indices));
}

// Own all coordinates, connectivity and exact area vectors. No Python-owned
// arrays or caller-supplied verification flags survive construction. This is
// arithmetic validation only, not a replacement for geometric/source proofs.
class PreparedBoundary {
    const Boundary mesh_;
public:
    PreparedBoundary(const std::vector<std::vector<std::string>> &vertices,
                     const py::array &faces) : mesh_(boundary(vertices, faces)) {}
    py::array_t<double> integrate(double density) const {
        return integrateBoundary(mesh_, densityValue(density));
    }
    py::array_t<double> transport(double density,
                                 const std::vector<long long> &indices) const {
        const auto rho = densityValue(density);
        return selectedTransport(mesh_, rho, selectedIndices(mesh_, indices));
    }
};

PYBIND11_MODULE(_cherchi_rational_moments, module) {
    module.doc() = "Exact divergence integrals and Reynolds coefficients, with binary64 nearest-even rounding";
    module.def("integrate_moments", &integrateMoments);
    module.def("transport_coefficients", &transportCoefficients);
    module.def("transport_coefficients_selected", &transportCoefficientsSelected);
    py::class_<PreparedBoundary>(module, "PreparedBoundary")
        .def(py::init<const std::vector<std::vector<std::string>> &, const py::array &>())
        .def("integrate_moments", &PreparedBoundary::integrate)
        .def("transport_coefficients_selected", &PreparedBoundary::transport);
    module.def("round_to_binary64", [](const std::string &token) { return RoundToDouble()(parseRational(token)); });
    module.attr("arithmetic_identity") = "DIV_SURFACE_V1_REYNOLDS_V1_GMP_BINARY64_NEAREST_EVEN_V1";
    module.attr("arithmetic_implementation") = "LOCAL_PLANAR_DIV_REDUCTION_V2";
    module.attr("transport_selection_identity") = "FULL_BOUNDARY_CHECKED_ORDERED_FACE_SELECTION_V1";
    module.attr("gmp_version") = gmp_version;
}
