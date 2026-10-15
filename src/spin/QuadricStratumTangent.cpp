#include "spin/QuadricStratumTangent.h"

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

namespace {
using Vector = spin::QuadricCoefficients10::Vector;

Eigen::Matrix<double, 1, 10> quadraticCoefficientRow(const Eigen::Vector3d& n) {
    Eigen::Matrix<double, 1, 10> row = Eigen::Matrix<double, 1, 10>::Zero();
    row(4) = n.x() * n.x();
    row(5) = n.y() * n.y();
    row(6) = n.z() * n.z();
    row(7) = 2.0 * n.y() * n.z();
    row(8) = 2.0 * n.x() * n.z();
    row(9) = 2.0 * n.x() * n.y();
    return row;
}

Eigen::Matrix<double, 1, 10> mixedQuadraticCoefficientRow(
    const Eigen::Vector3d& first,const Eigen::Vector3d& second) {
    Eigen::Matrix<double,1,10> row=Eigen::Matrix<double,1,10>::Zero();
    for(int axis=0;axis<3;++axis)row(4+axis)=first(axis)*second(axis);
    row(7)=first.y()*second.z()+first.z()*second.y();
    row(8)=first.x()*second.z()+first.z()*second.x();
    row(9)=first.x()*second.y()+first.y()*second.x();
    return row;
}

Eigen::Matrix3d coefficientMatrix(int index) {
    Eigen::Matrix3d out = Eigen::Matrix3d::Zero();
    switch (index) {
    case 4: out(0, 0) = 1.; break;
    case 5: out(1, 1) = 1.; break;
    case 6: out(2, 2) = 1.; break;
    // QuadricCoefficients10::quadraticMatrix stores each cross coefficient
    // directly in both symmetric off-diagonal entries. The polynomial then
    // contains 2*q_k*x_i*x_j, so dA/dq_k has two unit entries (not halves).
    case 7: out(1, 2) = out(2, 1) = 1.; break;
    case 8: out(0, 2) = out(2, 0) = 1.; break;
    case 9: out(0, 1) = out(1, 0) = 1.; break;
    default: break;
    }
    return out;
}

void normalizeRows(Eigen::MatrixXd* rows) {
    for (Eigen::Index i = 0; i < rows->rows(); ++i) {
        const double norm = rows->row(i).norm();
        if (norm > 0.0) rows->row(i) /= norm;
    }
}
}

namespace spin {

bool QuadricStratumTangent::contains(
    const QuadricCoefficients10::Vector& direction, double tolerance) const {
    return defined && direction.allFinite() &&
        (constraintRows * direction).norm() <= tolerance;
}

QuadricStratumTangent evaluateQuadricStratumTangent(
    const QuadricCoefficients10& quadric) {
    QuadricStratumTangent result;
    result.subtype = classifyQuadricSurfaceSubtype(quadric);
    if (result.subtype == QuadricSurfaceSubtype::Cone ||
        result.subtype == QuadricSurfaceSubtype::IntersectingPlanes) {
        throw std::invalid_argument(
            std::string(quadricSurfaceSubtypeToken(result.subtype)) +
            " fixed-stratum tangent is not implemented; refusing to return a "
            "normalization-only ambient basis");
    }
    Eigen::MatrixXd rows(1, 10);
    rows.row(0) = quadric.values().transpose();

    if (result.subtype == QuadricSurfaceSubtype::Plane) {
        Eigen::MatrixXd extended(7, 10);
        extended.setZero();
        extended.row(0) = rows.row(0);
        for (int i = 0; i < 6; ++i) extended(i + 1, i + 4) = 1.0;
        rows = extended;
    } else if (result.subtype == QuadricSurfaceSubtype::Paraboloid ||
               result.subtype == QuadricSurfaceSubtype::EllipticCylinder ||
               result.subtype == QuadricSurfaceSubtype::HyperbolicCylinder) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(
            quadric.quadraticMatrix());
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("quadric stratum eigendecomposition failed");
        const auto eigenvalues = solver.eigenvalues();
        int nullIndex = 0;
        eigenvalues.cwiseAbs().minCoeff(&nullIndex);
        const Eigen::Vector3d nullVector = solver.eigenvectors().col(nullIndex);
        Eigen::Matrix<double, 1, 10> rankRow = quadraticCoefficientRow(nullVector);
        Eigen::MatrixXd extended(2, 10);
        extended.row(0) = rows.row(0);
        extended.row(1) = rankRow;
        rows = extended;

        if (result.subtype == QuadricSurfaceSubtype::EllipticCylinder ||
            result.subtype == QuadricSurfaceSubtype::HyperbolicCylinder) {
            const Eigen::Vector3d b = quadric.linearVector();
            Eigen::Matrix<double, 1, 10> cylinderRow =
                Eigen::Matrix<double, 1, 10>::Zero();
            cylinderRow.segment<3>(1) = nullVector.transpose();
            for (int coefficient = 4; coefficient <= 9; ++coefficient) {
                const Eigen::Matrix3d dA = coefficientMatrix(coefficient);
                double derivative = 0.0;
                for (int i = 0; i < 3; ++i) {
                    if (i == nullIndex) continue;
                    const double lambda = eigenvalues(i);
                    if (std::abs(lambda) <= 1e-12)
                        throw std::runtime_error("elliptic-cylinder null eigenspace is not simple");
                    const Eigen::Vector3d dNull =
                        -solver.eigenvectors().col(i) *
                        (solver.eigenvectors().col(i).dot(dA * nullVector)) / lambda;
                    derivative += dNull.dot(b);
                }
                cylinderRow(coefficient) = derivative;
            }
            Eigen::MatrixXd withCylinder(3, 10);
            withCylinder.topRows(2) = rows;
            withCylinder.row(2) = cylinderRow;
            rows = withCylinder;
        }
    } else if (result.subtype == QuadricSurfaceSubtype::ParabolicCylinder) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(
            quadric.quadraticMatrix());
        if(solver.info()!=Eigen::Success)
            throw std::runtime_error("parabolic-cylinder eigendecomposition failed");
        const auto eigenvalues=solver.eigenvalues();
        int curvedIndex=0;
        eigenvalues.cwiseAbs().maxCoeff(&curvedIndex);
        const double scale=std::max(1.0,eigenvalues.cwiseAbs().maxCoeff());
        if(std::abs(eigenvalues(curvedIndex))<=1e-12*scale)
            throw std::invalid_argument("parabolic-cylinder curvature is singular");
        std::array<int,2> nullIndices{};int cursor=0;
        for(int i=0;i<3;++i)if(i!=curvedIndex)nullIndices[cursor++]=i;
        const Eigen::Vector3d first=solver.eigenvectors().col(nullIndices[0]);
        const Eigen::Vector3d second=solver.eigenvectors().col(nullIndices[1]);
        Eigen::MatrixXd extended(4,10);
        extended.row(0)=rows.row(0);
        extended.row(1)=quadraticCoefficientRow(first);
        extended.row(2)=quadraticCoefficientRow(second);
        extended.row(3)=mixedQuadraticCoefficientRow(first,second);
        rows=extended;
    }

    normalizeRows(&rows);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        rows, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const double scale = std::max(1.0, rows.norm());
    const int rank = static_cast<int>((svd.singularValues().array() >
                                       1e-11 * scale).count());
    result.constraintRows = rows;
    result.tangentDimension = 10 - rank;
    result.basis = svd.matrixV().rightCols(result.tangentDimension);
    result.defined = true;
    std::ostringstream message;
    message << "FIXED_STRATUM_TANGENT:" << quadricSurfaceSubtypeToken(result.subtype)
            << ":constraints=" << rank << ":dimension=" << result.tangentDimension;
    result.diagnostic = message.str();
    return result;
}

} // namespace spin
