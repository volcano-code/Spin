#include "reachability/CarvingReachability.h"

#include <Highs.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace reachability
{
namespace
{

struct LinearSystem
{
    Eigen::MatrixXd C;
    Eigen::VectorXd d;
};

struct LpResult
{
    HighsModelStatus status = HighsModelStatus::kNotset;
    std::string statusText;
    Eigen::VectorXd values;
    bool apiSucceeded = false;
};

bool finite(const Eigen::MatrixXd& value)
{
    return value.array().isFinite().all();
}

double signedDoubleArea(const Eigen::MatrixXd& polygon)
{
    double area = 0.0;
    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Index next = (index + 1) % polygon.rows();
        area += polygon(index, 0) * polygon(next, 1) -
            polygon(next, 0) * polygon(index, 1);
    }
    return area;
}

LinearSystem buildPolygonSystem(
    const CarvingReachabilityInput& input,
    const Halfspaces2D& halfspaces
)
{
    LinearSystem system;
    const Eigen::Index cellCount = static_cast<Eigen::Index>(
        input.variableCells.size()
    );
    system.C.resize(halfspaces.A.rows(), cellCount);
    system.d.resize(halfspaces.A.rows());
    for (Eigen::Index row = 0; row < halfspaces.A.rows(); ++row)
    {
        const Eigen::Vector2d normal = halfspaces.A.row(row).transpose();
        const double offset = halfspaces.b(row);
        system.d(row) = offset * input.fixedMass -
            normal.dot(input.fixedFirstMoment2D);
        for (Eigen::Index column = 0; column < cellCount; ++column)
        {
            const ReachabilityCell& cell =
                input.variableCells[static_cast<std::size_t>(column)];
            system.C(row, column) =
                normal.dot(cell.firstMoment2D) - offset * cell.mass;
        }
    }
    return system;
}

LinearSystem appendUnitBounds(const LinearSystem& polygonSystem)
{
    const Eigen::Index polygonRows = polygonSystem.C.rows();
    const Eigen::Index cellCount = polygonSystem.C.cols();
    LinearSystem system;
    system.C = Eigen::MatrixXd::Zero(
        polygonRows + 2 * cellCount,
        cellCount
    );
    system.d = Eigen::VectorXd::Zero(system.C.rows());
    system.C.topRows(polygonRows) = polygonSystem.C;
    system.d.head(polygonRows) = polygonSystem.d;
    for (Eigen::Index column = 0; column < cellCount; ++column)
    {
        system.C(polygonRows + column, column) = 1.0;
        system.d(polygonRows + column) = 1.0;
        system.C(polygonRows + cellCount + column, column) = -1.0;
    }
    return system;
}

LpResult solveLp(
    const Eigen::MatrixXd& matrix,
    const Eigen::VectorXd& rowLower,
    const Eigen::VectorXd& rowUpper,
    const Eigen::VectorXd& cost,
    const Eigen::VectorXd& columnLower,
    const Eigen::VectorXd& columnUpper
)
{
    LpResult result;
    HighsModel model;
    model.lp_.num_col_ = static_cast<HighsInt>(matrix.cols());
    model.lp_.num_row_ = static_cast<HighsInt>(matrix.rows());
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.offset_ = 0.0;
    model.lp_.col_cost_.assign(cost.data(), cost.data() + cost.size());
    model.lp_.col_lower_.assign(
        columnLower.data(),
        columnLower.data() + columnLower.size()
    );
    model.lp_.col_upper_.assign(
        columnUpper.data(),
        columnUpper.data() + columnUpper.size()
    );
    model.lp_.row_lower_.assign(
        rowLower.data(),
        rowLower.data() + rowLower.size()
    );
    model.lp_.row_upper_.assign(
        rowUpper.data(),
        rowUpper.data() + rowUpper.size()
    );
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    // HighsSparseMatrix starts with a sentinel zero. Rebuild the complete
    // start array here; appending to the default would shift every column.
    model.lp_.a_matrix_.start_.clear();
    model.lp_.a_matrix_.start_.reserve(
        static_cast<std::size_t>(matrix.cols() + 1)
    );
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
    {
        model.lp_.a_matrix_.start_.push_back(
            static_cast<HighsInt>(model.lp_.a_matrix_.value_.size())
        );
        for (Eigen::Index row = 0; row < matrix.rows(); ++row)
        {
            const double value = matrix(row, column);
            if (value != 0.0)
            {
                model.lp_.a_matrix_.index_.push_back(
                    static_cast<HighsInt>(row)
                );
                model.lp_.a_matrix_.value_.push_back(value);
            }
        }
    }
    model.lp_.a_matrix_.start_.push_back(
        static_cast<HighsInt>(model.lp_.a_matrix_.value_.size())
    );

    Highs solver;
    bool optionsOk = true;
    optionsOk = optionsOk &&
        solver.setOptionValue("output_flag", false) == HighsStatus::kOk;
    optionsOk = optionsOk &&
        solver.setOptionValue("presolve", "on") == HighsStatus::kOk;
    optionsOk = optionsOk &&
        solver.setOptionValue("threads", static_cast<HighsInt>(1)) ==
            HighsStatus::kOk;
    optionsOk = optionsOk &&
        solver.setOptionValue("parallel", "off") == HighsStatus::kOk;
    optionsOk = optionsOk &&
        solver.setOptionValue("random_seed", static_cast<HighsInt>(0)) ==
            HighsStatus::kOk;
    if (!optionsOk || solver.passModel(model) == HighsStatus::kError)
    {
        result.statusText = "API_ERROR";
        return result;
    }
    const HighsStatus runStatus = solver.run();
    result.status = solver.getModelStatus();
    result.statusText = solver.modelStatusToString(result.status);
    result.apiSucceeded = runStatus != HighsStatus::kError;
    if (result.status == HighsModelStatus::kOptimal)
    {
        const std::vector<double>& values = solver.getSolution().col_value;
        result.values = Eigen::Map<const Eigen::VectorXd>(
            values.data(),
            static_cast<Eigen::Index>(values.size())
        );
    }
    return result;
}

double maximumPrimalViolation(
    const LinearSystem& system,
    const Eigen::VectorXd& retained
)
{
    double violation = 0.0;
    if (system.C.rows() > 0)
    {
        violation = std::max(
            violation,
            (system.C * retained - system.d).maxCoeff()
        );
    }
    if (retained.size() > 0)
    {
        violation = std::max(violation, -retained.minCoeff());
        violation = std::max(violation, retained.maxCoeff() - 1.0);
    }
    return std::max(0.0, violation);
}

void populatePhysicalWitness(
    const CarvingReachabilityInput& input,
    const LinearSystem& system,
    const Eigen::VectorXd& retained,
    CarvingReachabilityResult& result
)
{
    result.retainedFractions = retained;
    result.finalMass = input.fixedMass;
    result.finalFirstMoment2D = input.fixedFirstMoment2D;
    for (Eigen::Index index = 0; index < retained.size(); ++index)
    {
        const ReachabilityCell& cell =
            input.variableCells[static_cast<std::size_t>(index)];
        result.finalMass += retained(index) * cell.mass;
        result.finalFirstMoment2D +=
            retained(index) * cell.firstMoment2D;
    }
    result.finalCenter2D = result.finalFirstMoment2D / result.finalMass;
    result.maximumPrimalViolation =
        maximumPrimalViolation(system, retained);
}

} // namespace

Halfspaces2D safePolygonHalfspaces(const Eigen::MatrixXd& polygon)
{
    if (polygon.cols() != 2 || polygon.rows() < 3)
    {
        throw std::invalid_argument(
            "safe polygon must contain at least three 2D vertices"
        );
    }
    if (!finite(polygon))
    {
        throw std::invalid_argument("safe polygon contains non-finite values");
    }
    if (!(signedDoubleArea(polygon) > 0.0))
    {
        throw std::invalid_argument(
            "safe polygon vertices must be strictly counter-clockwise"
        );
    }

    Halfspaces2D result;
    result.A.resize(polygon.rows(), 2);
    result.b.resize(polygon.rows());
    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Index next = (index + 1) % polygon.rows();
        const Eigen::Vector2d edge =
            polygon.row(next).transpose() - polygon.row(index).transpose();
        const double length = edge.norm();
        if (!(length > 0.0) || !std::isfinite(length))
        {
            throw std::invalid_argument(
                "safe polygon contains a duplicate or invalid edge"
            );
        }
        const Eigen::Vector2d outward(edge.y(), -edge.x());
        result.A.row(index) = (outward / length).transpose();
        result.b(index) = result.A.row(index).dot(polygon.row(index));
    }

    for (Eigen::Index vertex = 0; vertex < polygon.rows(); ++vertex)
    {
        const Eigen::VectorXd violations =
            result.A * polygon.row(vertex).transpose() - result.b;
        if (violations.maxCoeff() > 1e-12)
        {
            throw std::invalid_argument(
                "safe polygon must be convex and counter-clockwise"
            );
        }
    }
    return result;
}

bool contains(
    const Halfspaces2D& halfspaces,
    const Eigen::Vector2d& point,
    double tolerance
)
{
    if (
        halfspaces.A.cols() != 2 ||
        halfspaces.A.rows() != halfspaces.b.size() ||
        !finite(halfspaces.A) ||
        !halfspaces.b.array().isFinite().all() ||
        !point.array().isFinite().all() ||
        !std::isfinite(tolerance) ||
        tolerance < 0.0
    )
    {
        throw std::invalid_argument("invalid halfspace containment query");
    }
    if (halfspaces.A.rows() == 0)
    {
        return false;
    }
    return (
        halfspaces.A * point - halfspaces.b
    ).maxCoeff() <= tolerance;
}

void validateInput(const CarvingReachabilityInput& input)
{
    if (!(input.fixedMass > 0.0) || !std::isfinite(input.fixedMass))
    {
        throw std::invalid_argument("fixed shell mass must be positive");
    }
    if (!input.fixedFirstMoment2D.array().isFinite().all())
    {
        throw std::invalid_argument(
            "fixed shell first moment contains non-finite values"
        );
    }
    if (
        !(input.feasibilityTolerance > 0.0) ||
        !std::isfinite(input.feasibilityTolerance)
    )
    {
        throw std::invalid_argument("feasibility tolerance must be positive");
    }
    static_cast<void>(safePolygonHalfspaces(input.safePolygon2D));

    std::unordered_set<int> identifiers;
    for (const ReachabilityCell& cell : input.variableCells)
    {
        if (
            cell.cellIndex < 0 ||
            cell.mass < 0.0 ||
            !std::isfinite(cell.mass) ||
            !cell.firstMoment2D.array().isFinite().all()
        )
        {
            throw std::invalid_argument("invalid variable cell moments");
        }
        if (!identifiers.insert(cell.cellIndex).second)
        {
            throw std::invalid_argument("duplicate variable cell identifier");
        }
    }
}

ReachabilityCertificate validateFarkasCertificate(
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& d,
    const Eigen::VectorXd& weights,
    double tolerance
)
{
    if (
        C.rows() == 0 ||
        C.rows() != d.size() ||
        C.rows() != weights.size() ||
        !finite(C) ||
        !d.array().isFinite().all() ||
        !weights.array().isFinite().all() ||
        !(tolerance > 0.0) ||
        !std::isfinite(tolerance)
    )
    {
        throw std::invalid_argument("invalid Farkas certificate inputs");
    }

    ReachabilityCertificate certificate;
    certificate.weights = weights;
    certificate.minimumWeight = weights.minCoeff();
    certificate.normalizationResidual = std::abs(weights.sum() - 1.0);
    certificate.stationarityInfinityNorm = C.cols() == 0
        ? 0.0
        : (C.transpose() * weights).lpNorm<Eigen::Infinity>();
    certificate.contradiction = d.dot(weights);
    certificate.valid =
        certificate.minimumWeight >= -tolerance &&
        certificate.normalizationResidual <= tolerance &&
        certificate.stationarityInfinityNorm <= tolerance &&
        certificate.contradiction < -tolerance;
    return certificate;
}

const char* toString(ReachabilityStatus status) noexcept
{
    switch (status)
    {
    case ReachabilityStatus::AlreadySafe: return "ALREADY_SAFE";
    case ReachabilityStatus::ProvenUnreachable: return "PROVEN_UNREACHABLE";
    case ReachabilityStatus::RelaxationReachable:
        return "RELAXATION_REACHABLE";
    case ReachabilityStatus::Indeterminate: return "INDETERMINATE";
    }
    return "INDETERMINATE";
}

CarvingReachabilityResult CarvingReachabilityAnalyzer::analyze(
    const CarvingReachabilityInput& input
) const
{
    validateInput(input);
    const Halfspaces2D halfspaces =
        safePolygonHalfspaces(input.safePolygon2D);
    const LinearSystem system = buildPolygonSystem(input, halfspaces);
    const Eigen::Index cellCount = static_cast<Eigen::Index>(
        input.variableCells.size()
    );

    CarvingReachabilityResult result;
    result.inputCellCount = static_cast<int>(cellCount);
    const Eigen::VectorXd full = Eigen::VectorXd::Ones(cellCount);
    if (maximumPrimalViolation(system, full) <= input.feasibilityTolerance)
    {
        result.status = ReachabilityStatus::AlreadySafe;
        result.primalSolverStatus = "NOT_NEEDED_FULL_MATERIAL_SAFE";
        populatePhysicalWitness(input, system, full, result);
        return result;
    }

    const LpResult primal = solveLp(
        system.C,
        Eigen::VectorXd::Constant(system.C.rows(), -1e30),
        system.d,
        Eigen::VectorXd::Zero(cellCount),
        Eigen::VectorXd::Zero(cellCount),
        Eigen::VectorXd::Ones(cellCount)
    );
    result.primalSolverStatus = primal.statusText;
    if (
        primal.apiSucceeded &&
        primal.status == HighsModelStatus::kOptimal &&
        primal.values.size() == cellCount
    )
    {
        populatePhysicalWitness(input, system, primal.values, result);
        if (
            result.maximumPrimalViolation <= input.feasibilityTolerance &&
            contains(
                halfspaces,
                result.finalCenter2D,
                input.feasibilityTolerance
            )
        )
        {
            result.status = ReachabilityStatus::RelaxationReachable;
        }
        return result;
    }
    if (
        primal.apiSucceeded &&
        primal.status == HighsModelStatus::kInfeasible
    )
    {
        const LinearSystem fullSystem = appendUnitBounds(system);
        const Eigen::Index proofRows = fullSystem.C.rows();
        Eigen::MatrixXd certificateMatrix = Eigen::MatrixXd::Zero(
            cellCount + 1,
            proofRows
        );
        if (cellCount > 0)
        {
            certificateMatrix.topRows(cellCount) =
                fullSystem.C.transpose();
        }
        certificateMatrix.row(cellCount).setOnes();
        Eigen::VectorXd certificateRows = Eigen::VectorXd::Zero(
            cellCount + 1
        );
        certificateRows(cellCount) = 1.0;
        const LpResult certificateLp = solveLp(
            certificateMatrix,
            certificateRows,
            certificateRows,
            fullSystem.d,
            Eigen::VectorXd::Zero(proofRows),
            Eigen::VectorXd::Constant(proofRows, 1e30)
        );
        result.certificateSolverStatus = certificateLp.statusText;
        if (
            certificateLp.apiSucceeded &&
            certificateLp.status == HighsModelStatus::kOptimal &&
            certificateLp.values.size() == proofRows
        )
        {
            result.certificate = validateFarkasCertificate(
                fullSystem.C,
                fullSystem.d,
                certificateLp.values,
                input.feasibilityTolerance
            );
            if (result.certificate.valid)
            {
                result.status = ReachabilityStatus::ProvenUnreachable;
                return result;
            }
        }
    }
    result.status = ReachabilityStatus::Indeterminate;
    return result;
}

} // namespace reachability
