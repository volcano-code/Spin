#include "deformation/LaplacianEditingEnergy.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


namespace deformation
{
namespace
{


Eigen::VectorXd flattenVertexRows(
    const Eigen::MatrixXd& vertices
)
{
    Eigen::VectorXd flattened(3 * vertices.rows());
    for (Eigen::Index vertex = 0;
         vertex < vertices.rows();
         ++vertex)
    {
        flattened.segment<3>(3 * vertex) =
            vertices.row(vertex).transpose();
    }
    return flattened;
}


void validateMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    if (
        vertices.rows() == 0 ||
        vertices.cols() != 3 ||
        !vertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "Laplacian rest vertices must be a nonempty finite N-by-3 matrix"
        );
    }

    if (
        faces.rows() == 0 ||
        faces.cols() != 3 ||
        faces.minCoeff() < 0 ||
        faces.maxCoeff() >= vertices.rows()
    )
    {
        throw std::invalid_argument(
            "Laplacian faces must contain valid triangles"
        );
    }

    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        if (
            faces(face, 0) == faces(face, 1) ||
            faces(face, 1) == faces(face, 2) ||
            faces(face, 2) == faces(face, 0)
        )
        {
            throw std::invalid_argument(
                "Laplacian faces must not contain repeated vertices"
            );
        }
    }

    /* Each undirected edge must occur exactly twice for a closed 2-manifold. */
    std::map<std::pair<int, int>, int> edgeUses;
    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        for (int edge = 0; edge < 3; ++edge)
        {
            const int first = faces(face, edge);
            const int second = faces(face, (edge + 1) % 3);
            ++edgeUses[std::minmax(first, second)];
        }
    }
    for (const auto& edgeUse : edgeUses)
    {
        if (edgeUse.second != 2)
        {
            throw std::invalid_argument(
                "Laplacian editing requires a closed manifold triangle mesh"
            );
        }
    }
}


void assignSimilarityBlock(
    const Eigen::Vector3d& point,
    double translationWeight,
    Eigen::MatrixXd& matrix,
    Eigen::Index row
)
{
    matrix(row, 0) = point(0);
    matrix(row, 2) = point(2);
    matrix(row, 3) = -point(1);
    matrix(row, 4) = translationWeight;

    matrix(row + 1, 0) = point(1);
    matrix(row + 1, 1) = -point(2);
    matrix(row + 1, 3) = point(0);
    matrix(row + 1, 5) = translationWeight;

    matrix(row + 2, 0) = point(2);
    matrix(row + 2, 1) = point(1);
    matrix(row + 2, 2) = -point(0);
    matrix(row + 2, 6) = translationWeight;
}


Eigen::MatrixXd fullRankPseudoInverse(
    const Eigen::MatrixXd& matrix,
    Eigen::Index vertex
)
{
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        matrix,
        Eigen::ComputeThinU | Eigen::ComputeThinV
    );

    if (svd.singularValues().size() != 7)
    {
        throw std::invalid_argument(
            "implicit-similarity neighborhood is underdetermined"
        );
    }

    const double maximum = svd.singularValues()(0);
    const double tolerance =
        std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(matrix.rows(), matrix.cols())) *
        maximum;

    if (
        !std::isfinite(maximum) ||
        !(maximum > 0.0) ||
        svd.singularValues()(6) <= tolerance
    )
    {
        throw std::invalid_argument(
            "implicit-similarity neighborhood at vertex " +
            std::to_string(vertex) + " is rank deficient"
        );
    }

    Eigen::VectorXd reciprocals = svd.singularValues();
    for (Eigen::Index index = 0;
         index < reciprocals.size();
         ++index)
    {
        reciprocals(index) = 1.0 / reciprocals(index);
    }

    return
        svd.matrixV() * reciprocals.asDiagonal() *
        svd.matrixU().transpose();
}


} // namespace


LaplacianEditingEnergy::LaplacianEditingEnergy(
    const Eigen::MatrixXd& restVertices,
    const Eigen::MatrixXi& faces
)
    :
    vertexCount_(restVertices.rows()),
    residualOperator_(
        3 * restVertices.rows(),
        3 * restVertices.rows()
    )
{
    validateMesh(restVertices, faces);

    std::vector<std::set<int>> adjacency(
        static_cast<std::size_t>(restVertices.rows())
    );

    for (Eigen::Index face = 0; face < faces.rows(); ++face)
    {
        for (int edge = 0; edge < 3; ++edge)
        {
            const int first = faces(face, edge);
            const int second = faces(face, (edge + 1) % 3);
            adjacency[static_cast<std::size_t>(first)].insert(second);
            adjacency[static_cast<std::size_t>(second)].insert(first);
        }
    }

    /*
     * Centering is a conditioning improvement over the beta implementation.
     * Translation is already part of the local seven-parameter model, so it
     * leaves the mathematical residual operator unchanged up to roundoff.
     */
    Eigen::MatrixXd centered = restVertices;
    centered.rowwise() -= restVertices.colwise().mean();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(
        static_cast<std::size_t>(27 * restVertices.rows())
    );

    for (Eigen::Index vertex = 0;
         vertex < restVertices.rows();
         ++vertex)
    {
        const std::set<int>& neighborSet =
            adjacency[static_cast<std::size_t>(vertex)];

        if (neighborSet.size() < 3U)
        {
            throw std::invalid_argument(
                "implicit-similarity Laplacian requires at least three neighbors"
            );
        }

        std::vector<int> localVertices;
        localVertices.reserve(neighborSet.size() + 1U);
        localVertices.push_back(static_cast<int>(vertex));
        localVertices.insert(
            localVertices.end(),
            neighborSet.begin(),
            neighborSet.end()
        );

        const Eigen::Index localCount =
            static_cast<Eigen::Index>(localVertices.size());
        Eigen::MatrixXd similarity = Eigen::MatrixXd::Zero(
            3 * localCount,
            7
        );
        Eigen::MatrixXd umbrella = Eigen::MatrixXd::Zero(
            3,
            3 * localCount
        );

        assignSimilarityBlock(
            centered.row(vertex).transpose(),
            1.0,
            similarity,
            0
        );
        umbrella.block<3, 3>(0, 0).setIdentity();

        const double neighborWeight =
            1.0 / static_cast<double>(neighborSet.size());

        for (Eigen::Index local = 1;
             local < localCount;
             ++local)
        {
            const int neighbor =
                localVertices[static_cast<std::size_t>(local)];
            assignSimilarityBlock(
                centered.row(neighbor).transpose(),
                1.0,
                similarity,
                3 * local
            );
            umbrella.block<3, 3>(0, 3 * local) =
                -neighborWeight * Eigen::Matrix3d::Identity();
        }

        const Eigen::Vector3d differential =
            umbrella * similarity.col(0);
        Eigen::MatrixXd differentialSimilarity =
            Eigen::MatrixXd::Zero(3, 7);
        assignSimilarityBlock(
            differential,
            0.0,
            differentialSimilarity,
            0
        );

        const Eigen::MatrixXd pseudoInverse =
            fullRankPseudoInverse(similarity, vertex);
        const Eigen::MatrixXd localResidual =
            differentialSimilarity * pseudoInverse - umbrella;

        if (!localResidual.allFinite())
        {
            throw std::domain_error(
                "implicit-similarity Laplacian produced non-finite coefficients"
            );
        }

        for (Eigen::Index local = 0;
             local < localCount;
             ++local)
        {
            const int global =
                localVertices[static_cast<std::size_t>(local)];
            for (Eigen::Index row = 0; row < 3; ++row)
            {
                for (Eigen::Index column = 0; column < 3; ++column)
                {
                    const double value =
                        localResidual(row, 3 * local + column);
                    if (value != 0.0)
                    {
                        triplets.emplace_back(
                            3 * vertex + row,
                            3 * global + column,
                            value
                        );
                    }
                }
            }
        }
    }

    residualOperator_.setFromTriplets(
        triplets.begin(),
        triplets.end()
    );
    residualOperator_.makeCompressed();

    if (residualOperator_.nonZeros() == 0)
    {
        throw std::logic_error(
            "implicit-similarity Laplacian is unexpectedly empty"
        );
    }
}


const Eigen::SparseMatrix<double>&
LaplacianEditingEnergy::residualOperator() const noexcept
{
    return residualOperator_;
}


Eigen::SparseMatrix<double>
LaplacianEditingEnergy::quadraticMatrix() const
{
    Eigen::SparseMatrix<double> quadratic =
        residualOperator_.transpose() * residualOperator_;
    quadratic.makeCompressed();
    return quadratic;
}


double LaplacianEditingEnergy::energy(
    const Eigen::MatrixXd& vertices
) const
{
    if (
        vertices.rows() != vertexCount_ ||
        vertices.cols() != 3 ||
        !vertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "Laplacian evaluation vertices do not match the rest mesh"
        );
    }

    const Eigen::VectorXd residual =
        residualOperator_ * flattenVertexRows(vertices);
    const double value = 0.5 * residual.squaredNorm();

    if (!std::isfinite(value))
    {
        throw std::overflow_error(
            "Laplacian editing energy overflowed"
        );
    }
    return value;
}


Eigen::VectorXd LaplacianEditingEnergy::vertexGradient(
    const Eigen::MatrixXd& vertices
) const
{
    if (
        vertices.rows() != vertexCount_ ||
        vertices.cols() != 3 ||
        !vertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "Laplacian gradient vertices do not match the rest mesh"
        );
    }

    const Eigen::VectorXd positions = flattenVertexRows(vertices);
    const Eigen::VectorXd gradient =
        residualOperator_.transpose() *
        (residualOperator_ * positions);

    if (!gradient.allFinite())
    {
        throw std::overflow_error(
            "Laplacian editing gradient overflowed"
        );
    }
    return gradient;
}


} // namespace deformation
