#pragma once

#include "deformation/DeformationSystem.h"
#include "deformation/LaplacianEditingEnergy.h"
#include "deformation/LinearBlendSkinning.h"
#include "physics/BalanceObjective.h"

#include <Eigen/Core>

#include <stdexcept>


namespace deformation
{


/** A finite trial geometry that violates a physical parameter/mass domain. */
class InfeasibleDeformationState : public std::domain_error
{
public:
    using std::domain_error::domain_error;
};


/** One value-and-analytic-gradient evaluation of paper Eq. (11)--(14). */
struct DeformationEnergyEvaluation
{
    double totalEnergy = 0.0;
    double centerOfMassEnergy = 0.0;
    double shapeEnergy = 0.0;
    double unscaledShapeEnergy = 0.0;
    double mass = 0.0;
    Eigen::Vector3d centerOfMass = Eigen::Vector3d::Zero();
    Eigen::VectorXd gradient;
};


/**
 * Fixed-alpha deformation objective from paper Section 5.
 *
 * Voxel fill, cavity connectivity, targets, support, rotations, weights and
 * rest geometry are all frozen at construction. Only the unlocked-user
 * [translation, uniform-scale] vector q is varied. Carving/refill and target
 * refresh deliberately remain outside this R6 component.
 */
class DeformationObjective
{
public:
    DeformationObjective(
        const DeformationSystem& system,
        const VoxelGrid& fixedVoxelFill,
        const BalanceObjectives& objectives,
        double shapeScaleLambda,
        double density = 1.0,
        double minimumPositiveMass = 0.0
    );

    Eigen::VectorXd initialParameters() const;
    HandleSet makeHandles(const Eigen::VectorXd& parameters) const;

    DeformationEnergyEvaluation evaluate(
        const Eigen::VectorXd& parameters,
        double mu
    ) const;

    double shapeScaleLambda() const noexcept;
    const BalanceObjectives& objectives() const noexcept;

private:
    struct CompactCavityBinding
    {
        Eigen::MatrixXd restVertices;
        Eigen::MatrixXd weights;
        Eigen::MatrixXi faces;
    };

    static CompactCavityBinding makeCompactCavityBinding(
        const DeformationSystem& system,
        const VoxelGrid& fixedVoxelFill
    );

    DeformationObjective(
        const DeformationSystem& system,
        const VoxelGrid& fixedVoxelFill,
        const BalanceObjectives& objectives,
        double shapeScaleLambda,
        double density,
        double minimumPositiveMass,
        const CompactCavityBinding& compactCavity
    );

    HandleSet bindingHandles_;
    FreeParameterLbsModel outerModel_;
    FreeParameterLbsModel innerModel_;
    Eigen::MatrixXi outerFaces_;
    Eigen::MatrixXi innerFaces_;
    LaplacianEditingEnergy laplacianEnergy_;
    BalanceObjectives objectives_;
    double shapeScaleLambda_;
    double density_;
    double minimumPositiveMass_;
};


} // namespace deformation
