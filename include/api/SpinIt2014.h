#pragma once

#include "spin2014/AdaptiveVoxelSpinSolver.h"

#include <Eigen/Core>

#include <memory>
#include <string>


namespace api
{


enum class SpinIt2014Evidence
{
    PaperFormulaReproduction,
    ReconstructedAsset
};


const char* spinIt2014EvidenceToken(SpinIt2014Evidence evidence) noexcept;


/** Immutable, hash-bound public product for adaptive Spin-It hollowing. */
class SpinIt2014Artifact
{
public:
    SpinIt2014Artifact(const SpinIt2014Artifact&) noexcept = default;
    SpinIt2014Artifact(SpinIt2014Artifact&&) = delete;
    SpinIt2014Artifact& operator=(const SpinIt2014Artifact&) noexcept = default;
    SpinIt2014Artifact& operator=(SpinIt2014Artifact&&) = delete;
    ~SpinIt2014Artifact();

    std::string statusToken() const;
    std::string evidenceToken() const;
    bool succeeded() const;
    bool geometrySealed() const;
    const spin2014::AdaptiveVoxelSpinResult& result() const;
    const Eigen::MatrixXd& outerVertices() const;
    const Eigen::MatrixXi& outerFaces() const;
    const Eigen::MatrixXd& designVertices() const;
    const Eigen::MatrixXi& designFaces() const;
    const Eigen::MatrixXd& materialVertices() const;
    const Eigen::MatrixXi& materialFaces() const;
    const Eigen::MatrixXd& cavityVertices() const;
    const Eigen::MatrixXi& cavityFaces() const;
    const std::string& diagnosticsJson() const;

    void saveOuterMesh(const std::string& filename) const;
    void saveDesignMesh(const std::string& filename) const;
    void saveMaterialMesh(const std::string& filename) const;
    void saveCavityMesh(const std::string& filename) const;
    void saveDiagnostics(const std::string& filename) const;

private:
    struct Storage;
    std::shared_ptr<const Storage> storage_;
    explicit SpinIt2014Artifact(std::shared_ptr<const Storage> storage);

    friend SpinIt2014Artifact solveSpinIt2014(
        const spin2014::AdaptiveVoxelSpinRequest&,
        SpinIt2014Evidence
    );
};


SpinIt2014Artifact solveSpinIt2014(
    const spin2014::AdaptiveVoxelSpinRequest& request,
    SpinIt2014Evidence evidence = SpinIt2014Evidence::PaperFormulaReproduction
);


} // namespace api
