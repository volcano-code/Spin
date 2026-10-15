#include "manufacturing/ManufacturingGeometry.h"
#include "physics/MassProperties.h"
#include "spin/QuadricMomentEvaluator.h"
#include "spin/QuadricDiscreteMomentDifferential.h"
#include "spin/FullRankMomentEvaluator.h"
#include "spin/CanonicalFullRankHollowing.h"
#include "spin/QuadricSourceMotion.h"
#include "verification/ExactBoundaryAudit.h"
#include "verification/ExactRationalCoordinates.h"
#include "quadric_certified_empty_cases.h"

#include <igl/piecewise_constant_winding_number.h>
#include <igl/winding_number.h>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <chrono>


namespace
{


spin::TriangleMesh cube(double minimum, double maximum)
{
    spin::TriangleMesh mesh;
    mesh.vertices.resize(8, 3);
    mesh.vertices <<
        minimum, minimum, minimum,
        maximum, minimum, minimum,
        maximum, maximum, minimum,
        minimum, maximum, minimum,
        minimum, minimum, maximum,
        maximum, minimum, maximum,
        maximum, maximum, maximum,
        minimum, maximum, maximum;
    mesh.faces.resize(12, 3);
    mesh.faces <<
        0, 2, 1, 0, 3, 2,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5;
    return mesh;
}

spin::QuadricCoefficients10 hyperbolicCylinderRetraction(
    const spin::QuadricCoefficients10& base,
    const Eigen::Matrix<double, 10, 1>& tangent,
    double step)
{
    // A straight line in the 10-D coefficient space leaves the rank-two
    // cylinder stratum at second order for tangent directions that rotate
    // its null axis. Project back to rank two and enforce b·n=0 so the
    // central differences follow the same smooth geometric stratum.
    Eigen::Matrix<double, 10, 1> values = base.values() + step * tangent;
    Eigen::Matrix3d quadratic;
    quadratic << values(4), values(9), values(8),
        values(9), values(5), values(7),
        values(8), values(7), values(6);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(quadratic);
    if (eigensolver.info() != Eigen::Success)
        throw std::runtime_error("test hyperbolic-cylinder retraction eigensolve failed");
    int nullIndex = 0;
    eigensolver.eigenvalues().cwiseAbs().minCoeff(&nullIndex);
    Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
    eigenvalues(nullIndex) = 0.0;
    quadratic = eigensolver.eigenvectors() * eigenvalues.asDiagonal() *
        eigensolver.eigenvectors().transpose();
    const Eigen::Vector3d nullAxis = eigensolver.eigenvectors().col(nullIndex);
    Eigen::Vector3d linear = values.segment<3>(1);
    linear -= nullAxis.dot(linear) * nullAxis;
    values.segment<3>(1) = linear;
    values(4) = quadratic(0, 0);
    values(5) = quadratic(1, 1);
    values(6) = quadratic(2, 2);
    values(7) = quadratic(1, 2);
    values(8) = quadratic(0, 2);
    values(9) = quadratic(0, 1);
    return spin::QuadricCoefficients10::fromUnnormalized(values);
}


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-moment-evaluator "
              << message << '\n';
    return EXIT_FAILURE;
}

spin::SpinDesignRequest request();

int degenerateBooleanWitnessRegression()
{
    const spin::SpinDesignRequest input = request();
    spin::QuadricMomentEvaluatorSettings settings;
    settings.triangulation.icosphereSubdivisions = 2;
    settings.triangulation.truncationRadius = 8.0;
    settings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    spin::QuadricCoefficients10::Vector values;
    values << 5.8039689277516467e-15, 1.6851299752289185e-11,
        4.6761148234240197e-13, -1.5926745456717686e-11,
        0.40440691574567661, 0.40441404629432803,
        -0.82030745800395988, -2.1774792968580463e-12,
        7.0187201890690117e-11, -6.3251619031823912e-15;
    try
    {
        (void)spin::evaluateQuadricMaterialMoments(
            input,
            spin::QuadricCoefficients10::fromUnnormalized(values),
            settings
        );
        return fail("captured near-boundary quadric unexpectedly passed");
    }
    catch (const std::exception& error)
    {
        const std::string diagnostic = error.what();
        std::cout << "[WITNESS] " << diagnostic << '\n';
        if (diagnostic.find("nondegenerate=0") == std::string::npos ||
            diagnostic.find("minimum_area_face=") == std::string::npos ||
            diagnostic.find("area_witness_vertices=") == std::string::npos ||
            diagnostic.find("area_witness_points=") == std::string::npos)
        {
            return fail("captured invalid Boolean lacks degenerate-face witness: " +
                        diagnostic);
        }
    }
    std::cout << "[PASS] invalid Boolean preserves degenerate-face witness\n";
    return EXIT_SUCCESS;
}

int sourceWindingRegression()
{
    const spin::SpinDesignRequest input = request();
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions = 7;
    settings.maximumRadiusDoublings = 4;
    Eigen::Matrix<double, 10, 1> values;
    values << 0.00012209761677822807, 0.00013586216695654453,
        0.00004062482246953246, -0.015203731818126223,
        -0.0029325147146796226, 0.0011128695928818108,
        -0.9998444518008077, 0.0028733175213401523,
        0.007853543165293901, 0.0003339312833061543;
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(values);
    const auto triangulation = spin::triangulateQuadric(
        quadric, input.designDomain.vertices, settings);
    if (!igl::piecewise_constant_winding_number(
            triangulation.closedBoundary.faces))
    {
        return fail("captured full-rank quadric source is closed but not PWN");
    }
    const auto audit = manufacturing::auditTriangleBoundary(
        triangulation.closedBoundary.vertices,
        triangulation.closedBoundary.faces,
        manufacturing::BoundaryOrientation::OutwardPositive);
    if (!audit.strictBoundaryValid())
    {
        return fail(std::string("captured full-rank quadric source is not outward PWN: ") +
            "closed=" + std::to_string(audit.closed) +
            " consistent=" + std::to_string(audit.consistentlyOriented) +
            " components_valid=" + std::to_string(audit.componentOrientationsValid) +
            " intersections=" + std::to_string(audit.selfIntersections.intersectingPairCount));
    }
    spin::QuadricMomentEvaluatorSettings evaluatorSettings;
    evaluatorSettings.triangulation = settings;
    evaluatorSettings.capRadiusPolicy =
        spin::QuadricCapRadiusPolicy::GrowToCertifiedV1;
    evaluatorSettings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    try
    {
        const auto evaluation = spin::evaluateQuadricMaterialMoments(
            input, quadric, evaluatorSettings);
        const auto materialAudit = manufacturing::auditTriangleBoundary(
            evaluation.materialBoundary.vertices,
            evaluation.materialBoundary.faces,
            manufacturing::BoundaryOrientation::OutwardPositive);
        if (!evaluation.derivativeDefined || !materialAudit.strictBoundaryValid() ||
            evaluation.variableMomentsTarget.mass() <= 0.0 ||
            evaluation.geometryHash.size() != 64 || evaluation.topologyHash.size() != 64)
        {
            return fail("captured quadric Boolean returned an invalid material result");
        }
    }
    catch (const std::runtime_error& error)
    {
        const std::string diagnostic = error.what();
        if (diagnostic.find("quadric material Boolean is not a positive closed manifold:") != 0 ||
            diagnostic.find("degenerate_faces=") == std::string::npos ||
            diagnostic.find("BOOLEAN_ARRANGEMENT_PWN_FAILURE") != std::string::npos)
        {
            return fail("captured normalized quadric failed outside the strict geometric "
                "admission gate: " + diagnostic);
        }
        std::cout << "[WITNESS] normalized Boolean source passes PWN; output remains "
            "rejected by strict geometry admission: " << diagnostic << '\n';
    }
    std::cout << "[PASS] captured full-rank quadric source has outward PWN winding\n";
    return EXIT_SUCCESS;
}

int sourceInteriorSideWitnessRegression()
{
    // Captured Newton trial 58 from the deterministic Cherchi seed-0 run.
    // Its closed-source volume centroid is outside the source and has q>0,
    // while the consistently outward-oriented quadric interface bounds the
    // q<0 side. This guards against using a non-convex volume centroid as an
    // interior-side witness.
    Eigen::Matrix<double, 10, 1> values;
    values << 0.0005536453405000541, 0.00031719656565886746,
        2.5047430820147344e-05, 0.026524212134560434,
        0.022029948311259807, 0.022375274659397464,
        -0.9986765016162581, 0.02355956491378978,
        -0.019986897884358167, 0.0008959918010455094;
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(values);
    const auto input = request();
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions = 5;
    settings.truncationRadius = 6.928203230275509;
    const auto triangulation = spin::triangulateQuadric(
        quadric, input.designDomain.vertices, settings);
    const auto moments = computeSurfaceMassMoments(
        triangulation.closedBoundary.vertices,
        triangulation.closedBoundary.faces,
        1.0);
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    if (!moments.tryCenterOfMass(centroid, 1e-15) ||
        !(quadric.evaluate(centroid) > 0.0))
        return fail("captured source no longer reproduces its positive-q outside centroid");
    const double winding = igl::winding_number(
        triangulation.closedBoundary.vertices,
        triangulation.closedBoundary.faces,
        centroid);
    if (!std::isfinite(winding) || std::abs(winding) > 0.5 + 1e-6)
        return fail("captured source centroid is no longer outside by winding query");
    if (spin::quadricSolidInteriorIsPositiveRegion(quadric, triangulation))
        return fail("oriented interface failed to select the captured q-negative interior");
    std::cout << "[PASS] outward interface sign overrides the outside volume centroid\n";
    return EXIT_SUCCESS;
}

int sourceInteriorOffsetWitnessRegression()
{
    // Captured first deterministic seed-1 trial. The centroid is outside and
    // a small set of interface facets disagree because the finite triangulated
    // source is not perfectly aligned with the analytic quadric. A source-side
    // decision must be supported by points independently certified inside the
    // closed source, not by majority voting over boundary normals.
    Eigen::Matrix<double, 10, 1> values;
    values << -0.009715901854153618, 0.05490310208246309,
        0.005391325863425369, 0.10629342765972563,
        -0.02740799088056872, 0.32949278935694676,
        -0.9169842633806424, -0.014101472989526355,
        -0.164639263940708, 0.08990028823692814;
    const auto quadric = spin::QuadricCoefficients10::fromUnnormalized(values);
    const auto input = request();
    spin::QuadricTriangulationSettings settings;
    settings.icosphereSubdivisions = 5;
    settings.truncationRadius = 6.928203230275509;
    settings.maximumRadiusDoublings = 4;
    const auto triangulation = spin::triangulateQuadric(
        quadric, input.designDomain.vertices, settings);

    Eigen::Vector3d lower = triangulation.closedBoundary.vertices.colwise().minCoeff();
    Eigen::Vector3d upper = triangulation.closedBoundary.vertices.colwise().maxCoeff();
    const double length = (upper - lower).norm();
    int positiveInteriorSamples = 0;
    int negativeInteriorSamples = 0;
    double minimumCertifiedWinding = std::numeric_limits<double>::infinity();
    double maximumCertifiedWinding = -std::numeric_limits<double>::infinity();
    std::vector<Eigen::Index> interfaceFaces;
    for (Eigen::Index face = 0; face < triangulation.closedBoundary.faces.rows(); ++face)
        if (triangulation.faceIsCap(face) == 0) interfaceFaces.push_back(face);
    const int sampleCount = std::min<int>(32, interfaceFaces.size());
    for (int sample = 0; sample < sampleCount; ++sample)
    {
        const Eigen::Index face = interfaceFaces[
            static_cast<std::size_t>(sample) * (interfaceFaces.size() - 1) /
            std::max(1, sampleCount - 1)];
        const Eigen::Vector3i ids = triangulation.closedBoundary.faces.row(face);
        const Eigen::Vector3d a = triangulation.closedBoundary.vertices.row(ids[0]);
        const Eigen::Vector3d b = triangulation.closedBoundary.vertices.row(ids[1]);
        const Eigen::Vector3d c = triangulation.closedBoundary.vertices.row(ids[2]);
        Eigen::Vector3d normal = (b - a).cross(c - a);
        if (normal.norm() <= 1e-15) continue;
        normal.normalize();
        const Eigen::Vector3d surfacePoint = (a + b + c) / 3.0;
        for (int exponent = 6; exponent <= 14; ++exponent)
        {
            const Eigen::Vector3d point = surfacePoint -
                std::ldexp(length, -exponent) * normal;
            const double winding = igl::winding_number(
                triangulation.closedBoundary.vertices,
                triangulation.closedBoundary.faces,
                point);
            if (!std::isfinite(winding) || std::abs(winding) <= 0.5 + 1e-6)
                continue;
            const double q = quadric.evaluate(point);
            if (!std::isfinite(q) || std::abs(q) <= 1e-9) continue;
            if (q > 0.0) ++positiveInteriorSamples;
            else ++negativeInteriorSamples;
            minimumCertifiedWinding = std::min(minimumCertifiedWinding, winding);
            maximumCertifiedWinding = std::max(maximumCertifiedWinding, winding);
            break;
        }
    }
    if (positiveInteriorSamples + negativeInteriorSamples == 0)
        return fail("captured source has no certified interior interface-offset witness");

    int positiveVotes = 0;
    int negativeVotes = 0;
    int positiveMinorityInteriorSamples = 0;
    int negativeMinorityInteriorSamples = 0;
    double positiveArea = 0.0;
    double negativeArea = 0.0;
    double maximumVertexResidual = 0.0;
    for (Eigen::Index face = 0; face < triangulation.closedBoundary.faces.rows(); ++face)
    {
        if (triangulation.faceIsCap(face) != 0) continue;
        const Eigen::Vector3i ids = triangulation.closedBoundary.faces.row(face);
        const Eigen::Vector3d a = triangulation.closedBoundary.vertices.row(ids[0]);
        const Eigen::Vector3d b = triangulation.closedBoundary.vertices.row(ids[1]);
        const Eigen::Vector3d c = triangulation.closedBoundary.vertices.row(ids[2]);
        const Eigen::Vector3d area = (b - a).cross(c - a);
        const double areaLength = area.norm() * 0.5;
        const Eigen::Vector3d gradient = quadric.gradient((a + b + c) / 3.0);
        const double alignment = area.dot(gradient) /
            std::max(1e-300, area.norm() * gradient.norm());
        if (alignment > 0.0)
        {
            ++positiveVotes;
            positiveArea += areaLength;
            Eigen::Vector3d normal = area.normalized();
            const Eigen::Vector3d point = (a + b + c) / 3.0;
            for (int exponent = 6; exponent <= 14; ++exponent)
            {
                const Eigen::Vector3d inside = point - std::ldexp(length, -exponent) * normal;
                const double winding = igl::winding_number(
                    triangulation.closedBoundary.vertices,
                    triangulation.closedBoundary.faces,
                    inside);
                const double q = quadric.evaluate(inside);
                if (std::isfinite(winding) && std::abs(winding) > 0.5 + 1e-6 &&
                    std::isfinite(q) && std::abs(q) > 1e-9)
                {
                    if (q > 0.0) ++positiveMinorityInteriorSamples;
                    else ++negativeMinorityInteriorSamples;
                    break;
                }
            }
        }
        else { ++negativeVotes; negativeArea += areaLength; }
        maximumVertexResidual = std::max(maximumVertexResidual,
            std::max({std::abs(quadric.evaluate(a)), std::abs(quadric.evaluate(b)),
                      std::abs(quadric.evaluate(c))}));
    }
    std::cout << "[DIAGNOSTIC] positive_interior_samples=" << positiveInteriorSamples
              << " negative_interior_samples=" << negativeInteriorSamples
              << " winding_range=[" << minimumCertifiedWinding << ','
              << maximumCertifiedWinding << ']'
              << " outward_votes_positive=" << positiveVotes
              << " outward_votes_negative=" << negativeVotes
              << " minority_interior_positive=" << positiveMinorityInteriorSamples
              << " minority_interior_negative=" << negativeMinorityInteriorSamples
              << " minority_area_fraction="
              << positiveArea / std::max(1e-300, positiveArea + negativeArea)
              << " max_interface_vertex_abs_q=" << maximumVertexResidual << '\n';

    bool selectedPositive = false;
    try
    {
        selectedPositive = spin::quadricSolidInteriorIsPositiveRegion(
            quadric, triangulation);
    }
    catch (const std::exception& error)
    {
        if (positiveInteriorSamples == 0 || negativeMinorityInteriorSamples == 0)
            return fail("captured source did not establish both quadric sides inside its components: " +
                std::string(error.what()));
        std::cout << "[PASS] source with certified interior points on both quadric sides is rejected: "
                  << error.what() << '\n';
        return EXIT_SUCCESS;
    }
    (void)selectedPositive;
    return fail("source-side selector guessed despite certified interior points on both quadric sides");
}

spin::SpinDesignRequest request()
{
    spin::SpinDesignRequest result;
    result.outerSurface = cube(-2.0, 2.0);
    result.designDomain = cube(-1.0, 1.0);
    result.materialDensity = 1.0;
    result.minimumWallThickness = 0.5;
    return result;
}


spin::QuadricCoefficients10 sphere()
{
    spin::QuadricCoefficients10::Vector p =
        spin::QuadricCoefficients10::Vector::Zero();
    p(0) = 0.25;
    p(4) = -1.0;
    p(5) = -1.0;
    p(6) = -1.0;
    return spin::QuadricCoefficients10::fromUnnormalized(p);
}

int algebraicDiscreteEllipsoidRegression()
{
    // Isolate the source-anchor contract from the independent geometry guard.
    // This ill-conditioned source is NOT admitted as a material fixture:
    // its actual Boolean was rejected for degenerate faces (preserved log).
    spin::QuadricCoefficients10::Vector values;
    values<<-.125,0,0,-1,-1,-2,-1e-6,0,0,0;
    const auto q=spin::QuadricCoefficients10::fromUnnormalized(values);
    spin::QuadricMomentEvaluatorSettings settings;
    settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
    settings.triangulation.icosphereSubdivisions=1;
    settings.triangulation.algebraicFiniteMap=true;
    settings.triangulation.retainFixedNormalProvenance=true;
    const auto design=cube(-2e6,2e6);
    const auto source=spin::triangulateQuadric(q,design.vertices,settings.triangulation);
    spin::TriangleMesh empty;
    empty.vertices.resize(0,3);empty.faces.resize(0,3);
    const auto anchor=spin::differentiateDiscreteQuadricMoments(
        design,source,q,empty,Eigen::VectorXi(0),1.0,settings.triangulation);
    if(!anchor.defined)
        return fail("empty-output source-anchor contract: "+anchor.diagnostic);
    if(anchor.jacobian.norm()!=0.0)
        return fail("empty-output differential must be zero");

    // Separate well-conditioned, nonempty actual Boolean integration test.
    values<<.36,.04,-.03,.02,-1,-1.4,-1.8,.01,-.02,.03;
    const auto ordinary=spin::QuadricCoefficients10::fromUnnormalized(values);
    const auto result=spin::evaluateQuadricMaterialMoments(request(),ordinary,settings);
    std::cout<<"algebraic_discrete_ellipsoid: "<<result.derivativeDiagnostic<<'\n';
    if(!result.derivativeDefined)
        return fail("algebraic ellipsoid material derivative still uses legacy point: "+result.derivativeDiagnostic);
    if(!result.momentJacobian.allFinite() || result.materialBoundaryTarget.faces.rows()==0)
        return fail("algebraic ellipsoid lacks actual material or finite derivative");
    std::cout<<"[PASS] actual algebraic ellipsoid Boolean/discrete moment path\n";
    return EXIT_SUCCESS;
}

int projectiveRequestIdentityRegression()
{
    const auto input=request();
    for(int variant=0;variant<3;++variant)
    {
        spin::QuadricMomentEvaluatorSettings settings;
        settings.triangulation.maximumRadiusDoublings=8; // Historical fast-path tuple.
        if(variant==1)settings.triangulation.truncationRadius=16.;
        if(variant==2)settings.derivativeMethod=
            spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
        const auto legacy=spin::computeQuadricEvaluatorRequestHash(input,settings);
        if(variant==0 && legacy!=spin::computeSpinDesignRequestHash(input))
            return fail("legacy default request identity changed");
        settings.triangulation.paperProjectiveHyperboloid=true;
        const auto paper=spin::computeQuadricEvaluatorRequestHash(input,settings);
        if(legacy==paper)
            return fail("paper projective source map omitted from request identity, variant="+
                        std::to_string(variant));
        if(paper!=spin::computeQuadricEvaluatorRequestHash(input,settings))
            return fail("paper request identity is nondeterministic");
        settings.triangulation.paperProjectiveHyperboloid=false;
        if(legacy!=spin::computeQuadricEvaluatorRequestHash(input,settings))
            return fail("explicit legacy mode does not restore identity");
        settings.triangulation.algebraicFiniteMap=true;
        if(legacy==spin::computeQuadricEvaluatorRequestHash(input,settings))
            return fail("algebraic finite map omitted from evaluator request identity");
        settings.triangulation.algebraicFiniteMap=false;
        if(legacy!=spin::computeQuadricEvaluatorRequestHash(input,settings))
            return fail("algebraic policy changed legacy identity when disabled");
    }
    std::cout<<"[PASS] projective map identity distinguishes default and configured requests\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
#if MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
    if(argc==2 && std::string(argv[1])=="--rational-hyperboloid") {
        auto input=request();
        spin::QuadricMomentEvaluatorSettings settings;
        settings.derivativeMethod=spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
        settings.triangulation.icosphereSubdivisions=1;
        settings.triangulation.truncationRadius=8;
        settings.triangulation.rationalSourceCoordinates=true;
        const auto started=std::chrono::steady_clock::now();
        settings.diagnosticStageObserver=[&](const std::string& name) {
            std::cerr<<"[STAGE] "<<name<<" elapsed="<<std::chrono::duration<double>(
                std::chrono::steady_clock::now()-started).count()<<std::endl;
        };
        spin::QuadricCoefficients10::Vector p;
        p<<-.2,.07,-.04,.03,1.2,.9,-1.1,.02,.03,-.04;
        const auto result=spin::evaluateCanonicalFullRankMaterialMoments(input,
            spin::QuadricCoefficients10::fromUnnormalized(p),settings);
        if(!result.coefficientJacobian || result.variableMomentsExactTarget.size()!=10 ||
           result.materialExactVerticesTarget.empty() || result.sourceExactVerticesTarget.empty())
            return fail("actual rational hyperboloid Boolean/transport unavailable: "+result.diagnostic);
        std::cout<<"[PASS] actual rational hyperboloid Boolean, moments and same-source transport\n";
        return EXIT_SUCCESS;
    }
    if (argc == 2 && std::string(argv[1]) == "--rational-source")
    {
        auto input = request();
        spin::QuadricMomentEvaluatorSettings settings;
        settings.derivativeMethod = spin::QuadricMomentDerivativeMethod::DiscreteFullRankQuadric;
        settings.triangulation.icosphereSubdivisions = 1;
        const auto legacyRequestHash = spin::computeQuadricEvaluatorRequestHash(input, settings);
        settings.triangulation.rationalSourceCoordinates = true;
        const auto uninstrumentedRequestHash=spin::computeQuadricEvaluatorRequestHash(input,settings);
        std::vector<std::string> stages;
        settings.diagnosticStageObserver=[&](const std::string& name){stages.push_back(name);};
        const auto q = sphere();
        const auto result = spin::evaluateQuadricMaterialMoments(input, q, settings);
        if (spin::computeQuadricEvaluatorRequestHash(input, settings) == legacyRequestHash ||
            result.quadricTriangulation.exactVertices.empty() ||
            result.materialExactVerticesTarget.size() != static_cast<std::size_t>(result.materialBoundaryTarget.vertices.rows()) ||
            result.variableMomentsExactTarget.size() != 10 || !result.derivativeDefined)
            return fail("rational profile failed same-operation values, actual discrete derivative, or request identity: " + result.derivativeDiagnostic);
        if (!verification::auditExactTriangleBoundaryRational(result.materialExactVerticesTarget,
                result.materialBoundaryTarget.faces,
                verification::ExactBoundaryOrientation::OutwardPositive).facts().boundaryValid)
            return fail("actual rational Boolean material is invalid");
        if(stages.empty() || stages.front()!="source-construction-begin" ||
           stages.back()!="material-transport-end" ||
           spin::computeQuadricEvaluatorRequestHash(input,settings)!=uninstrumentedRequestHash)
            return fail("diagnostic stages are incomplete or changed scientific request identity");
        if(std::count(stages.begin(),stages.end(),"transport-material-authority-reused")!=1 ||
           std::count(stages.begin(),stages.end(),"transport-material-authority-audit-begin")!=0)
            return fail("production transport rebuilt the already authenticated SAME material boundary");
        settings.diagnosticStageObserver={};
        const auto full = spin::evaluateCanonicalFullRankMaterialMoments(input, q, settings);
        if (full.materialExactVerticesTarget != result.materialExactVerticesTarget ||
            full.variableMomentsExactTarget != result.variableMomentsExactTarget || !full.coefficientJacobian ||
            full.geometryHash!=result.geometryHash ||
            !full.coefficientJacobian->isApprox(result.momentJacobian,0))
            return fail("full-rank adapter lost same-operation authority coordinates or exact moments");
        std::cout << "[MAPPING] " << full.scientificHash << '\n';
        const auto motion=spin::evaluateQuadricSourceMotion(q,input.designDomain.vertices,settings.triangulation);
        if(!motion.defined)return fail("rational source motion cannot be reproduced");
        auto wrongSource=result.quadricTriangulation.exactVertices;
        wrongSource[0][0]="999/1";
        const auto invalidSource=spin::differentiatePolyhedralMaterialMoments(input.designDomain,
            result.quadricTriangulation.closedBoundary,motion.velocity,result.materialBoundaryTarget,
            result.materialBirthFaces,input.materialDensity,false,false,false,
            result.materialExactVerticesTarget,wrongSource,result.designExactVerticesTarget);
        if(invalidSource.defined)return fail("wrong authority source tape was accepted");
        auto wrongBirth=result.materialBirthFaces;wrongBirth(0)=0;
        const auto invalidBirth=spin::differentiatePolyhedralMaterialMoments(input.designDomain,
            result.quadricTriangulation.closedBoundary,motion.velocity,result.materialBoundaryTarget,
            wrongBirth,input.materialDensity,false,false,false,
            result.materialExactVerticesTarget,result.quadricTriangulation.exactVertices,result.designExactVerticesTarget);
        if(invalidBirth.defined)return fail("wrong Boolean source provenance was accepted");
        input.materialDensity = 2;
        const auto doubled = spin::evaluateQuadricMaterialMoments(input, q, settings);
        for (int i = 0; i < 10; ++i)
        {
            const auto original = verification::parseExactRationalCoordinate(result.variableMomentsExactTarget[i]);
            const auto twice = verification::parseExactRationalCoordinate(doubled.variableMomentsExactTarget[i]);
            if (twice != 2 * original || std::abs(CGAL::to_double(original) - result.variableMomentsTarget.values()(i)) > 1e-14)
                return fail("same-operation exact moments lost density or value consistency");
        }
        // Use actual intersecting facets, not just a source fully within the cube.
        auto p = q.values(); p(0) *= 6;
        const auto cut = spin::evaluateQuadricMaterialMoments(input,
            spin::QuadricCoefficients10::fromUnnormalized(p), settings);
        if (!cut.derivativeDefined || cut.constructedIntersectionVertexCount == 0)
            return fail("rational Boolean intersection transport failed: " + cut.derivativeDiagnostic);
        const double h = 1e-5;
        auto plus = q.values(), minus = q.values(); plus(1) += h; minus(1) -= h;
        const auto a = spin::evaluateQuadricMaterialMoments(request(), spin::QuadricCoefficients10::fromUnnormalized(plus), settings);
        const auto b = spin::evaluateQuadricMaterialMoments(request(), spin::QuadricCoefficients10::fromUnnormalized(minus), settings);
        const auto chord = (spin::QuadricCoefficients10::fromUnnormalized(plus).values() -
            spin::QuadricCoefficients10::fromUnnormalized(minus).values()) / (2*h);
        for (int i = 0; i < 10; ++i)
        {
            const auto exact = (verification::parseExactRationalCoordinate(a.variableMomentsExactTarget[i]) -
                verification::parseExactRationalCoordinate(b.variableMomentsExactTarget[i])) / CGAL::Gmpq(2*h);
            const double expected = (result.momentJacobian * chord)(i);
            if (std::abs(CGAL::to_double(exact) - expected) / std::max(1.0, std::abs(expected)) > 1e-6)
                return fail("rational source/Boolean/value/derivative local secant is inconsistent");
        }
        std::cout << "[PASS] actual rational source Boolean, exact density moments, source transport and full-rank payload\n";
        return EXIT_SUCCESS;
    }
#endif
    if (argc > 1 && std::string(argv[1]) == "--algebraic-discrete")
        return algebraicDiscreteEllipsoidRegression();
    if (argc > 1 && std::string(argv[1]) == "--paper-request-identity")
        return projectiveRequestIdentityRegression();
    if (argc > 1 && std::string(argv[1]) == "--certified-empty")
        return certifiedEmptyRegression();
    if (argc > 1 && std::string(argv[1]) == "--degenerate-boolean-witness")
        return degenerateBooleanWitnessRegression();
    if (argc > 1 && std::string(argv[1]) == "--source-winding-regression")
        return sourceWindingRegression();
    if (argc > 1 && std::string(argv[1]) == "--source-interior-side-witness")
        return sourceInteriorSideWitnessRegression();
    if (argc > 1 && std::string(argv[1]) == "--source-interior-offset-witness")
        return sourceInteriorOffsetWitnessRegression();
    if(projectiveRequestIdentityRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    if(algebraicDiscreteEllipsoidRegression()!=EXIT_SUCCESS)return EXIT_FAILURE;
    const spin::SpinDesignRequest input = request();
    spin::QuadricMomentEvaluatorSettings settings;
    settings.triangulation.icosphereSubdivisions = 3;

    const std::string defaultRequestHash =
        spin::computeQuadricEvaluatorRequestHash(input, settings);
    for (const auto& changed : {
             std::pair<const char*, spin::QuadricMomentEvaluatorSettings>{
                 "hyperbolic", [&] { auto v = settings; v.triangulation.hyperbolicCylinderGridResolution = 25; return v; }()},
                 std::pair<const char*, spin::QuadricMomentEvaluatorSettings>{
                 "parabolic", [&] { auto v = settings; v.triangulation.parabolicCylinderGridResolution = 25; return v; }()}})
    {
        if (spin::computeQuadricEvaluatorRequestHash(input, changed.second) == defaultRequestHash)
            return fail(std::string(changed.first) + " cylinder resolution was omitted from evaluator request identity");
    }
    for (const auto& invalid : {
             std::pair<const char*, spin::QuadricMomentEvaluatorSettings>{
                 "hyperbolic", [&] { auto v = settings; v.triangulation.hyperbolicCylinderGridResolution = 1; return v; }()},
                 std::pair<const char*, spin::QuadricMomentEvaluatorSettings>{
                 "parabolic", [&] { auto v = settings; v.triangulation.parabolicCylinderGridResolution = 0; return v; }()}})
    {
        try
        {
            (void)spin::computeQuadricEvaluatorRequestHash(input, invalid.second);
            return fail(std::string("invalid ") + invalid.first + " cylinder resolution was accepted");
        }
        catch (const std::invalid_argument&)
        {
        }
    }

    const spin::QuadricMomentEvaluation inside =
        spin::evaluateQuadricMaterialMoments(input, sphere(), settings);
    const std::string expectedBackendIdentity =
        spin::quadricBooleanBackendIdentity();
    if (expectedBackendIdentity.empty() ||
        inside.booleanBackendIdentity != expectedBackendIdentity)
    {
        return fail("material evaluation does not bind its actual Boolean backend identity");
    }
    const manufacturing::BoundaryAudit insideAudit =
        manufacturing::auditTriangleBoundary(
            inside.materialBoundary.vertices,
            inside.materialBoundary.faces,
            manufacturing::BoundaryOrientation::AnyNonzero
        );
    if (
        !inside.derivativeDefined || !insideAudit.strictBoundaryValid() ||
        inside.variableMomentsTarget.mass() <= 0.0 ||
        inside.interfaceTriangleCount <= 0 ||
        inside.momentJacobian.rows() != 10 ||
        inside.momentJacobian.cols() != 10 ||
        !inside.momentJacobian.allFinite() ||
        (inside.momentJacobian * sphere().values()).norm() > 2e-10 ||
        inside.geometryHash.size() != 64 || inside.topologyHash.size() != 64
    )
    {
        return fail("contained sphere material artifact is invalid");
    }

    const spin::QuadricMomentEvaluation outside =
        spin::evaluateQuadricMaterialMoments(
            input,
            sphere().inverted(),
            settings
        );
    const double volumeSum = inside.variableMomentsTarget.mass() +
        outside.variableMomentsTarget.mass();
    if (
        std::abs(volumeSum - 8.0) > 2e-10 ||
        outside.variableMomentsTarget.mass() <=
            inside.variableMomentsTarget.mass()
    )
    {
        return fail(
            "quadric sign inversion did not complement the domain: inside=" +
            std::to_string(inside.variableMomentsTarget.mass()) +
            " outside=" +
            std::to_string(outside.variableMomentsTarget.mass()) +
            " sum=" + std::to_string(volumeSum)
        );
    }

    spin::QuadricCoefficients10::Vector planeValues =
        spin::QuadricCoefficients10::Vector::Zero();
    planeValues(1) = 1.0;
    const spin::QuadricMomentEvaluation half =
        spin::evaluateQuadricMaterialMoments(
            input,
            spin::QuadricCoefficients10::fromUnnormalized(planeValues),
            settings
        );
    if (
        std::abs(half.variableMomentsTarget.mass() - 4.0) > 1e-12 ||
        std::abs(half.variableMomentsTarget.centerOfMass().x() - 0.5) >
            1e-12
    )
    {
        return fail("planar positive half-space moments are incorrect");
    }

    spin::QuadricMomentEvaluatorSettings reducedPlaneSettings = settings;
    reducedPlaneSettings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum;
    reducedPlaneSettings.triangulation.retainFixedNormalProvenance = true;
    spin::QuadricCoefficients10::Vector genericPlaneValues =
        spin::QuadricCoefficients10::Vector::Zero();
    genericPlaneValues << 0.14, 0.31, -0.47, 0.826, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    const auto genericPlane =
        spin::QuadricCoefficients10::fromUnnormalized(genericPlaneValues);
    const spin::QuadricMomentEvaluation genericPlaneMoments =
        spin::evaluateQuadricMaterialMoments(input, genericPlane, settings);
    const spin::QuadricMomentEvaluation reducedPlane =
        spin::evaluateQuadricMaterialMoments(
            input,
            genericPlane,
            reducedPlaneSettings
        );
    if (
        reducedPlane.derivativeMethod !=
            spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum ||
        !reducedPlane.stratumDifferential.defined ||
        reducedPlane.stratumDifferential.jacobian.rows() != 10 ||
        reducedPlane.stratumDifferential.jacobian.cols() != 3 ||
        reducedPlane.stratumDifferential.coefficientBasis.rows() != 10 ||
        reducedPlane.stratumDifferential.coefficientBasis.cols() != 3 ||
        reducedPlane.derivativeDefined ||
        reducedPlane.momentJacobian.norm() != 0.0 ||
        !reducedPlane.stratumDifferential.jacobian.allFinite() ||
        reducedPlane.geometryHash.size() != 64
    )
    {
        return fail("unified evaluator did not preserve the reduced plane differential contract: " +
            reducedPlane.stratumDifferential.diagnostic +
            " basis=" + std::to_string(reducedPlane.stratumDifferential.coefficientBasis.rows()) +
            "x" + std::to_string(reducedPlane.stratumDifferential.coefficientBasis.cols()) +
            " jacobian=" + std::to_string(reducedPlane.stratumDifferential.jacobian.rows()) +
            "x" + std::to_string(reducedPlane.stratumDifferential.jacobian.cols()) +
            " ambient_defined=" + std::to_string(reducedPlane.derivativeDefined));
    }
    if (reducedPlane.variableMomentsTarget.values().isApprox(
            genericPlaneMoments.variableMomentsTarget.values(), 0.0) == false ||
        spin::computeQuadricEvaluatorRequestHash(input, reducedPlaneSettings) ==
            spin::computeQuadricEvaluatorRequestHash(input, settings))
    {
        return fail("reduced plane evaluation changed moments or reused the legacy request identity");
    }
    double previousReducedError = std::numeric_limits<double>::infinity();
    for (double step : {1e-4, 3e-5, 1e-5})
    {
        Eigen::Matrix<double, 10, 3> finiteDifference;
        for (int column = 0; column < 3; ++column)
        {
            const Eigen::Matrix<double, 10, 1> direction =
                reducedPlane.stratumDifferential.coefficientBasis.col(column);
            const auto plus = spin::QuadricCoefficients10::fromUnnormalized(
                genericPlane.values() + step * direction);
            const auto minus = spin::QuadricCoefficients10::fromUnnormalized(
                genericPlane.values() - step * direction);
            finiteDifference.col(column) =
                (spin::evaluateQuadricMaterialMoments(input, plus, settings)
                    .variableMomentsTarget.values() -
                 spin::evaluateQuadricMaterialMoments(input, minus, settings)
                    .variableMomentsTarget.values()) / (2.0 * step);
        }
        const double error =
            (reducedPlane.stratumDifferential.jacobian - finiteDifference).norm() /
            std::max(1.0, finiteDifference.norm());
        if (!std::isfinite(error) || error > previousReducedError * 1.05 ||
            error > 1e-6)
            return fail("unified reduced-plane Jacobian disagrees with fixed-map finite differences: " +
                std::to_string(error));
        previousReducedError = error;
    }

    const auto verifyReducedStratum = [&](const spin::SpinDesignRequest& stratumRequest,
                                          const spin::QuadricCoefficients10& stratumQuadric,
                                          spin::QuadricMomentEvaluatorSettings stratumSettings,
                                          Eigen::Index expectedDimension,
                                          const std::string& label) {
        stratumSettings.derivativeMethod =
            spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum;
        stratumSettings.triangulation.retainFixedNormalProvenance = true;
        const auto evaluation = spin::evaluateQuadricMaterialMoments(
            stratumRequest, stratumQuadric, stratumSettings);
        const auto independent = spin::differentiateDiscreteQuadricStratumMoments(
            stratumRequest.designDomain, evaluation.quadricTriangulation,
            stratumQuadric, evaluation.materialBoundaryTarget,
            evaluation.materialBirthFaces, stratumRequest.materialDensity,
            stratumSettings.triangulation);
        if (!evaluation.stratumDifferential.defined || !independent.defined ||
            evaluation.derivativeDefined ||
            evaluation.stratumDifferential.jacobian.rows() != 10 ||
            evaluation.stratumDifferential.jacobian.cols() != expectedDimension ||
            evaluation.stratumDifferential.coefficientBasis.rows() != 10 ||
            evaluation.stratumDifferential.coefficientBasis.cols() != expectedDimension ||
            !evaluation.stratumDifferential.jacobian.isApprox(independent.jacobian, 0.0) ||
            !evaluation.stratumDifferential.coefficientBasis.isApprox(
                independent.coefficientBasis, 0.0) ||
            evaluation.provenanceTopologyHash.size() != 64 ||
            evaluation.derivativeScientificHash.size() != 64)
        {
            std::cerr << "[FAIL] unified evaluator lost " << label
                      << " reduced differential: "
                      << evaluation.stratumDifferential.diagnostic << '\n';
            return false;
        }
        std::cout << "[PASS] unified evaluator " << label
                  << " reduced Jacobian 10x" << expectedDimension << '\n';
        return true;
    };

    spin::QuadricCoefficients10::Vector cylinderValues =
        spin::QuadricCoefficients10::Vector::Zero();
    cylinderValues << -1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0;
    spin::QuadricMomentEvaluatorSettings cylinderSettings;
    cylinderSettings.triangulation.cylinderAngularResolution = 24;
    cylinderSettings.triangulation.truncationRadius = 8.0;
    if (!verifyReducedStratum(input,
            spin::QuadricCoefficients10::fromUnnormalized(cylinderValues),
            cylinderSettings, 7, "elliptic-cylinder"))
        return EXIT_FAILURE;

    spin::SpinDesignRequest hyperbolicCylinderRequest;
    hyperbolicCylinderRequest.outerSurface = cube(-2.0, 2.0);
    hyperbolicCylinderRequest.designDomain = cube(-0.8, 0.8);
    const Eigen::Vector3d hyperbolicDomainOffset(0.113, -0.071, 0.057);
    hyperbolicCylinderRequest.designDomain.vertices.rowwise() +=
        hyperbolicDomainOffset.transpose();
    hyperbolicCylinderRequest.materialDensity = 1.0;
    spin::QuadricCoefficients10::Vector hyperbolicCylinderValues =
        spin::QuadricCoefficients10::Vector::Zero();
    hyperbolicCylinderValues << -0.2, 0.0, 0.0, 0.0, 1.0, -0.7, 0.0,
        0.0, 0.0, 0.0;
    spin::QuadricMomentEvaluatorSettings hyperbolicCylinderSettings;
    hyperbolicCylinderSettings.triangulation.hyperbolicCylinderGridResolution = 16;
    hyperbolicCylinderSettings.triangulation.truncationRadius = 8.0;
    const auto hyperbolicCylinderQuadric =
        spin::QuadricCoefficients10::fromUnnormalized(hyperbolicCylinderValues);
    if (!verifyReducedStratum(hyperbolicCylinderRequest,
            hyperbolicCylinderQuadric, hyperbolicCylinderSettings, 7,
            "hyperbolic-cylinder"))
        return EXIT_FAILURE;
    hyperbolicCylinderSettings.derivativeMethod =
        spin::QuadricMomentDerivativeMethod::DiscreteReducedStratum;
    hyperbolicCylinderSettings.triangulation.retainFixedNormalProvenance = true;
    const auto hyperbolicCylinderEvaluation =
        spin::evaluateQuadricMaterialMoments(hyperbolicCylinderRequest,
            hyperbolicCylinderQuadric, hyperbolicCylinderSettings);
    double previousHyperbolicCylinderError =
        std::numeric_limits<double>::infinity();
    for (double step : {1e-4, 3e-5, 1e-5})
    {
        Eigen::MatrixXd finiteDifference(10, 7);
        for (int column = 0; column < 7; ++column)
        {
            const Eigen::Matrix<double, 10, 1> direction =
                hyperbolicCylinderEvaluation.stratumDifferential.coefficientBasis.col(column);
            const auto plus = hyperbolicCylinderRetraction(
                hyperbolicCylinderQuadric, direction, step);
            const auto minus = hyperbolicCylinderRetraction(
                hyperbolicCylinderQuadric, direction, -step);
            try
            {
                finiteDifference.col(column) =
                    (spin::evaluateQuadricMaterialMoments(hyperbolicCylinderRequest,
                        plus, hyperbolicCylinderSettings).variableMomentsTarget.values() -
                     spin::evaluateQuadricMaterialMoments(hyperbolicCylinderRequest,
                        minus, hyperbolicCylinderSettings).variableMomentsTarget.values()) /
                    (2.0 * step);
            }
            catch (const std::exception& error)
            {
                return fail("hyperbolic-cylinder finite-difference sample failed at h=" +
                    std::to_string(step) + " column=" + std::to_string(column) +
                    ": " + error.what());
            }
        }
        const double error =
            (hyperbolicCylinderEvaluation.stratumDifferential.jacobian -
             finiteDifference).norm() / std::max(1.0, finiteDifference.norm());
        std::cout << "[INFO] hyperbolic-cylinder stratum FD h=" << step
                  << " relative_error=" << error << '\n';
        if (!std::isfinite(error) || error > previousHyperbolicCylinderError * 1.05 ||
            (step == 1e-5 && error > 1e-6))
            return fail("unified hyperbolic-cylinder reduced Jacobian disagrees with "
                "multi-step Boolean finite differences: " + std::to_string(error));
        previousHyperbolicCylinderError = error;
    }

    spin::SpinDesignRequest paraboloidRequest;
    paraboloidRequest.outerSurface = cube(-2.0, 2.0);
    paraboloidRequest.designDomain = cube(-0.8, 0.8);
    const Eigen::Vector3d domainOffset(0.113, -0.071, 0.057);
    paraboloidRequest.designDomain.vertices.rowwise() += domainOffset.transpose();
    paraboloidRequest.materialDensity = 1.0;
    spin::QuadricCoefficients10::Vector paraboloidValues =
        spin::QuadricCoefficients10::Vector::Zero();
    paraboloidValues << -0.2, 0.0, 0.0, -1.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0;
    spin::QuadricMomentEvaluatorSettings paraboloidSettings;
    paraboloidSettings.triangulation.paraboloidGridResolution = 16;
    paraboloidSettings.triangulation.truncationRadius = 2.0;
    if (!verifyReducedStratum(paraboloidRequest,
            spin::QuadricCoefficients10::fromUnnormalized(paraboloidValues),
            paraboloidSettings, 8, "paraboloid"))
        return EXIT_FAILURE;

    try
    {
        (void)spin::solveQuadricHollowing(
            input, reducedPlaneSettings, {genericPlane});
        return fail("reduced stratum derivative was accepted by the 10-D Eq. 14 solver");
    }
    catch (const std::invalid_argument& error)
    {
        if (std::string(error.what()).find("evaluation-only") == std::string::npos)
            return fail("reduced-stratum solver rejection was not explicit: " +
                std::string(error.what()));
    }
    std::cout << "[PASS] raw ten-moment volume and first-moment convention\n";

    constexpr double step = 2e-5;
    Eigen::Matrix<double, 10, 10> finiteDifference;
    for (int coordinate = 0; coordinate < 10; ++coordinate)
    {
        spin::QuadricCoefficients10::Vector plus = sphere().values();
        spin::QuadricCoefficients10::Vector minus = sphere().values();
        plus(coordinate) += step;
        minus(coordinate) -= step;
        const Eigen::Matrix<double, 10, 1> plusMoments =
            spin::evaluateQuadricMaterialMoments(
                input,
                spin::QuadricCoefficients10::fromUnnormalized(plus),
                settings
            ).variableMomentsTarget.values();
        const Eigen::Matrix<double, 10, 1> minusMoments =
            spin::evaluateQuadricMaterialMoments(
                input,
                spin::QuadricCoefficients10::fromUnnormalized(minus),
                settings
            ).variableMomentsTarget.values();
        finiteDifference.col(coordinate) =
            (plusMoments - minusMoments) / (2.0 * step);
    }
    const double relativeDerivativeError =
        (inside.momentJacobian - finiteDifference).norm() /
        std::max(1.0, finiteDifference.norm());
    // Opt-in negative control preserves the historical fixture and legacy gate.
    const bool strictNegativeControl = argc > 1 &&
        std::string(argv[1]) == "--strict-discrete-negative-control";
    if (relativeDerivativeError > (strictNegativeControl ? 1e-6 : 0.08))
    {
        return fail(
            "shape Jacobian disagrees with symmetric finite differences: " +
            std::to_string(relativeDerivativeError)
        );
    }

    const std::string outerBefore =
        spin::computeSpinDesignRequestHash(input);
    const std::string outerAfter =
        spin::computeSpinDesignRequestHash(input);
    if (outerBefore != outerAfter ||
        inside.geometryHash != spin::evaluateQuadricMaterialMoments(
            input,
            sphere(),
            settings
        ).geometryHash)
    {
        return fail("evaluation mutated input or is non-deterministic");
    }

    std::cout << "[PASS] scenario=quadric-moment-evaluator derivative_error="
              << relativeDerivativeError << " geometry="
              << inside.geometryHash << '\n';
    return certifiedEmptyRegression();
}
