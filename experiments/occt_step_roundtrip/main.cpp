#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeHalfSpace.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Hyperbola.hxx>
#include <GProp_GProps.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Precision.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <StepData_StepModel.hxx>
#include <DESTEP_Parameters.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <gp_Pnt.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Circ.hxx>
#include <gp_Mat.hxx>
#include <gp_Pln.hxx>
#include <gp_Hypr.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>
#include <gmpxx.h>

#include <cmath>
#include <array>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Moments {
  double volume;
  double x;
  double y;
  double z;
  double raw[10];
};

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "occt_step_roundtrip: " << message << '\n';
  std::exit(1);
}

// Optional diagnostic output for the *rejected* original shape. This never
// heals geometry or relaxes the analyzer and is not a material/CAD approval.
std::string diagnostic_json_string(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (const unsigned char byte:value) {
    if (byte=='"' || byte=='\\') out << '\\' << static_cast<char>(byte);
    else if (byte<0x20)
      out << "\\u00" << std::hex << std::setw(2) << std::setfill('0')
          << static_cast<unsigned>(byte) << std::dec;
    else out << static_cast<char>(byte);
  }
  out << '"';
  return out.str();
}

void require_valid_topology(const TopoDS_Shape& shape, const std::string& stage,
                            const std::string& failure) {
  BRepCheck_Analyzer analyzer(shape, Standard_True);
  if (analyzer.IsValid()) return;
  const char* requested = std::getenv("MAKE_IT_STAND_OCCT_FAILURE_DIR");
  if (requested && *requested) {
    try {
      const std::filesystem::path directory(requested);
      std::error_code error;
      if (!std::filesystem::create_directory(directory,error)) {
        std::cerr << "topology witness: refusing to overwrite or create directory: "
                  << directory << ": " << error.message() << '\n';
      } else {
        if (!BRepTools::Write(shape,(directory/"invalid.brep").string().c_str()))
          throw std::runtime_error("cannot save rejected B-rep");
        TopTools_IndexedMapOfShape shapes;
        TopExp::MapShapes(shape,shapes);
        std::ofstream out(directory/"topology.json");
        out << std::setprecision(17)
            << "{\"schema\":\"occt-rejected-topology-v1\",\"valid\":false,"
               "\"stage\":" << diagnostic_json_string(stage)
            << ",\"occt_version\":" << diagnostic_json_string(OCC_VERSION_COMPLETE)
            << ",\"subshape_count\":" << shapes.Extent() << ",\"defects\":[";
        bool first = true;
        for (int i=1;i<=shapes.Extent();++i) {
          const TopoDS_Shape& subshape=shapes(i);
          const Handle(BRepCheck_Result)& result=analyzer.Result(subshape);
          if (result.IsNull()) continue;
          const auto append = [&](const BRepCheck_ListOfStatus& statuses,int context) {
            for (BRepCheck_ListOfStatus::Iterator it(statuses);it.More();it.Next()) {
              const BRepCheck_Status status=it.Value();
              if (status==BRepCheck_NoError) continue;
              std::ostringstream name;
              BRepCheck::Print(status,name);
              std::string label=name.str();
              while (!label.empty() && (label.back()=='\n' || label.back()=='\r'))
                label.pop_back();
              if (!first) out << ',';
              first=false;
              out << "{\"subshape_index\":" << i
                  << ",\"shape_type\":" << static_cast<int>(subshape.ShapeType())
                  << ",\"context_index\":" << context
                  << ",\"status_code\":" << static_cast<int>(status)
                  << ",\"status\":" << diagnostic_json_string(label) << '}';
            }
          };
          append(result->Status(),0);
          for (result->InitContextIterator();result->MoreShapeInContext();
               result->NextShapeInContext())
            append(result->StatusOnShape(),shapes.FindIndex(result->ContextualShape()));
        }
        out << "]}\n";
        out.close();
        if (!out) throw std::runtime_error("cannot save topology status report");
      }
    } catch (const Standard_Failure& error) {
      std::cerr << "topology witness failed: " << error.GetMessageString() << '\n';
    } catch (const std::exception& error) {
      std::cerr << "topology witness failed: " << error.what() << '\n';
    }
  }
  fail(failure);
}

void trace_carrier_topology(const TopoDS_Shape& shape,const char* stage) {
  const char* requested=std::getenv("MAKE_IT_STAND_OCCT_FAILURE_DIR");
  if (!requested || !*requested) return;
  TopTools_IndexedMapOfShape shapes;
  TopExp::MapShapes(shape,shapes);
  std::array<int,9> counts{};
  for (int i=1;i<=shapes.Extent();++i)
    ++counts.at(static_cast<unsigned>(shapes(i).ShapeType()));
  std::cerr << "OCCT_CARRIER_STAGE {\"stage\":" << diagnostic_json_string(stage)
            << ",\"valid\":"
            << (BRepCheck_Analyzer(shape,Standard_True).IsValid()?"true":"false")
            << ",\"subshape_count\":" << shapes.Extent() << ",\"type_counts\":[";
  for (std::size_t i=0;i<counts.size();++i) {
    if (i) std::cerr << ',';
    std::cerr << counts[i];
  }
  std::cerr << "]}\n";
}

Moments volume_moments(const TopoDS_Shape& shape) {
  if (shape.IsNull()) {
    fail("received a null shape");
  }
  GProp_GProps properties;
  const double integration_error = BRepGProp::VolumeProperties(
      shape, properties, 1e-12, Standard_True, Standard_False);
  if (!std::isfinite(integration_error) || integration_error < 0.0 ||
      integration_error > 1e-10) {
    std::ostringstream message;
    message << std::setprecision(17)
            << "OCCT adaptive volume integration did not meet its error contract (reported "
            << integration_error << ", required <= 1e-10)";
    fail(message.str());
  }
  const gp_Pnt center = properties.CentreOfMass();
  const double mass = properties.Mass();
  const gp_Mat inertia = properties.MatrixOfInertia();
  const double trace = inertia.Value(1, 1) + inertia.Value(2, 2) + inertia.Value(3, 3);
  // The inertia matrix is central. Recover the central second-product matrix
  // Q = integral (x-c)(x-c)^T dV, then translate it to the origin.
  double q[3][3] = {};
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      const double identity = row == column ? 1.0 : 0.0;
      q[row][column] = 0.5 * trace * identity - inertia.Value(row + 1, column + 1);
    }
  }
  const double c[3] = {center.X(), center.Y(), center.Z()};
  double raw_q[3][3] = {};
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      raw_q[row][column] = q[row][column] + mass * c[row] * c[column];
    }
  }
  Moments out{mass, center.X(), center.Y(), center.Z(), {}};
  out.raw[0] = mass;
  out.raw[1] = mass * center.X();
  out.raw[2] = mass * center.Y();
  out.raw[3] = mass * center.Z();
  out.raw[4] = raw_q[0][0];
  out.raw[5] = raw_q[1][1];
  out.raw[6] = raw_q[2][2];
  out.raw[7] = raw_q[1][2];
  out.raw[8] = raw_q[0][2];
  out.raw[9] = raw_q[0][1];
  return out;
}

void emit_moments(const Moments& value) {
  std::cout << "{\"volume\":" << value.volume << ",\"center\":[" << value.x
            << ',' << value.y << ',' << value.z << "]}";
}

void emit_raw_moments(const Moments& value) {
  std::cout << '[';
  for (int i = 0; i < 10; ++i) {
    if (i != 0) std::cout << ',';
    std::cout << value.raw[i];
  }
  std::cout << ']';
}

double parse_coordinate_scale(const std::string& text) {
  errno = 0;
  char* end = nullptr;
  const double scale = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE ||
      !std::isfinite(scale) || !(scale > 0.0)) {
    fail("coordinate scale must be a finite positive number");
  }
  return scale;
}

double parse_finite_double(const std::string& text, const std::string& label) {
  errno = 0;
  char* end = nullptr;
  const double value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || *end != '\0' || errno == ERANGE || !std::isfinite(value))
    fail(label + " must be a finite number");
  return value;
}

struct SphereQuadric {
  gp_Pnt center;
  double radius;
};

SphereQuadric sphere_from_solver_coefficients(const std::array<double, 10>& q) {
  double norm_squared = 0.0;
  for (const double value : q) norm_squared += value * value;
  const double norm = std::sqrt(norm_squared);
  if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1e-10)
    fail("solver quadric coefficients must have unit Euclidean norm");

  // Solver convention: c + b.x + x^T A x, with coefficient order
  // [c,bx,by,bz,Axx,Ayy,Azz,Ayz,Axz,Axy].
  const double lambda = (q[4] + q[5] + q[6]) / 3.0;
  const double coefficient_scale = std::max(1.0, std::abs(lambda));
  const double isotropic_tolerance = 1e-10 * coefficient_scale;
  if (!(lambda > 0.0) || std::abs(q[4] - lambda) > isotropic_tolerance ||
      std::abs(q[5] - lambda) > isotropic_tolerance ||
      std::abs(q[6] - lambda) > isotropic_tolerance ||
      std::abs(q[7]) > isotropic_tolerance || std::abs(q[8]) > isotropic_tolerance ||
      std::abs(q[9]) > isotropic_tolerance)
    fail("solver quadric is not an isotropic sphere quadric");

  const double cx = -q[1] / (2.0 * lambda);
  const double cy = -q[2] / (2.0 * lambda);
  const double cz = -q[3] / (2.0 * lambda);
  const double radius_squared = cx * cx + cy * cy + cz * cz - q[0] / lambda;
  if (!std::isfinite(radius_squared) || !(radius_squared > 0.0))
    fail("solver sphere quadric must have a positive finite radius");
  return {gp_Pnt(cx, cy, cz), std::sqrt(radius_squared)};
}

struct EllipsoidQuadric {
  gp_Pnt center;
  gp_GTrsf unit_sphere_map;
  double volume;
};

EllipsoidQuadric ellipsoid_from_solver_coefficients(const std::array<double, 10>& q) {
  double norm_squared = 0.0;
  for (const double value : q) norm_squared += value * value;
  const double norm = std::sqrt(norm_squared);
  if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1e-10)
    fail("solver quadric coefficients must have unit Euclidean norm");

  // Solver convention: c + b.x + x^T A x, order
  // [c,bx,by,bz,Axx,Ayy,Azz,Ayz,Axz,Axy].
  const double a[3][3] = {{q[4], q[9], q[8]},
                          {q[9], q[5], q[7]},
                          {q[8], q[7], q[6]}};
  const double b[3] = {q[1], q[2], q[3]};
  double lower[3][3] = {};
  const double matrix_scale = std::max({std::abs(a[0][0]), std::abs(a[1][1]),
                                        std::abs(a[2][2]), 1e-300});
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column <= row; ++column) {
      double value = a[row][column];
      for (int k = 0; k < column; ++k) value -= lower[row][k] * lower[column][k];
      if (row == column) {
        if (!std::isfinite(value) || value <= 1e-12 * matrix_scale)
          fail("solver quadric is not strictly positive definite");
        lower[row][column] = std::sqrt(value);
      } else {
        lower[row][column] = value / lower[column][column];
      }
    }
  }

  // Solve A*center=-b/2 through the Cholesky factors.
  double intermediate[3] = {}, center[3] = {};
  for (int row = 0; row < 3; ++row) {
    double value = -0.5 * b[row];
    for (int k = 0; k < row; ++k) value -= lower[row][k] * intermediate[k];
    intermediate[row] = value / lower[row][row];
  }
  for (int row = 2; row >= 0; --row) {
    double value = intermediate[row];
    for (int k = row + 1; k < 3; ++k) value -= lower[k][row] * center[k];
    center[row] = value / lower[row][row];
  }
  double center_level = -q[0];
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column)
      center_level += center[row] * a[row][column] * center[column];
  if (!std::isfinite(center_level) || !(center_level > 0.0))
    fail("solver ellipsoid quadric must enclose a positive finite volume");

  // Map the unit sphere by M=sqrt(center_level)*inverse(L^T), so
  // (M*y)^T A (M*y)=center_level for every unit vector y.
  double inverse_upper[3][3] = {};
  for (int column = 0; column < 3; ++column) {
    for (int row = 2; row >= 0; --row) {
      double value = row == column ? 1.0 : 0.0;
      for (int k = row + 1; k < 3; ++k)
        value -= lower[k][row] * inverse_upper[k][column];
      inverse_upper[row][column] = value / lower[row][row];
    }
  }
  const double scale = std::sqrt(center_level);
  double map[3][3] = {};
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column)
      map[row][column] = scale * inverse_upper[row][column];
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  gp_GTrsf transform(vectorial, gp_XYZ(center[0], center[1], center[2]));
  const double determinant = vectorial.Determinant();
  const double volume = 4.0 * std::acos(-1.0) * std::abs(determinant) / 3.0;
  if (!std::isfinite(volume) || !(volume > 0.0))
    fail("solver ellipsoid has nonpositive or non-finite volume");
  return {gp_Pnt(center[0], center[1], center[2]), transform, volume};
}

struct RankTwoQuadricFrame {
  gp_Pnt axis_point;
  gp_Dir axis;
  gp_Dir radial_x;
  gp_Dir radial_y;
  double lambda_x;
  double lambda_y;
  double center_level;
  double axial_linear;
};

RankTwoQuadricFrame rank_two_quadric_from_solver_coefficients(
    const std::array<double, 10>& q, bool requireSaddle = false) {
  double norm_squared = 0.0;
  for (const double value : q) norm_squared += value * value;
  const double norm = std::sqrt(norm_squared);
  if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1e-10)
    fail("solver quadric coefficients must have unit Euclidean norm");

  // Solver convention: c + b.x + x^T A x, order
  // [c,bx,by,bz,Axx,Ayy,Azz,Ayz,Axz,Axy].
  const double a[3][3] = {{q[4], q[9], q[8]},
                          {q[9], q[5], q[7]},
                          {q[8], q[7], q[6]}};
  const double b[3] = {q[1], q[2], q[3]};
  double matrix_scale = 0.0;
  for (const auto& row : a)
    for (const double value : row) matrix_scale = std::max(matrix_scale, std::abs(value));
  if (!(matrix_scale > 0.0)) fail("rank-two quadric has a zero quadratic block");

  // An elliptic cylinder has a rank-two positive-semidefinite quadratic
  // block. Recover the axis from the two independent rows, then diagonalize
  // the restriction to the perpendicular plane and verify the full matrix.
  std::array<double, 3> axis{};
  double best_cross_squared = 0.0;
  for (int first = 0; first < 3; ++first) {
    for (int second = first + 1; second < 3; ++second) {
      const std::array<double, 3> cross{
          a[first][1] * a[second][2] - a[first][2] * a[second][1],
          a[first][2] * a[second][0] - a[first][0] * a[second][2],
          a[first][0] * a[second][1] - a[first][1] * a[second][0]};
      const double cross_squared = cross[0] * cross[0] + cross[1] * cross[1] +
                                   cross[2] * cross[2];
      if (cross_squared > best_cross_squared) {
        best_cross_squared = cross_squared;
        axis = cross;
      }
    }
  }
  const double axis_norm = std::sqrt(best_cross_squared);
  if (!(axis_norm > 1e-10 * matrix_scale * matrix_scale))
    fail("quadric is not a rank-two cylinder quadratic");
  for (double& value : axis) value /= axis_norm;
  for (const double value : axis) {
    if (std::abs(value) > 1e-12) {
      if (value < 0.0) for (double& component : axis) component = -component;
      break;
    }
  }

  std::array<double, 3> radial_x{};
  int smallest_axis_component = 0;
  for (int i = 1; i < 3; ++i)
    if (std::abs(axis[i]) < std::abs(axis[smallest_axis_component]))
      smallest_axis_component = i;
  radial_x[smallest_axis_component] = 1.0;
  const double projection = radial_x[0] * axis[0] + radial_x[1] * axis[1] +
                            radial_x[2] * axis[2];
  for (int i = 0; i < 3; ++i) radial_x[i] -= projection * axis[i];
  const double radial_norm = std::sqrt(radial_x[0] * radial_x[0] +
                                       radial_x[1] * radial_x[1] +
                                       radial_x[2] * radial_x[2]);
  for (double& value : radial_x) value /= radial_norm;
  const std::array<double, 3> radial_y{
      axis[1] * radial_x[2] - axis[2] * radial_x[1],
      axis[2] * radial_x[0] - axis[0] * radial_x[2],
      axis[0] * radial_x[1] - axis[1] * radial_x[0]};

  const auto bilinear = [&](const std::array<double, 3>& left,
                            const std::array<double, 3>& right) {
    double value = 0.0;
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column)
        value += left[row] * a[row][column] * right[column];
    return value;
  };
  const double b00 = bilinear(radial_x, radial_x);
  const double b01 = bilinear(radial_x, radial_y);
  const double b11 = bilinear(radial_y, radial_y);
  const double half_difference = 0.5 * (b00 - b11);
  const double off_diagonal = std::hypot(half_difference, b01);
  const double lambda_x = 0.5 * (b00 + b11) + off_diagonal;
  const double lambda_y = 0.5 * (b00 + b11) - off_diagonal;
  if (!std::isfinite(lambda_x) || !std::isfinite(lambda_y))
    fail("rank-two quadric has non-finite transverse eigenvalues");
  if (requireSaddle) {
    if (lambda_x <= 1e-12 * matrix_scale ||
        lambda_y >= -1e-12 * matrix_scale)
      fail("hyperbolic-paraboloid quadratic block must have one positive and one negative eigenvalue");
  } else if (lambda_y <= 1e-12 * matrix_scale) {
    fail("elliptic-cylinder quadratic block must have rank two and positive transverse eigenvalues");
  }
  const double angle = 0.5 * std::atan2(2.0 * b01, b00 - b11);
  const std::array<double, 3> eigen_x{
      std::cos(angle) * radial_x[0] + std::sin(angle) * radial_y[0],
      std::cos(angle) * radial_x[1] + std::sin(angle) * radial_y[1],
      std::cos(angle) * radial_x[2] + std::sin(angle) * radial_y[2]};
  const std::array<double, 3> eigen_y{
      -std::sin(angle) * radial_x[0] + std::cos(angle) * radial_y[0],
      -std::sin(angle) * radial_x[1] + std::cos(angle) * radial_y[1],
      -std::sin(angle) * radial_x[2] + std::cos(angle) * radial_y[2]};
  const double matrix_tolerance = 1e-10 * matrix_scale;
  for (int row = 0; row < 3; ++row) {
    for (int column = 0; column < 3; ++column) {
      const double expected = lambda_x * eigen_x[row] * eigen_x[column] +
                              lambda_y * eigen_y[row] * eigen_y[column];
      if (std::abs(a[row][column] - expected) > matrix_tolerance)
        fail("quadric is not a rank-two elliptic-cylinder quadratic");
    }
  }
  const double axial_linear = b[0] * axis[0] + b[1] * axis[1] + b[2] * axis[2];
  const double target[3] = {-0.5 * b[0], -0.5 * b[1], -0.5 * b[2]};
  const double coordinate_x = (target[0] * eigen_x[0] + target[1] * eigen_x[1] +
                               target[2] * eigen_x[2]) / lambda_x;
  const double coordinate_y = (target[0] * eigen_y[0] + target[1] * eigen_y[1] +
                               target[2] * eigen_y[2]) / lambda_y;
  const double center[3] = {coordinate_x * eigen_x[0] + coordinate_y * eigen_y[0],
                            coordinate_x * eigen_x[1] + coordinate_y * eigen_y[1],
                            coordinate_x * eigen_x[2] + coordinate_y * eigen_y[2]};
  double center_quadratic = 0.0;
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column)
      center_quadratic += center[row] * a[row][column] * center[column];
  const double center_level = center_quadratic - q[0];
  if (!std::isfinite(center_level) || !std::isfinite(axial_linear))
    fail("rank-two quadric has non-finite axial level data");
  return {gp_Pnt(center[0], center[1], center[2]),
          gp_Dir(axis[0], axis[1], axis[2]),
          gp_Dir(eigen_x[0], eigen_x[1], eigen_x[2]),
          gp_Dir(eigen_y[0], eigen_y[1], eigen_y[2]),
          lambda_x, lambda_y, center_level, axial_linear};
}

TopoDS_Shape finite_elliptic_cylinder_carrier(
    const std::array<double, 10>& coefficients, const TopoDS_Shape& design) {
  const RankTwoQuadricFrame cylinder =
      rank_two_quadric_from_solver_coefficients(coefficients);
  double linear_scale = 0.0;
  for (int index = 1; index <= 6; ++index)
    linear_scale = std::max(linear_scale, std::abs(coefficients[index]));
  if (std::abs(cylinder.axial_linear) >
      1e-10 * std::max(linear_scale, 1e-300))
    fail("elliptic-cylinder quadric has a linear term along its axis");
  if (!(cylinder.center_level > 0.0))
    fail("elliptic-cylinder quadric has no positive finite cross-section");
  const double radius_x = std::sqrt(cylinder.center_level / cylinder.lambda_x);
  const double radius_y = std::sqrt(cylinder.center_level / cylinder.lambda_y);
  if (!std::isfinite(radius_x) || !std::isfinite(radius_y) ||
      !(radius_x > 0.0) || !(radius_y > 0.0))
    fail("elliptic-cylinder quadric has invalid radii");

  Bnd_Box bounds;
  BRepBndLib::AddOptimal(design, bounds, Standard_False, Standard_True);
  if (bounds.IsVoid() || bounds.IsOpenXmin() || bounds.IsOpenXmax() ||
      bounds.IsOpenYmin() || bounds.IsOpenYmax() ||
      bounds.IsOpenZmin() || bounds.IsOpenZmax())
    fail("elliptic-cylinder Boolean requires a finite design-domain box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const gp_Vec axis(cylinder.axis);
  const gp_Pnt center = cylinder.axis_point;
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const double projection = gp_Vec(center, corner).Dot(axis);
        minimum = std::min(minimum, projection);
        maximum = std::max(maximum, projection);
      }
  const double projected_span = maximum - minimum;
  if (!std::isfinite(projected_span) || !(projected_span > 0.0))
    fail("elliptic-cylinder axis has no finite design-domain span");
  // Both planar caps lie outside the complete design AABB projection; they
  // therefore cannot affect the represented unbounded-cylinder intersection.
  const double cap_margin = std::max(1e-7, 0.1 * projected_span);
  const double carrier_height = projected_span + 2.0 * cap_margin;
  const double map[3][3] = {
      {radius_x * cylinder.radial_x.X(), radius_y * cylinder.radial_y.X(),
       carrier_height * cylinder.axis.X()},
      {radius_x * cylinder.radial_x.Y(), radius_y * cylinder.radial_y.Y(),
       carrier_height * cylinder.axis.Y()},
      {radius_x * cylinder.radial_x.Z(), radius_y * cylinder.radial_y.Z(),
       carrier_height * cylinder.axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_Pnt start = center.Translated(axis * (minimum - cap_margin));
  const gp_GTrsf transform(vectorial,
      gp_XYZ(start.X(), start.Y(), start.Z()));
  const TopoDS_Shape unit = BRepPrimAPI_MakeCylinder(1.0, 1.0).Shape();
  BRepBuilderAPI_GTransform affine(unit, transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not construct the finite elliptic-cylinder carrier");
  return affine.Shape();
}

struct HyperboloidQuadricFrame {
  gp_Pnt center;
  gp_Dir radial_x;
  gp_Dir radial_y;
  gp_Dir axis;
  double lambda_x;
  double lambda_y;
  double lambda_axis;
  double center_level;
};

HyperboloidQuadricFrame hyperboloid_from_solver_coefficients(
    const std::array<double, 10>& q, bool expectTwoSheet = false) {
  double coefficient_norm_squared = 0.0;
  for (const double value : q) coefficient_norm_squared += value * value;
  const double coefficient_norm = std::sqrt(coefficient_norm_squared);
  if (!std::isfinite(coefficient_norm) ||
      std::abs(coefficient_norm - 1.0) > 1e-10)
    fail("solver quadric coefficients must have unit Euclidean norm");

  // Symmetric Jacobi diagonalization for the solver's quadratic block.
  // Keeping this small implementation local avoids pulling a second numeric
  // ABI into the pinned OCCT-only probe.
  double a[3][3] = {{q[4], q[9], q[8]},
                    {q[9], q[5], q[7]},
                    {q[8], q[7], q[6]}};
  double eigenvectors[3][3] = {{1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}};
  double matrix_scale = 0.0;
  for (const auto& row : a)
    for (const double value : row)
      matrix_scale = std::max(matrix_scale, std::abs(value));
  if (!(matrix_scale > 0.0)) fail("hyperboloid quadratic block is zero");
  for (int iteration = 0; iteration < 64; ++iteration) {
    int p = 0, r = 1;
    double largest = std::abs(a[p][r]);
    for (const auto pair : {std::array<int, 2>{0, 2},
                            std::array<int, 2>{1, 2}}) {
      if (std::abs(a[pair[0]][pair[1]]) > largest) {
        p = pair[0];
        r = pair[1];
        largest = std::abs(a[p][r]);
      }
    }
    if (largest <= 2e-16 * matrix_scale) break;
    const double tau = (a[r][r] - a[p][p]) / (2.0 * a[p][r]);
    const double t = std::copysign(1.0, tau) /
                     (std::abs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;
    const double app = a[p][p], arr = a[r][r], apr = a[p][r];
    a[p][p] = app - t * apr;
    a[r][r] = arr + t * apr;
    a[p][r] = a[r][p] = 0.0;
    for (int k = 0; k < 3; ++k) {
      if (k == p || k == r) continue;
      const double akp = a[k][p], akr = a[k][r];
      a[k][p] = a[p][k] = c * akp - s * akr;
      a[k][r] = a[r][k] = s * akp + c * akr;
    }
    for (int k = 0; k < 3; ++k) {
      const double vkp = eigenvectors[k][p], vkr = eigenvectors[k][r];
      eigenvectors[k][p] = c * vkp - s * vkr;
      eigenvectors[k][r] = s * vkp + c * vkr;
    }
  }

  std::array<int, 2> positive_indices{};
  int positive_count = 0, negative_index = -1;
  const double eigenvalue_tolerance = 1e-10 * matrix_scale;
  for (int column = 0; column < 3; ++column) {
    if (a[column][column] > eigenvalue_tolerance) {
      if (positive_count >= 2)
        fail("quadric does not have the one-sheet hyperboloid signature (+,+,-)");
      positive_indices[positive_count++] = column;
    } else if (a[column][column] < -eigenvalue_tolerance) {
      if (negative_index >= 0)
        fail("quadric does not have the one-sheet hyperboloid signature (+,+,-)");
      negative_index = column;
    } else {
      fail("one-sheet hyperboloid quadratic block must have full rank");
    }
  }
  if (positive_count != 2 || negative_index < 0)
    fail("quadric does not have the one-sheet hyperboloid signature (+,+,-)");
  if (a[positive_indices[1]][positive_indices[1]] >
      a[positive_indices[0]][positive_indices[0]])
    std::swap(positive_indices[0], positive_indices[1]);

  auto vector_for = [&](int index) {
    std::array<double, 3> value{eigenvectors[0][index],
                                eigenvectors[1][index],
                                eigenvectors[2][index]};
    int largest_component = 0;
    for (int i = 1; i < 3; ++i)
      if (std::abs(value[i]) > std::abs(value[largest_component]))
        largest_component = i;
    if (value[largest_component] < 0.0)
      for (double& component : value) component = -component;
    return value;
  };
  auto ex = vector_for(positive_indices[0]);
  auto ey = vector_for(positive_indices[1]);
  auto ez = vector_for(negative_index);
  const std::array<double, 3> ex_cross_ey{
      ex[1] * ey[2] - ex[2] * ey[1],
      ex[2] * ey[0] - ex[0] * ey[2],
      ex[0] * ey[1] - ex[1] * ey[0]};
  if (ex_cross_ey[0] * ez[0] + ex_cross_ey[1] * ez[1] +
          ex_cross_ey[2] * ez[2] < 0.0)
    for (double& component : ez) component = -component;

  // Reject an incomplete/unstable diagonalization rather than build a carrier
  // from eigenvectors that do not reconstruct the requested solver quadric.
  const double eigenvalues[3] = {a[positive_indices[0]][positive_indices[0]],
                                 a[positive_indices[1]][positive_indices[1]],
                                 a[negative_index][negative_index]};
  const std::array<std::array<double, 3>, 3> basis{ex, ey, ez};
  const double original_a[3][3] = {{q[4], q[9], q[8]},
                                    {q[9], q[5], q[7]},
                                    {q[8], q[7], q[6]}};
  for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column) {
      double reconstructed = 0.0;
      for (int k = 0; k < 3; ++k)
        reconstructed += eigenvalues[k] * basis[k][row] * basis[k][column];
      if (std::abs(reconstructed - original_a[row][column]) >
          2e-10 * matrix_scale)
        fail("hyperboloid eigensystem failed quadratic reconstruction check");
    }

  // Classify the actual represented polynomial, not a rounded center or an
  // absolute "near zero" cutoff. A positive 1e-36 level is not a cone.
  // Rational elimination preserves cancellation in q(c); the CAD geometry
  // and integration acceptance tolerances are not changed by this predicate.
  mpq_class augmented[3][4];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) augmented[i][j] = original_a[i][j];
    augmented[i][3] = -mpq_class(q[i+1]) / 2;
  }
  for (int i = 0; i < 3; ++i) {
    int pivot = i;
    while (pivot < 3 && augmented[pivot][i] == 0) ++pivot;
    if (pivot == 3) fail("hyperboloid quadratic block is exactly singular");
    for (int j = 0; j < 4; ++j) std::swap(augmented[i][j], augmented[pivot][j]);
    const mpq_class divisor = augmented[i][i];
    for (int j = 0; j < 4; ++j) augmented[i][j] /= divisor;
    for (int row = 0; row < 3; ++row) if (row != i) {
      const mpq_class factor = augmented[row][i];
      for (int j = 0; j < 4; ++j) augmented[row][j] -= factor * augmented[i][j];
    }
  }
  // A*c=-b/2 gives -q(c)=-q0-b.c/2 exactly.
  mpq_class exact_level = -mpq_class(q[0]);
  std::array<double, 3> center{};
  for (int i = 0; i < 3; ++i) {
    exact_level -= mpq_class(q[i+1]) * augmented[i][3] / 2;
    center[i] = augmented[i][3].get_d();
    if (!std::isfinite(center[i])) fail("hyperboloid center is not representable");
  }
  if (expectTwoSheet ? exact_level >= 0 : exact_level <= 0)
    fail(expectTwoSheet
             ? "two-sheet hyperboloid sublevel must exclude its center with negative level"
             : "one-sheet hyperboloid sublevel must contain its center with positive level");
  const double center_level = exact_level.get_d();
  if (!std::isfinite(center_level) || center_level == 0.0)
    fail("nonzero hyperboloid center level is not representable by OCCT");
  return {gp_Pnt(center[0], center[1], center[2]),
          gp_Dir(ex[0], ex[1], ex[2]), gp_Dir(ey[0], ey[1], ey[2]),
          gp_Dir(ez[0], ez[1], ez[2]), eigenvalues[0], eigenvalues[1],
          eigenvalues[2], center_level};
}

TopoDS_Shape mesh_to_solid(const std::string& path, double coordinate_scale) {
  std::ifstream input(path);
  if (!input) fail("cannot open input OBJ mesh");
  std::vector<gp_Pnt> vertices;
  std::vector<std::array<int, 3>> triangles;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    std::istringstream fields(line);
    std::string tag;
    fields >> tag;
    if (tag.empty() || tag[0] == '#') continue;
    if (tag == "v") {
      double x, y, z;
      if (!(fields >> x >> y >> z) || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        fail("malformed or non-finite OBJ vertex at line " + std::to_string(line_number));
      x *= coordinate_scale;
      y *= coordinate_scale;
      z *= coordinate_scale;
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        fail("coordinate scale produced a non-finite OBJ vertex at line " + std::to_string(line_number));
      vertices.emplace_back(x, y, z);
    } else if (tag == "f") {
      std::array<int, 3> face{};
      std::string token;
      for (int corner = 0; corner < 3; ++corner) {
        if (!(fields >> token)) fail("non-triangular OBJ face at line " + std::to_string(line_number));
        const std::size_t slash = token.find('/');
        const std::string index_text = token.substr(0, slash);
        std::size_t consumed = 0;
        const long parsed = std::stol(index_text, &consumed);
        if (consumed != index_text.size() || parsed <= 0 || parsed > std::numeric_limits<int>::max())
          fail("unsupported OBJ face index at line " + std::to_string(line_number));
        face[corner] = static_cast<int>(parsed - 1);
      }
      std::string extra;
      if (fields >> extra) fail("non-triangular OBJ face at line " + std::to_string(line_number));
      triangles.push_back(face);
    } else {
      fail("unsupported OBJ record '" + tag + "' at line " + std::to_string(line_number));
    }
  }
  if (vertices.empty() || triangles.empty()) fail("OBJ mesh is empty");

  BRepBuilderAPI_Sewing sewing(0.0, Standard_True, Standard_True, Standard_True, Standard_False);
  for (const auto& triangle : triangles) {
    for (const int index : triangle) {
      if (index < 0 || static_cast<std::size_t>(index) >= vertices.size())
        fail("OBJ face index is outside the vertex array");
    }
    const gp_Pnt& a = vertices[triangle[0]];
    const gp_Pnt& b = vertices[triangle[1]];
    const gp_Pnt& c = vertices[triangle[2]];
    if (a.Distance(b) == 0.0 || b.Distance(c) == 0.0 || c.Distance(a) == 0.0 ||
        (b.XYZ() - a.XYZ()).CrossSquareMagnitude(c.XYZ() - a.XYZ()) == 0.0)
      fail("OBJ contains a degenerate triangle");
    BRepBuilderAPI_MakePolygon polygon;
    polygon.Add(a);
    polygon.Add(b);
    polygon.Add(c);
    polygon.Close();
    if (!polygon.IsDone()) fail("OCCT could not build an OBJ triangle wire");
    BRepBuilderAPI_MakeFace face(polygon.Wire());
    if (!face.IsDone()) fail("OCCT could not build a planar OBJ triangle face");
    sewing.Add(face.Face());
  }
  sewing.Perform();
  if (sewing.NbFreeEdges() != 0 || sewing.NbMultipleEdges() != 0 ||
      sewing.NbDegeneratedShapes() != 0 || sewing.NbDeletedFaces() != 0)
    fail("OBJ must preserve every closed shell: sewing found open/nonmanifold edges or removed degenerate faces");
  const TopoDS_Shape sewn = sewing.SewedShape();
  std::vector<TopoDS_Shell> shells;
  for (TopExp_Explorer it(sewn, TopAbs_SHELL); it.More(); it.Next()) {
    shells.push_back(TopoDS::Shell(it.Current()));
  }
  if (shells.empty()) fail("OBJ did not sew to any closed shell");
  std::vector<TopoDS_Solid> bounded;
  std::vector<double> signed_volumes;
  std::vector<gp_Pnt> witnesses;
  for (const auto& shell : shells) {
    if (!BRep_Tool::IsClosed(shell)) fail("sewn OBJ shell is not closed");
    BRepBuilderAPI_MakeSolid make_solid(shell);
    if (!make_solid.IsDone()) fail("OCCT could not convert sewn OBJ shell to a solid");
    TopoDS_Solid solid = make_solid.Solid();
    signed_volumes.push_back(volume_moments(solid).volume);
    if (!std::isfinite(signed_volumes.back()) || signed_volumes.back() == 0.)
      fail("OBJ shell has zero or nonfinite signed volume");
    if (!BRepLib::OrientClosedSolid(solid)) fail("OCCT could not orient the closed OBJ solid");
    if (!BRepCheck_Analyzer(solid, Standard_True).IsValid())
      fail("sewn OBJ solid failed OCCT topology validation");
    if (!(volume_moments(solid).volume > 0.0)) fail("sewn OBJ solid has non-positive volume");
    bounded.push_back(solid);
    TopExp_Explorer vertex(shell, TopAbs_VERTEX);
    if (!vertex.More()) fail("OBJ shell has no point for containment classification");
    witnesses.push_back(BRep_Tool::Pnt(TopoDS::Vertex(vertex.Current())));
  }
  // Preserve the historical single-shell orientation normalization. For a
  // multi-shell material, orientation carries cavity semantics and must agree
  // with the containment forest; never fill a hole by orienting every shell
  // independently and returning their union.
  if (shells.size() == 1) return bounded.front();
  const std::size_t count = shells.size();
  std::vector<std::vector<bool>> inside(count, std::vector<bool>(count, false));
  std::vector<int> depth(count, 0), parent(count, -1);
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = i+1; j < count; ++j) {
      BRepExtrema_DistShapeShape distance(shells[i], shells[j]);
      if (!distance.IsDone() || !std::isfinite(distance.Value()) ||
          distance.Value() <= Precision::Confusion())
        fail("OBJ shell boundaries intersect, touch or are unresolved at OCCT precision");
      // Once disjoint connected boundaries are established, a boundary point
      // determines their containment. ON/UNKNOWN is not silently classified.
      BRepClass3d_SolidClassifier i_in_j(bounded[j], witnesses[i], Precision::Confusion());
      BRepClass3d_SolidClassifier j_in_i(bounded[i], witnesses[j], Precision::Confusion());
      if ((i_in_j.State() != TopAbs_IN && i_in_j.State() != TopAbs_OUT) ||
          (j_in_i.State() != TopAbs_IN && j_in_i.State() != TopAbs_OUT))
        fail("OBJ shell containment is unresolved at OCCT precision");
      inside[i][j] = i_in_j.State() == TopAbs_IN;
      inside[j][i] = j_in_i.State() == TopAbs_IN;
      if (inside[i][j] && inside[j][i]) fail("OBJ shells have inconsistent mutual containment");
    }
  }
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = 0; j < count; ++j) if (inside[i][j]) ++depth[i];
    if ((signed_volumes[i] > 0.) != (depth[i] % 2 == 0))
      fail("OBJ shell winding disagrees with nested material/cavity parity");
  }
  for (std::size_t i = 0; i < count; ++i) {
    for (std::size_t j = 0; j < count; ++j) {
      if (inside[i][j] && depth[j] == depth[i]-1) {
        if (parent[i] != -1) fail("OBJ containment does not form a shell forest");
        parent[i] = static_cast<int>(j);
      }
    }
    if (depth[i] > 0 && parent[i] < 0) fail("OBJ shell has no immediate container");
  }
  BRep_Builder builder;
  TopoDS_Compound material;
  builder.MakeCompound(material);
  for (std::size_t i = 0; i < count; ++i) if (depth[i] % 2 == 0) {
    BRepBuilderAPI_MakeSolid part(shells[i]);
    for (std::size_t j = 0; j < count; ++j)
      if (parent[j] == static_cast<int>(i)) part.Add(shells[j]);
    TopoDS_Solid solid = part.Solid();
    if (!BRepCheck_Analyzer(solid, Standard_True).IsValid() ||
        !(volume_moments(solid).volume > 0.))
      fail("nested OBJ material failed OCCT solid validation");
    builder.Add(material, solid);
  }
  double expected_volume = 0.;
  for (double volume : signed_volumes) expected_volume += volume;
  const double actual_volume = volume_moments(material).volume;
  if (!(expected_volume > 0.) || !std::isfinite(actual_volume) ||
      std::abs(actual_volume-expected_volume) > 1e-10*expected_volume)
    fail("nested OBJ material volume differs from signed shell contribution");
  return material;
}

void write_step(const TopoDS_Shape& shape, const std::string& step_path) {
  STEPControl_Writer writer;
  if (writer.Transfer(shape, STEPControl_AsIs) != IFSelect_RetDone)
    fail("failed to transfer material artifact to STEP writer");
  std::ofstream step_stream(step_path, std::ios::out | std::ios::binary | std::ios::trunc);
  if (!step_stream) fail("failed to open STEP material artifact for writing");
  if (writer.WriteStream(step_stream) != IFSelect_RetDone || !step_stream)
    fail("failed to write STEP material artifact");
}

TopoDS_Shape read_brep(const std::string& path) {
  BRep_Builder builder;
  TopoDS_Shape shape;
  if (!BRepTools::Read(shape, path.c_str(), builder)) fail("failed to read B-rep material artifact");
  return shape;
}

TopoDS_Shape read_step(const std::string& path) {
  STEPControl_Reader reader;
  if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) fail("failed to read STEP material artifact");
  if (reader.TransferRoots() <= 0) fail("STEP reader did not transfer any roots");
  return reader.OneShape();
}

void validate_roundtrip_solid(const TopoDS_Shape& shape, const std::string& stage) {
  require_valid_topology(shape,stage,stage + " CAD shape failed OCCT topology validation");
  const Moments moments = volume_moments(shape);
  if (!(moments.volume > 0.0) || !std::isfinite(moments.volume))
    fail(stage + " CAD shape has non-positive or non-finite volume");
  for (double value : moments.raw) {
    if (!std::isfinite(value)) fail(stage + " CAD shape has non-finite raw moments");
  }
}

int run_mesh_roundtrip(int argc, char** argv) {
  if (argc != 6) fail("usage: occt_step_roundtrip --mesh INPUT.obj COORDINATE_SCALE MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  const std::string brep_path(argv[4]);
  const std::string step_path(argv[5]);
  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape original = mesh_to_solid(obj_path, coordinate_scale);
  if (!BRepTools::Write(original, brep_path.c_str())) fail("failed to write B-rep material artifact");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(original, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(original, "source mesh");
  validate_roundtrip_solid(brep_shape, "B-rep reload");
  validate_roundtrip_solid(step_shape, "STEP reload");
  const Moments original_moments = volume_moments(original);
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);
  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_METHOD_VALIDATION_NOT_AUTHOR_CAD\","
               "\"status\":\"SOLVER_MESH_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED\","
               "\"coordinate_scale\":"
            << coordinate_scale << ",\"raw_moments\":{";
  const char* names[] = {"original", "brep", "step"};
  const Moments* moments[] = {&original_moments, &brep_moments, &step_moments};
  for (int stage = 0; stage < 3; ++stage) {
    if (stage != 0) std::cout << ',';
    std::cout << '\"' << names[stage] << "\":";
    emit_raw_moments(*moments[stage]);
  }
  std::cout << "}}\n";
  return 0;
}

int run_solver_hyperbolic_paraboloid_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-hyperbolic-paraboloid INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const RankTwoQuadricFrame quadric =
      rank_two_quadric_from_solver_coefficients(coefficients, true);
  const double matrix_scale = std::max(std::abs(quadric.lambda_x),
                                       std::abs(quadric.lambda_y));
  if (std::abs(quadric.axial_linear) <= 1e-10 * matrix_scale)
    fail("hyperbolic paraboloid requires a nonzero linear component along its null axis");
  const std::filesystem::path brep_path =
      std::filesystem::absolute(argv[14]).lexically_normal();
  const std::filesystem::path step_path =
      std::filesystem::absolute(argv[15]).lexically_normal();
  if (brep_path == step_path)
    fail("B-rep and STEP output paths must differ");
  std::error_code filesystem_error;
  if (std::filesystem::exists(brep_path, filesystem_error) || filesystem_error ||
      std::filesystem::exists(step_path, filesystem_error) || filesystem_error)
    fail("CAD output path already exists; refusing to overwrite existing evidence");
  if (!std::filesystem::is_directory(brep_path.parent_path(), filesystem_error) ||
      filesystem_error ||
      !std::filesystem::is_directory(step_path.parent_path(), filesystem_error) ||
      filesystem_error)
    fail("CAD output parent directories must already exist");

  const double scale_x = std::sqrt(std::abs(quadric.axial_linear) /
                                   std::abs(quadric.lambda_x));
  const double scale_y = std::sqrt(std::abs(quadric.axial_linear) /
                                   std::abs(quadric.lambda_y));
  if (!std::isfinite(scale_x) || !std::isfinite(scale_y) ||
      !(scale_x > 0.0) || !(scale_y > 0.0))
    fail("hyperbolic-paraboloid transverse scales are invalid");

  // After this sign choice the implicit equation is u^2-v^2+w=0 in the
  // canonical frame. Its q<=0 component is the region below w=-u^2+v^2.
  const gp_Dir canonical_axis = quadric.axial_linear < 0.0
      ? gp_Dir(-quadric.axis.X(), -quadric.axis.Y(), -quadric.axis.Z())
      : quadric.axis;
  const double vertex_offset = quadric.center_level / quadric.axial_linear;
  const gp_Pnt vertex = quadric.axis_point.Translated(
      gp_Vec(quadric.axis) * vertex_offset);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  Bnd_Box outer_bounds;
  BRepBndLib::AddOptimal(outer, outer_bounds, Standard_False, Standard_True);
  if (outer_bounds.IsVoid() || outer_bounds.IsOpenXmin() || outer_bounds.IsOpenXmax() ||
      outer_bounds.IsOpenYmin() || outer_bounds.IsOpenYmax() ||
      outer_bounds.IsOpenZmin() || outer_bounds.IsOpenZmax())
    fail("hyperbolic-paraboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  outer_bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);

  double u_extent = 0.0, v_extent = 0.0;
  double minimum_w = std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const gp_Vec relative(vertex, corner);
        const double u = relative.Dot(gp_Vec(quadric.radial_x)) / scale_x;
        const double v = relative.Dot(gp_Vec(quadric.radial_y)) / scale_y;
        const double w = relative.Dot(gp_Vec(canonical_axis));
        u_extent = std::max(u_extent, std::abs(u));
        v_extent = std::max(v_extent, std::abs(v));
        minimum_w = std::min(minimum_w, w);
      }
    }
  }
  const double domain_diagonal = std::sqrt((xmax-xmin)*(xmax-xmin) +
                                           (ymax-ymin)*(ymax-ymin) +
                                           (zmax-zmin)*(zmax-zmin));
  const double margin = std::max(1e-7, 0.1 * domain_diagonal);
  u_extent += std::max(1e-7, 0.1 * u_extent);
  v_extent += std::max(1e-7, 0.1 * v_extent);
  const double carrier_depth = v_extent*v_extent - minimum_w + margin;
  if (!std::isfinite(u_extent) || !std::isfinite(v_extent) ||
      !std::isfinite(carrier_depth) || !(u_extent > 0.0) ||
      !(v_extent > 0.0) || !(carrier_depth > 0.0))
    fail("could not construct a finite hyperbolic-paraboloid carrier outside the domain");

  // A tensor-product quadratic Bezier patch represents the saddle exactly:
  // u and v are linear in their parameters, and its height is -u^2+v^2.
  TColgp_Array2OfPnt poles(1, 3, 1, 3);
  const double u_poles[3] = {-u_extent, 0.0, u_extent};
  const double u_square_poles[3] = {-u_extent*u_extent,
                                     u_extent*u_extent,
                                    -u_extent*u_extent};
  const double v_poles[3] = {-v_extent, 0.0, v_extent};
  const double v_square_poles[3] = {v_extent*v_extent,
                                    -v_extent*v_extent,
                                     v_extent*v_extent};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      poles.SetValue(i + 1, j + 1,
                     gp_Pnt(u_poles[i], v_poles[j],
                            u_square_poles[i] + v_square_poles[j]));
  Handle(Geom_BezierSurface) saddle_surface = new Geom_BezierSurface(poles);
  BRepBuilderAPI_MakeFace saddle_face(saddle_surface, 1e-9);
  if (!saddle_face.IsDone())
    fail("failed to construct exact hyperbolic-paraboloid surface face");
  BRepPrimAPI_MakePrism prism(saddle_face.Face(), gp_Vec(0.0, 0.0, -carrier_depth),
                              Standard_True, Standard_True);
  if (!prism.IsDone() || prism.Shape().IsNull())
    fail("failed to construct finite hyperbolic-paraboloid carrier");

  const double map[3][3] = {
      {scale_x * quadric.radial_x.X(), scale_y * quadric.radial_y.X(), canonical_axis.X()},
      {scale_x * quadric.radial_x.Y(), scale_y * quadric.radial_y.Y(), canonical_axis.Y()},
      {scale_x * quadric.radial_x.Z(), scale_y * quadric.radial_y.Z(), canonical_axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_GTrsf transform(vectorial, gp_XYZ(vertex.X(), vertex.Y(), vertex.Z()));
  BRepBuilderAPI_GTransform affine(prism.Shape(), transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not place the saddle carrier in solver coordinates");
  const TopoDS_Shape cavity = affine.Shape();
  if (!BRepCheck_Analyzer(cavity, Standard_True).IsValid())
    fail("hyperbolic-paraboloid carrier failed OCCT topology validation");

  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("solver hyperbolic-paraboloid Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "hyperbolic-paraboloid Boolean result");
  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write hyperbolic-paraboloid B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "hyperbolic-paraboloid B-rep reload");
  validate_roundtrip_solid(step_shape, "hyperbolic-paraboloid STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);
  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_HYPERBOLIC_PARABOLOID_QUADRIC_STEP_BREP_VALIDATED\","
               "\"surface_representation\":\"exact_tensor_product_quadratic_bezier\","
               "\"vertex\":[" << vertex.X() << ',' << vertex.Y() << ',' << vertex.Z()
            << "],\"axis\":[" << canonical_axis.X() << ',' << canonical_axis.Y()
            << ',' << canonical_axis.Z()
            << "],\"carrier_depth\":" << carrier_depth
            << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

int run_solver_sphere_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-sphere INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const SphereQuadric sphere = sphere_from_solver_coefficients(coefficients);
  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  const TopoDS_Shape cavity = BRepPrimAPI_MakeSphere(sphere.center, sphere.radius).Shape();
  if (!BRepCheck_Analyzer(cavity, Standard_True).IsValid())
    fail("analytic OCCT sphere solid failed topology validation");
  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("analytic solver-sphere Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "analytic sphere Boolean result");

  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write analytic sphere B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "analytic sphere B-rep reload");
  validate_roundtrip_solid(step_shape, "analytic sphere STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);

  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_SPHERE_QUADRIC_STEP_BREP_VALIDATED\","
               "\"coordinate_scale\":" << coordinate_scale
            << ",\"sphere_center\":[" << sphere.center.X() << ','
            << sphere.center.Y() << ',' << sphere.center.Z()
            << "],\"sphere_radius\":" << sphere.radius << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

int run_solver_ellipsoid_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-ellipsoid INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const EllipsoidQuadric ellipsoid = ellipsoid_from_solver_coefficients(coefficients);
  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  const TopoDS_Shape unit_sphere = BRepPrimAPI_MakeSphere(gp_Pnt(0., 0., 0.), 1.).Shape();
  BRepBuilderAPI_GTransform affine(unit_sphere, ellipsoid.unit_sphere_map, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not construct the exact affine ellipsoid from the solver quadric");
  const TopoDS_Shape cavity = affine.Shape();
  if (!BRepCheck_Analyzer(cavity, Standard_True).IsValid())
    fail("affine ellipsoid B-rep failed OCCT topology validation");
  const double emitted_ellipsoid_volume = volume_moments(cavity).volume;
  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("solver ellipsoid Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "ellipsoid Boolean result");

  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write ellipsoid B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "ellipsoid B-rep reload");
  validate_roundtrip_solid(step_shape, "ellipsoid STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);

  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_ELLIPSOID_QUADRIC_STEP_BREP_VALIDATED\","
               "\"coordinate_scale\":" << coordinate_scale
            << ",\"ellipsoid_center\":[" << ellipsoid.center.X() << ','
            << ellipsoid.center.Y() << ',' << ellipsoid.center.Z()
            << "],\"ellipsoid_volume\":" << ellipsoid.volume
            << ",\"emitted_ellipsoid_volume\":" << emitted_ellipsoid_volume
            << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

int run_solver_cylinder_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-cylinder INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const RankTwoQuadricFrame cylinder = rank_two_quadric_from_solver_coefficients(coefficients);
  const double linear_scale = std::max({std::abs(coefficients[1]),
                                        std::abs(coefficients[2]),
                                        std::abs(coefficients[3]),
                                        std::abs(coefficients[4]),
                                        std::abs(coefficients[5]),
                                        std::abs(coefficients[6])});
  if (std::abs(cylinder.axial_linear) > 1e-10 * linear_scale)
    fail("quadric has a linear component along its cylinder axis");
  if (!(cylinder.center_level > 0.0))
    fail("solver elliptic-cylinder quadric must enclose a positive finite cross-section");
  const double radius_x = std::sqrt(cylinder.center_level / cylinder.lambda_x);
  const double radius_y = std::sqrt(cylinder.center_level / cylinder.lambda_y);
  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  Bnd_Box outer_bounds;
  BRepBndLib::AddOptimal(outer, outer_bounds, Standard_False, Standard_True);
  if (outer_bounds.IsVoid() || outer_bounds.IsOpenXmin() || outer_bounds.IsOpenXmax() ||
      outer_bounds.IsOpenYmin() || outer_bounds.IsOpenYmax() ||
      outer_bounds.IsOpenZmin() || outer_bounds.IsOpenZmax())
    fail("cylinder Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  outer_bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);

  const gp_Vec axis(cylinder.axis);
  const gp_Pnt center = cylinder.axis_point;
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const double projection = gp_Vec(center, corner).Dot(axis);
        minimum = std::min(minimum, projection);
        maximum = std::max(maximum, projection);
      }
    }
  }
  const double projected_span = maximum - minimum;
  if (!std::isfinite(projected_span) || !(projected_span > 0.0))
    fail("cylinder axis has no finite overlap interval with the design domain");
  // The finite carrier's planar caps are placed strictly outside the design
  // AABB projection, so the Boolean retains the unbounded cylinder's exact
  // intersection with every point in the domain.
  const double cap_margin = std::max(1e-7, 0.1 * projected_span);
  const double carrier_height = projected_span + 2.0 * cap_margin;
  const gp_Pnt carrier_start = center.Translated(axis * (minimum - cap_margin));
  const double map[3][3] = {
      {radius_x * cylinder.radial_x.X(),
       radius_y * cylinder.radial_y.X(), carrier_height * cylinder.axis.X()},
      {radius_x * cylinder.radial_x.Y(),
       radius_y * cylinder.radial_y.Y(), carrier_height * cylinder.axis.Y()},
      {radius_x * cylinder.radial_x.Z(),
       radius_y * cylinder.radial_y.Z(), carrier_height * cylinder.axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_GTrsf transform(vectorial, gp_XYZ(carrier_start.X(), carrier_start.Y(),
                                              carrier_start.Z()));
  const TopoDS_Shape unit_cylinder = BRepPrimAPI_MakeCylinder(1.0, 1.0).Shape();
  BRepBuilderAPI_GTransform affine(unit_cylinder, transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not construct the elliptic-cylinder carrier from the solver quadric");
  const TopoDS_Shape cavity = affine.Shape();
  if (!BRepCheck_Analyzer(cavity, Standard_True).IsValid())
    fail("analytic OCCT finite elliptic-cylinder carrier failed topology validation");
  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("solver cylinder Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "cylinder Boolean result");

  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write cylinder B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "cylinder B-rep reload");
  validate_roundtrip_solid(step_shape, "cylinder STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);

  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_ELLIPTIC_CYLINDER_QUADRIC_STEP_BREP_VALIDATED\","
               "\"coordinate_scale\":" << coordinate_scale
            << ",\"cylinder_center\":[" << cylinder.axis_point.X() << ','
            << cylinder.axis_point.Y() << ',' << cylinder.axis_point.Z()
            << "],\"cylinder_axis\":[" << cylinder.axis.X() << ','
            << cylinder.axis.Y() << ',' << cylinder.axis.Z()
            << "],\"cylinder_radii\":[" << radius_x << ','
            << radius_y << ']'
            << ",\"carrier_height\":" << carrier_height
            << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

int run_solver_paraboloid_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-paraboloid INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const RankTwoQuadricFrame quadric = rank_two_quadric_from_solver_coefficients(coefficients);
  if (std::abs(quadric.axial_linear) <= 1e-10 * std::max({
          std::abs(coefficients[1]), std::abs(coefficients[2]), std::abs(coefficients[3]),
          std::abs(coefficients[4]), std::abs(coefficients[5]), std::abs(coefficients[6])}))
    fail("paraboloid quadric requires a nonzero linear component along its axis");
  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);

  const double transverse_scale_x =
      std::sqrt(std::abs(quadric.axial_linear) / quadric.lambda_x);
  const double transverse_scale_y =
      std::sqrt(std::abs(quadric.axial_linear) / quadric.lambda_y);
  if (!std::isfinite(transverse_scale_x) || !std::isfinite(transverse_scale_y) ||
      !(transverse_scale_x > 0.0) || !(transverse_scale_y > 0.0))
    fail("paraboloid transverse scaling is non-finite or degenerate");

  const gp_Dir canonical_axis = quadric.axial_linear < 0.0
                                   ? quadric.axis
                                   : gp_Dir(-quadric.axis.X(), -quadric.axis.Y(),
                                            -quadric.axis.Z());
  const double vertex_offset = quadric.center_level / quadric.axial_linear;
  const gp_Pnt vertex = quadric.axis_point.Translated(gp_Vec(quadric.axis) * vertex_offset);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  Bnd_Box outer_bounds;
  BRepBndLib::AddOptimal(outer, outer_bounds, Standard_False, Standard_True);
  if (outer_bounds.IsVoid() || outer_bounds.IsOpenXmin() || outer_bounds.IsOpenXmax() ||
      outer_bounds.IsOpenYmin() || outer_bounds.IsOpenYmax() ||
      outer_bounds.IsOpenZmin() || outer_bounds.IsOpenZmax())
    fail("paraboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  outer_bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);

  double maximum_radius_squared = 0.0;
  double maximum_axial = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x) {
    for (int y = 0; y < 2; ++y) {
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const gp_Vec relative(vertex, corner);
        const double local_x = relative.Dot(gp_Vec(quadric.radial_x)) / transverse_scale_x;
        const double local_y = relative.Dot(gp_Vec(quadric.radial_y)) / transverse_scale_y;
        const double local_z = relative.Dot(gp_Vec(canonical_axis));
        maximum_radius_squared = std::max(maximum_radius_squared,
                                          local_x * local_x + local_y * local_y);
        maximum_axial = std::max(maximum_axial, local_z);
      }
    }
  }
  const double domain_diagonal = std::sqrt((xmax - xmin) * (xmax - xmin) +
                                           (ymax - ymin) * (ymax - ymin) +
                                           (zmax - zmin) * (zmax - zmin));
  const double geometric_margin = std::max(1e-7, 0.1 * domain_diagonal);
  const double radial_limit = std::sqrt(maximum_radius_squared) +
                              std::max(1e-7, 0.1 * std::sqrt(maximum_radius_squared));
  const double cap_height = std::max(maximum_axial, radial_limit * radial_limit) +
                            geometric_margin;
  if (!std::isfinite(radial_limit) || !std::isfinite(cap_height) ||
      !(radial_limit > 0.0) || !(cap_height > maximum_axial))
    fail("could not construct a finite paraboloid carrier outside the design domain");

  const gp_Pnt p0(0.0, 0.0, 0.0);
  const gp_Pnt p1(0.0, 0.0, cap_height);
  const gp_Pnt p2(radial_limit, 0.0, cap_height);
  const gp_Pnt p3(radial_limit, 0.0, radial_limit * radial_limit);
  BRepBuilderAPI_MakeWire profile_wire;
  profile_wire.Add(BRepBuilderAPI_MakeEdge(p0, p1).Edge());
  profile_wire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
  profile_wire.Add(BRepBuilderAPI_MakeEdge(p2, p3).Edge());
  TColgp_Array1OfPnt poles(1, 3);
  poles.SetValue(1, p0);
  poles.SetValue(2, gp_Pnt(0.5 * radial_limit, 0.0, 0.0));
  poles.SetValue(3, p3);
  Handle(Geom_BezierCurve) parabola = new Geom_BezierCurve(poles);
  TopoDS_Edge parabola_edge = BRepBuilderAPI_MakeEdge(parabola).Edge();
  parabola_edge.Reverse();
  profile_wire.Add(parabola_edge);
  if (!profile_wire.IsDone()) fail("failed to close the paraboloid meridian profile");
  const gp_Pln profile_plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 1.0, 0.0));
  BRepBuilderAPI_MakeFace profile_face(profile_plane, profile_wire.Wire(), Standard_True);
  if (!profile_face.IsDone()) fail("failed to construct the paraboloid meridian face");
  const gp_Ax1 revolution_axis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  BRepPrimAPI_MakeRevol revolution(profile_face.Face(), revolution_axis,
                                    2.0 * std::acos(-1.0), Standard_True);
  if (!revolution.IsDone() || revolution.Shape().IsNull())
    fail("failed to construct the finite paraboloid carrier");

  const double map[3][3] = {
      {transverse_scale_x * quadric.radial_x.X(),
       transverse_scale_y * quadric.radial_y.X(), canonical_axis.X()},
      {transverse_scale_x * quadric.radial_x.Y(),
       transverse_scale_y * quadric.radial_y.Y(), canonical_axis.Y()},
      {transverse_scale_x * quadric.radial_x.Z(),
       transverse_scale_y * quadric.radial_y.Z(), canonical_axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_GTrsf transform(vectorial, gp_XYZ(vertex.X(), vertex.Y(), vertex.Z()));
  BRepBuilderAPI_GTransform affine(revolution.Shape(), transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not place the finite paraboloid carrier in solver coordinates");
  const TopoDS_Shape cavity = affine.Shape();
  if (!BRepCheck_Analyzer(cavity, Standard_True).IsValid())
    fail("elliptic-paraboloid carrier failed OCCT topology validation");

  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("solver paraboloid Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "paraboloid Boolean result");
  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write paraboloid B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "paraboloid B-rep reload");
  validate_roundtrip_solid(step_shape, "paraboloid STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);

  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_ELLIPTIC_PARABOLOID_QUADRIC_STEP_BREP_VALIDATED\","
               "\"paraboloid_vertex\":[" << vertex.X() << ',' << vertex.Y() << ','
            << vertex.Z() << "],\"paraboloid_axis\":[" << canonical_axis.X() << ','
            << canonical_axis.Y() << ',' << canonical_axis.Z()
            << "],\"carrier_radial_limit\":" << radial_limit
            << ",\"carrier_height\":" << cap_height
            << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

TopoDS_Shape checked_one_sheet_affine_carrier(const TopoDS_Shape& analytic,
                                             const gp_GTrsf& transform) {
  // GTransform internally converts conics to NURBS before transforming them.
  // Near a cone that intermediate can already be invalid, while the final
  // shrunken B-rep reports valid and Boolean silently returns the whole domain.
  // Validate the actual intermediate before using it. No healing/tolerance
  // change or approximate cone substitution is an acceptable replacement.
  BRepBuilderAPI_NurbsConvert converted(analytic, Standard_True);
  if (!converted.IsDone() || converted.Shape().IsNull())
    fail("one-sheet hyperboloid NURBS conversion did not complete");
  require_valid_topology(converted.Shape(),"one-sheet hyperboloid NURBS conversion",
                        "one-sheet hyperboloid NURBS conversion failed topology validation");
  BRepBuilderAPI_GTransform affine(converted.Shape(),transform,Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not place the one-sheet hyperboloid carrier");
  return affine.Shape();
}

TopoDS_Shape finite_one_sheet_hyperboloid_carrier(
    const std::array<double, 10>& coefficients, const TopoDS_Shape& domain) {
  const HyperboloidQuadricFrame quadric =
      hyperboloid_from_solver_coefficients(coefficients);
  const double radius_x = std::sqrt(quadric.center_level / quadric.lambda_x);
  const double radius_y = std::sqrt(quadric.center_level / quadric.lambda_y);
  const double axial_scale =
      std::sqrt(quadric.center_level / std::abs(quadric.lambda_axis));
  if (!std::isfinite(radius_x) || !std::isfinite(radius_y) ||
      !std::isfinite(axial_scale) || !(radius_x > 0.0) ||
      !(radius_y > 0.0) || !(axial_scale > 0.0))
    fail("one-sheet hyperboloid canonical scales are invalid");

  Bnd_Box bounds_box;
  BRepBndLib::AddOptimal(domain, bounds_box, Standard_False, Standard_True);
  if (bounds_box.IsVoid() || bounds_box.IsOpenXmin() ||
      bounds_box.IsOpenXmax() || bounds_box.IsOpenYmin() ||
      bounds_box.IsOpenYmax() || bounds_box.IsOpenZmin() ||
      bounds_box.IsOpenZmax())
    fail("hyperboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds_box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double bounds[2][3] = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
  double minimum_axis = std::numeric_limits<double>::infinity();
  double maximum_axis = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(bounds[x][0], bounds[y][1], bounds[z][2]);
        const double projection =
            gp_Vec(quadric.center, corner).Dot(gp_Vec(quadric.axis)) /
            axial_scale;
        minimum_axis = std::min(minimum_axis, projection);
        maximum_axis = std::max(maximum_axis, projection);
      }
  const double diagonal = std::sqrt(
      (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
      (zmax - zmin) * (zmax - zmin));
  const double margin = std::max(1e-7, 0.1 * diagonal) / axial_scale;
  const double z_min = minimum_axis - margin;
  const double z_max = maximum_axis + margin;
  if (!std::isfinite(z_min) || !std::isfinite(z_max) || !(z_max > z_min))
    fail("hyperboloid Boolean has no finite axial interval");
  const double u_min = std::asinh(z_min), u_max = std::asinh(z_max);
  const double r_min = std::cosh(u_min), r_max = std::cosh(u_max);
  if (!std::isfinite(r_min) || !std::isfinite(r_max))
    fail("hyperboloid carrier cap radius is non-finite");

  const gp_Ax2 hyperbola_axes(gp_Pnt(0.0, 0.0, 0.0),
                              gp_Dir(0.0, -1.0, 0.0),
                              gp_Dir(1.0, 0.0, 0.0));
  const gp_Hypr canonical_hyperbola(hyperbola_axes, 1.0, 1.0);
  Handle(Geom_Hyperbola) curve = new Geom_Hyperbola(canonical_hyperbola);
  // Use the curve's own endpoints for adjacent cap edges. At large canonical
  // coordinates sinh(asinh(z)) differs from z by more than wire tolerance.
  const gp_Pnt bottom_surface=curve->Value(u_min), top_surface=curve->Value(u_max);
  if (!(bottom_surface.Z()<minimum_axis && top_surface.Z()>maximum_axis))
    fail("hyperboloid carrier caps no longer lie outside the design domain");
  const gp_Pnt bottom_axis(0.0,0.0,bottom_surface.Z());
  const gp_Pnt top_axis(0.0,0.0,top_surface.Z());
  BRepBuilderAPI_MakeWire profile_wire;
  profile_wire.Add(BRepBuilderAPI_MakeEdge(bottom_axis,top_axis).Edge());
  if (!profile_wire.IsDone()) fail("hyperboloid meridian axis edge was rejected");
  profile_wire.Add(BRepBuilderAPI_MakeEdge(top_axis,top_surface).Edge());
  if (!profile_wire.IsDone()) fail("hyperboloid meridian upper cap edge was rejected");
  TopoDS_Edge hyperbola_edge =
      BRepBuilderAPI_MakeEdge(curve, u_min, u_max).Edge();
  hyperbola_edge.Reverse();
  profile_wire.Add(hyperbola_edge);
  if (!profile_wire.IsDone()) fail("hyperboloid meridian curved edge was rejected");
  profile_wire.Add(BRepBuilderAPI_MakeEdge(bottom_surface, bottom_axis).Edge());
  if (!profile_wire.IsDone() || !BRep_Tool::IsClosed(profile_wire.Wire()))
    fail("failed to close the one-sheet hyperboloid meridian profile");
  const gp_Pln profile_plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 1.0, 0.0));
  BRepBuilderAPI_MakeFace profile_face(profile_plane, profile_wire.Wire(),
                                        Standard_True);
  if (!profile_face.IsDone())
    fail("failed to construct the hyperboloid meridian face");
  trace_carrier_topology(profile_face.Face(),"meridian");
  const gp_Ax1 revolution_axis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  BRepPrimAPI_MakeRevol revolution(profile_face.Face(), revolution_axis,
                                    2.0 * std::acos(-1.0), Standard_True);
  if (!revolution.IsDone() || revolution.Shape().IsNull())
    fail("failed to construct the finite one-sheet hyperboloid carrier");
  trace_carrier_topology(revolution.Shape(),"revolution");

  const double map[3][3] = {
      {radius_x * quadric.radial_x.X(), radius_y * quadric.radial_y.X(),
       axial_scale * quadric.axis.X()},
      {radius_x * quadric.radial_x.Y(), radius_y * quadric.radial_y.Y(),
       axial_scale * quadric.axis.Y()},
      {radius_x * quadric.radial_x.Z(), radius_y * quadric.radial_y.Z(),
       axial_scale * quadric.axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2], map[1][0],
                         map[1][1], map[1][2], map[2][0], map[2][1],
                         map[2][2]);
  const gp_GTrsf transform(
      vectorial, gp_XYZ(quadric.center.X(), quadric.center.Y(),
                        quadric.center.Z()));
  const TopoDS_Shape carrier = checked_one_sheet_affine_carrier(revolution.Shape(),transform);
  trace_carrier_topology(carrier,"affine");
  require_valid_topology(carrier,"one-sheet hyperboloid carrier",
                        "one-sheet hyperboloid carrier failed OCCT topology validation");
  return carrier;
}

TopoDS_Shape finite_two_sheet_hyperboloid_carrier(
    const std::array<double, 10>& negative_inside_coefficients,
    const TopoDS_Shape& domain) {
  const HyperboloidQuadricFrame quadric =
      hyperboloid_from_solver_coefficients(negative_inside_coefficients, true);
  const double level = -quadric.center_level;
  const double radial_scale_x = std::sqrt(level / quadric.lambda_x);
  const double radial_scale_y = std::sqrt(level / quadric.lambda_y);
  const double axial_scale = std::sqrt(level / std::abs(quadric.lambda_axis));
  if (!std::isfinite(radial_scale_x) || !std::isfinite(radial_scale_y) ||
      !std::isfinite(axial_scale) || !(radial_scale_x > 0.0) ||
      !(radial_scale_y > 0.0) || !(axial_scale > 0.0))
    fail("two-sheet hyperboloid canonical scales are invalid");

  Bnd_Box bounds_box;
  BRepBndLib::AddOptimal(domain, bounds_box, Standard_False, Standard_True);
  if (bounds_box.IsVoid() || bounds_box.IsOpenXmin() ||
      bounds_box.IsOpenXmax() || bounds_box.IsOpenYmin() ||
      bounds_box.IsOpenYmax() || bounds_box.IsOpenZmin() ||
      bounds_box.IsOpenZmax())
    fail("two-sheet hyperboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds_box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double bounds[2][3] = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
  double minimum_axis = std::numeric_limits<double>::infinity();
  double maximum_axis = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(bounds[x][0], bounds[y][1], bounds[z][2]);
        const double projection =
            gp_Vec(quadric.center, corner).Dot(gp_Vec(quadric.axis)) /
            axial_scale;
        minimum_axis = std::min(minimum_axis, projection);
        maximum_axis = std::max(maximum_axis, projection);
      }
  const double diagonal = std::sqrt(
      (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
      (zmax - zmin) * (zmax - zmin));
  const double margin = std::max(1e-7, 0.1 * diagonal) / axial_scale;

  const gp_Ax2 hyperbola_axes(gp_Pnt(0.0, 0.0, 0.0),
                              gp_Dir(0.0, 0.0, 1.0),
                              gp_Dir(1.0, 0.0, 0.0));
  const gp_Hypr canonical_hyperbola(hyperbola_axes, 1.0, 1.0);
  Handle(Geom_Hyperbola) curve = new Geom_Hyperbola(canonical_hyperbola);
  const gp_Pln profile_plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const gp_Ax1 revolution_axis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(1.0, 0.0, 0.0));
  const auto build_lobe = [&](double cap) {
    if (!(cap > 1.0) || !std::isfinite(cap))
      fail("two-sheet hyperboloid lobe cap does not lie beyond its apex");
    const double parameter = std::acosh(cap);
    const double radius = std::sinh(parameter);
    if (!std::isfinite(parameter) || !std::isfinite(radius) || !(radius > 0.0))
      fail("two-sheet hyperboloid carrier profile is non-finite");
    const gp_Pnt tip(1.0, 0.0, 0.0);
    const gp_Pnt cap_axis(cap, 0.0, 0.0);
    const gp_Pnt cap_surface(cap, radius, 0.0);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(tip, cap_axis).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(cap_axis, cap_surface).Edge());
    TopoDS_Edge boundary = BRepBuilderAPI_MakeEdge(curve, 0.0, parameter).Edge();
    boundary.Reverse();
    wire.Add(boundary);
    if (!wire.IsDone())
      fail("failed to close a two-sheet hyperboloid meridian lobe");
    BRepBuilderAPI_MakeFace face(profile_plane, wire.Wire(), Standard_True);
    if (!face.IsDone())
      fail("failed to construct a two-sheet hyperboloid lobe face");
    BRepPrimAPI_MakeRevol revolution(face.Face(), revolution_axis,
                                      2.0 * std::acos(-1.0), Standard_True);
    if (!revolution.IsDone() || revolution.Shape().IsNull())
      fail("failed to revolve a two-sheet hyperboloid lobe");
    return revolution.Shape();
  };

  BRep_Builder builder;
  TopoDS_Compound lobes;
  builder.MakeCompound(lobes);
  int lobe_count = 0;
  const auto add_lobe = [&](double sign, double cap) {
    if (!(cap > 1.0)) return;
    const TopoDS_Shape local = build_lobe(std::max(1.0 + 1e-7, cap + margin));
    const gp_Vec axis(quadric.axis);
    const gp_Vec radial_x(quadric.radial_x);
    const gp_Vec radial_y(quadric.radial_y);
    const gp_Vec mapped_axis = axis * sign;
    const gp_Vec mapped_radial_y = radial_y * sign;
    const double map[3][3] = {
        {axial_scale * mapped_axis.X(), radial_scale_x * radial_x.X(),
         radial_scale_y * mapped_radial_y.X()},
        {axial_scale * mapped_axis.Y(), radial_scale_x * radial_x.Y(),
         radial_scale_y * mapped_radial_y.Y()},
        {axial_scale * mapped_axis.Z(), radial_scale_x * radial_x.Z(),
         radial_scale_y * mapped_radial_y.Z()}};
    const gp_Mat matrix(map[0][0], map[0][1], map[0][2], map[1][0],
                        map[1][1], map[1][2], map[2][0], map[2][1],
                        map[2][2]);
    const gp_GTrsf transform(matrix, gp_XYZ(quadric.center.X(),
                                            quadric.center.Y(),
                                            quadric.center.Z()));
    BRepBuilderAPI_GTransform affine(local, transform, Standard_True);
    if (!affine.IsDone() || affine.Shape().IsNull())
      fail("OCCT could not place a two-sheet hyperboloid lobe");
    const TopoDS_Shape carrier = affine.Shape();
    if (!BRepCheck_Analyzer(carrier, Standard_True).IsValid())
      fail("two-sheet hyperboloid lobe carrier failed topology validation");
    builder.Add(lobes, carrier);
    ++lobe_count;
  };
  add_lobe(+1.0, maximum_axis - margin);
  add_lobe(-1.0, -minimum_axis - margin);
  if (lobe_count == 0)
    fail("two-sheet hyperboloid has no lobe intersecting the design-domain bounds");
  return lobes;
}

int run_solver_one_sheet_hyperboloid_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-hyperboloid-one-sheet INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  for (int index = 0; index < 10; ++index)
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
  const HyperboloidQuadricFrame quadric =
      hyperboloid_from_solver_coefficients(coefficients);
  const double radius_x = std::sqrt(quadric.center_level / quadric.lambda_x);
  const double radius_y = std::sqrt(quadric.center_level / quadric.lambda_y);
  const double axial_scale =
      std::sqrt(quadric.center_level / std::abs(quadric.lambda_axis));
  if (!std::isfinite(radius_x) || !std::isfinite(radius_y) ||
      !std::isfinite(axial_scale) || !(radius_x > 0.0) ||
      !(radius_y > 0.0) || !(axial_scale > 0.0))
    fail("one-sheet hyperboloid canonical scales are invalid");
  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);

  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  Bnd_Box outer_bounds;
  BRepBndLib::AddOptimal(outer, outer_bounds, Standard_False, Standard_True);
  if (outer_bounds.IsVoid() || outer_bounds.IsOpenXmin() ||
      outer_bounds.IsOpenXmax() || outer_bounds.IsOpenYmin() ||
      outer_bounds.IsOpenYmax() || outer_bounds.IsOpenZmin() ||
      outer_bounds.IsOpenZmax())
    fail("hyperboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  outer_bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double bounds[2][3] = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
  double minimum_axis = std::numeric_limits<double>::infinity();
  double maximum_axis = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(bounds[x][0], bounds[y][1], bounds[z][2]);
        const double projection =
            gp_Vec(quadric.center, corner).Dot(gp_Vec(quadric.axis)) /
            axial_scale;
        minimum_axis = std::min(minimum_axis, projection);
        maximum_axis = std::max(maximum_axis, projection);
      }
  const double domain_diagonal = std::sqrt(
      (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
      (zmax - zmin) * (zmax - zmin));
  const double parameter_margin = std::max(1e-7, 0.1 * domain_diagonal) /
                                  axial_scale;
  const double z_min = minimum_axis - parameter_margin;
  const double z_max = maximum_axis + parameter_margin;
  if (!std::isfinite(z_min) || !std::isfinite(z_max) || !(z_max > z_min))
    fail("hyperboloid Boolean has no finite axial interval");
  const double u_min = std::asinh(z_min);
  const double u_max = std::asinh(z_max);
  const double r_min = std::cosh(u_min);
  const double r_max = std::cosh(u_max);
  if (!std::isfinite(r_min) || !std::isfinite(r_max))
    fail("hyperboloid carrier cap radius is non-finite");

  // Exact meridian section X=cosh(u), Z=sinh(u) of the canonical
  // one-sheet sublevel X^2+Y^2-Z^2 <= 1. Its planar caps lie strictly beyond
  // the design-domain projection, so truncation cannot alter the Boolean.
  const gp_Ax2 hyperbola_axes(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, -1.0, 0.0),
                             gp_Dir(1.0, 0.0, 0.0));
  const gp_Hypr canonical_hyperbola(hyperbola_axes, 1.0, 1.0);
  Handle(Geom_Hyperbola) curve = new Geom_Hyperbola(canonical_hyperbola);
  const gp_Pnt bottom_surface=curve->Value(u_min), top_surface=curve->Value(u_max);
  if (!(bottom_surface.Z()<minimum_axis && top_surface.Z()>maximum_axis))
    fail("hyperboloid carrier caps no longer lie outside the design domain");
  const gp_Pnt bottom_axis(0.0,0.0,bottom_surface.Z());
  const gp_Pnt top_axis(0.0,0.0,top_surface.Z());
  BRepBuilderAPI_MakeWire profile_wire;
  profile_wire.Add(BRepBuilderAPI_MakeEdge(bottom_axis,top_axis).Edge());
  if (!profile_wire.IsDone()) fail("hyperboloid meridian axis edge was rejected");
  profile_wire.Add(BRepBuilderAPI_MakeEdge(top_axis,top_surface).Edge());
  if (!profile_wire.IsDone()) fail("hyperboloid meridian upper cap edge was rejected");
  TopoDS_Edge hyperbola_edge =
      BRepBuilderAPI_MakeEdge(curve, u_min, u_max).Edge();
  hyperbola_edge.Reverse();
  profile_wire.Add(hyperbola_edge);
  if (!profile_wire.IsDone()) fail("hyperboloid meridian curved edge was rejected");
  profile_wire.Add(BRepBuilderAPI_MakeEdge(bottom_surface, bottom_axis).Edge());
  if (!profile_wire.IsDone() || !BRep_Tool::IsClosed(profile_wire.Wire()))
    fail("failed to close the one-sheet hyperboloid meridian profile");
  const gp_Pln profile_plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 1.0, 0.0));
  BRepBuilderAPI_MakeFace profile_face(profile_plane, profile_wire.Wire(),
                                        Standard_True);
  if (!profile_face.IsDone())
    fail("failed to construct the hyperboloid meridian face");
  trace_carrier_topology(profile_face.Face(),"meridian");
  const gp_Ax1 revolution_axis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  BRepPrimAPI_MakeRevol revolution(profile_face.Face(), revolution_axis,
                                    2.0 * std::acos(-1.0), Standard_True);
  if (!revolution.IsDone() || revolution.Shape().IsNull())
    fail("failed to construct the finite one-sheet hyperboloid carrier");
  trace_carrier_topology(revolution.Shape(),"revolution");

  const double map[3][3] = {
      {radius_x * quadric.radial_x.X(), radius_y * quadric.radial_y.X(),
       axial_scale * quadric.axis.X()},
      {radius_x * quadric.radial_x.Y(), radius_y * quadric.radial_y.Y(),
       axial_scale * quadric.axis.Y()},
      {radius_x * quadric.radial_x.Z(), radius_y * quadric.radial_y.Z(),
       axial_scale * quadric.axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2], map[1][0],
                         map[1][1], map[1][2], map[2][0], map[2][1],
                         map[2][2]);
  const gp_GTrsf transform(
      vectorial, gp_XYZ(quadric.center.X(), quadric.center.Y(),
                        quadric.center.Z()));
  const TopoDS_Shape cavity = checked_one_sheet_affine_carrier(revolution.Shape(),transform);
  trace_carrier_topology(cavity,"affine");
  require_valid_topology(cavity,"one-sheet hyperboloid carrier",
                        "one-sheet hyperboloid carrier failed OCCT topology validation");
  BRepAlgoAPI_Cut cut(outer, cavity);
  cut.SetRunParallel(Standard_False);
  cut.Build();
  if (!cut.IsDone() || cut.Shape().IsNull())
    fail("one-sheet hyperboloid Boolean cut did not complete");
  const TopoDS_Shape material = cut.Shape();
  validate_roundtrip_solid(material, "one-sheet hyperboloid Boolean result");
  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write one-sheet hyperboloid B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "one-sheet hyperboloid B-rep reload");
  validate_roundtrip_solid(step_shape, "one-sheet hyperboloid STEP reload");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);

  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_ONE_SHEET_HYPERBOLOID_QUADRIC_STEP_BREP_VALIDATED\","
               "\"hyperboloid_center\":[" << quadric.center.X() << ','
            << quadric.center.Y() << ',' << quadric.center.Z()
            << "],\"hyperboloid_axis\":[" << quadric.axis.X() << ','
            << quadric.axis.Y() << ',' << quadric.axis.Z()
            << "],\"canonical_scales\":[" << radius_x << ',' << radius_y
            << ',' << axial_scale << "],\"carrier_axis_interval\":["
            << z_min << ',' << z_max << "],\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

int run_solver_two_sheet_hyperboloid_quadric(int argc, char** argv) {
  if (argc != 16)
    fail("usage: occt_step_roundtrip --quadric-hyperboloid-two-sheet INPUT.obj COORDINATE_SCALE c bx by bz qxx qyy qzz qyz qxz qxy MATERIAL.brep MATERIAL.step");
  const std::string obj_path(argv[2]);
  const double coordinate_scale = parse_coordinate_scale(argv[3]);
  std::array<double, 10> coefficients{};
  std::array<double, 10> complemented_coefficients{};
  for (int index = 0; index < 10; ++index) {
    coefficients[index] = parse_finite_double(argv[4 + index], "quadric coefficient");
    complemented_coefficients[index] = -coefficients[index];
  }
  // For the (+,-,-) solver quadric, retained material is q>=0. Negating its
  // coefficients gives two analytic q'<=0 lobes with the (+,+,-) parser's
  // radial eigendirections and a negative center level.
  const HyperboloidQuadricFrame quadric =
      hyperboloid_from_solver_coefficients(complemented_coefficients, true);
  const double level = -quadric.center_level;
  const double radial_scale_x = std::sqrt(level / quadric.lambda_x);
  const double radial_scale_y = std::sqrt(level / quadric.lambda_y);
  const double axial_scale = std::sqrt(level / std::abs(quadric.lambda_axis));
  if (!std::isfinite(radial_scale_x) || !std::isfinite(radial_scale_y) ||
      !std::isfinite(axial_scale) || !(radial_scale_x > 0.0) ||
      !(radial_scale_y > 0.0) || !(axial_scale > 0.0))
    fail("two-sheet hyperboloid canonical scales are invalid");

  const std::string brep_path(argv[14]);
  const std::string step_path(argv[15]);
  Message::DefaultMessenger()->ChangePrinters().Clear();
  const TopoDS_Shape outer = mesh_to_solid(obj_path, coordinate_scale);
  Bnd_Box outer_bounds;
  BRepBndLib::AddOptimal(outer, outer_bounds, Standard_False, Standard_True);
  if (outer_bounds.IsVoid() || outer_bounds.IsOpenXmin() ||
      outer_bounds.IsOpenXmax() || outer_bounds.IsOpenYmin() ||
      outer_bounds.IsOpenYmax() || outer_bounds.IsOpenZmin() ||
      outer_bounds.IsOpenZmax())
    fail("two-sheet hyperboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  outer_bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const double bounds[2][3] = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
  double minimum_axis = std::numeric_limits<double>::infinity();
  double maximum_axis = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(bounds[x][0], bounds[y][1], bounds[z][2]);
        const double projection =
            gp_Vec(quadric.center, corner).Dot(gp_Vec(quadric.axis)) /
            axial_scale;
        minimum_axis = std::min(minimum_axis, projection);
        maximum_axis = std::max(maximum_axis, projection);
      }
  const double domain_diagonal = std::sqrt(
      (xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) +
      (zmax - zmin) * (zmax - zmin));
  const double margin = std::max(1e-7, 0.1 * domain_diagonal) / axial_scale;
  const gp_Ax2 hyperbola_axes(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0),
                             gp_Dir(1.0, 0.0, 0.0));
  const gp_Hypr canonical_hyperbola(hyperbola_axes, 1.0, 1.0);
  Handle(Geom_Hyperbola) curve = new Geom_Hyperbola(canonical_hyperbola);
  const gp_Pln profile_plane(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
  const gp_Ax1 revolution_axis(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(1.0, 0.0, 0.0));
  const auto build_lobe = [&](double cap) {
    if (!(cap > 1.0) || !std::isfinite(cap))
      fail("two-sheet hyperboloid lobe cap does not lie beyond its apex");
    const double parameter = std::acosh(cap);
    const double radius = std::sinh(parameter);
    if (!std::isfinite(parameter) || !std::isfinite(radius) || !(radius > 0.0))
      fail("two-sheet hyperboloid carrier profile is non-finite");
    const gp_Pnt tip(1.0, 0.0, 0.0);
    const gp_Pnt cap_axis(cap, 0.0, 0.0);
    const gp_Pnt cap_surface(cap, radius, 0.0);
    BRepBuilderAPI_MakeWire wire;
    wire.Add(BRepBuilderAPI_MakeEdge(tip, cap_axis).Edge());
    wire.Add(BRepBuilderAPI_MakeEdge(cap_axis, cap_surface).Edge());
    TopoDS_Edge boundary = BRepBuilderAPI_MakeEdge(curve, 0.0, parameter).Edge();
    boundary.Reverse();
    wire.Add(boundary);
    if (!wire.IsDone())
      fail("failed to close a two-sheet hyperboloid meridian lobe");
    BRepBuilderAPI_MakeFace face(profile_plane, wire.Wire(), Standard_True);
    if (!face.IsDone())
      fail("failed to construct a two-sheet hyperboloid lobe face");
    BRepPrimAPI_MakeRevol revolution(face.Face(), revolution_axis,
                                      2.0 * std::acos(-1.0), Standard_True);
    if (!revolution.IsDone() || revolution.Shape().IsNull())
      fail("failed to revolve a two-sheet hyperboloid lobe");
    return revolution.Shape();
  };
  BRep_Builder compound_builder;
  TopoDS_Compound material_compound;
  compound_builder.MakeCompound(material_compound);
  int material_solid_count = 0;
  const auto intersect_lobe = [&](double axial_sign, double cap) {
    if (!(cap > 1.0)) return;
    const TopoDS_Shape local_lobe = build_lobe(cap);
    const gp_Vec radial_x(quadric.radial_x);
    const gp_Vec radial_y(quadric.radial_y);
    const gp_Vec axis(quadric.axis);
    // Local axes are (axial, radial-x, radial-y). A paired sign flip of the
    // axial and radial-y columns keeps the affine map right-handed.
    const double sign = axial_sign;
    const gp_Vec mapped_axis = axis * sign;
    const gp_Vec mapped_radial_y = radial_y * sign;
    const double map[3][3] = {
        {axial_scale * mapped_axis.X(), radial_scale_x * radial_x.X(),
         radial_scale_y * mapped_radial_y.X()},
        {axial_scale * mapped_axis.Y(), radial_scale_x * radial_x.Y(),
         radial_scale_y * mapped_radial_y.Y()},
        {axial_scale * mapped_axis.Z(), radial_scale_x * radial_x.Z(),
         radial_scale_y * mapped_radial_y.Z()}};
    const gp_Mat matrix(map[0][0], map[0][1], map[0][2], map[1][0],
                        map[1][1], map[1][2], map[2][0], map[2][1],
                        map[2][2]);
    const gp_GTrsf transform(matrix, gp_XYZ(quadric.center.X(),
                                           quadric.center.Y(),
                                           quadric.center.Z()));
    BRepBuilderAPI_GTransform affine(local_lobe, transform, Standard_True);
    if (!affine.IsDone() || affine.Shape().IsNull())
      fail("OCCT could not place a two-sheet hyperboloid lobe");
    const TopoDS_Shape carrier = affine.Shape();
    if (!BRepCheck_Analyzer(carrier, Standard_True).IsValid())
      fail("two-sheet hyperboloid lobe carrier failed topology validation");
    BRepAlgoAPI_Common common(outer, carrier);
    common.SetRunParallel(Standard_False);
    common.SetFuzzyValue(1e-9);
    common.Build();
    if (!common.IsDone())
      fail("two-sheet hyperboloid/domain intersection did not complete");
    if (common.Shape().IsNull()) return;
    int lobe_solids = 0;
    for (TopExp_Explorer it(common.Shape(), TopAbs_SOLID); it.More(); it.Next()) {
      const TopoDS_Shape solid_shape = it.Current();
      if (!(volume_moments(solid_shape).volume > 0.0))
        fail("two-sheet hyperboloid Boolean produced a non-positive material lobe");
      compound_builder.Add(material_compound, solid_shape);
      ++lobe_solids;
    }
    material_solid_count += lobe_solids;
  };
  const double positive_cap = maximum_axis - margin;
  const double negative_cap = -minimum_axis - margin;
  if (positive_cap > 1.0)
    intersect_lobe(+1.0, std::max(1.0 + 1e-7, positive_cap + 2.0 * margin));
  if (negative_cap > 1.0)
    intersect_lobe(-1.0, std::max(1.0 + 1e-7, negative_cap + 2.0 * margin));
  if (material_solid_count == 0)
    fail("two-sheet hyperboloid has no positive-volume retained material in the design domain");
  validate_roundtrip_solid(material_compound, "two-sheet hyperboloid material compound");
  if (!BRepTools::Write(material_compound, brep_path.c_str()))
    fail("failed to write two-sheet hyperboloid B-rep compound");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material_compound, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "two-sheet hyperboloid B-rep reload");
  validate_roundtrip_solid(step_shape, "two-sheet hyperboloid STEP reload");
  const auto count_solids = [](const TopoDS_Shape& shape) {
    int count = 0;
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next()) ++count;
    return count;
  };
  if (count_solids(brep_shape) != material_solid_count ||
      count_solids(step_shape) != material_solid_count)
    fail("two-sheet B-rep/STEP round-trip changed the material solid count");
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);
  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD\","
               "\"status\":\"ANALYTIC_TWO_SHEET_HYPERBOLOID_QUADRIC_STEP_BREP_VALIDATED\","
               "\"hyperboloid_center\":[" << quadric.center.X() << ','
            << quadric.center.Y() << ',' << quadric.center.Z()
            << "],\"hyperboloid_axis\":[" << quadric.axis.X() << ','
            << quadric.axis.Y() << ',' << quadric.axis.Z()
            << "],\"canonical_scales\":[" << radial_scale_x << ','
            << radial_scale_y << ',' << axial_scale
            << "],\"material_solid_count\":" << material_solid_count
            << ",\"raw_moments\":{\"brep\":";
  emit_raw_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_raw_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}

struct PlanarCadPiece {
  char kind = 'L';
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0;
};

struct PlanarCadEntity {
  std::string role;
  std::string identity;
  double density = 0.0;
  double z_min = 0.0, z_max = 0.0;
  std::vector<std::vector<PlanarCadPiece>> cycles;
};

std::string read_token(std::istream& input, const std::string& label) {
  std::string value;
  if (!(input >> value)) fail("truncated layered-extrusion input while reading " + label);
  return value;
}

std::size_t read_count(std::istream& input, const std::string& label,
                       std::size_t minimum, std::size_t maximum) {
  const std::string token = read_token(input, label);
  errno = 0;
  char* end = nullptr;
  const unsigned long long value = std::strtoull(token.c_str(), &end, 10);
  if (end == token.c_str() || *end != '\0' || errno == ERANGE ||
      value < minimum || value > maximum)
    fail(label + " is outside its accepted integer range");
  return static_cast<std::size_t>(value);
}

double read_real(std::istream& input, const std::string& label) {
  return parse_finite_double(read_token(input, label), label);
}

double signed_cycle_area(const std::vector<PlanarCadPiece>& cycle) {
  const double pi = std::acos(-1.0);
  double twice_area = 0.0;
  for (const PlanarCadPiece& piece : cycle) {
    if (piece.kind == 'L') {
      twice_area += piece.a * piece.d - piece.b * piece.c;
    } else {
      const double start = piece.d;
      const double finish = start + piece.e;
      twice_area += piece.c * piece.c * piece.e +
          piece.c * piece.a * (std::sin(finish) - std::sin(start)) -
          piece.c * piece.b * (std::cos(finish) - std::cos(start));
    }
  }
  (void)pi;
  return 0.5 * twice_area;
}

TopoDS_Wire make_planar_wire(const std::vector<PlanarCadPiece>& cycle,
                             double z) {
  BRepBuilderAPI_MakeWire wire_builder;
  for (const PlanarCadPiece& piece : cycle) {
    TopoDS_Edge edge;
    if (piece.kind == 'L') {
      BRepBuilderAPI_MakeEdge maker(gp_Pnt(piece.a, piece.b, z),
                                    gp_Pnt(piece.c, piece.d, z));
      if (!maker.IsDone()) fail("failed to construct a planar CAD line edge");
      edge = maker.Edge();
    } else {
      const gp_Circ circle(gp_Ax2(gp_Pnt(piece.a, piece.b, z),
                                  gp_Dir(0.0, 0.0, 1.0)), piece.c);
      const double full_turn = 2.0 * std::acos(-1.0);
      if (std::abs(piece.e) >= full_turn - 1e-12) {
        BRepBuilderAPI_MakeEdge maker(circle, 0.0, full_turn);
        if (!maker.IsDone()) fail("failed to construct a full-circle CAD edge");
        edge = maker.Edge();
        if (piece.e < 0.0) edge.Reverse();
      } else {
        double start = std::fmod(piece.d, full_turn);
        if (start < 0.0) start += full_turn;
        double finish = start + piece.e;
        if (piece.e > 0.0 && finish > full_turn) {
          start -= full_turn;
          finish -= full_turn;
        } else if (piece.e < 0.0 && finish < 0.0) {
          start += full_turn;
          finish += full_turn;
        }
        const double lower = std::min(start, finish);
        const double upper = std::max(start, finish);
        BRepBuilderAPI_MakeEdge maker(circle, lower, upper);
        if (!maker.IsDone()) fail("failed to construct a circular CAD arc");
        edge = maker.Edge();
        if (piece.e < 0.0) edge.Reverse();
      }
    }
    wire_builder.Add(edge);
    if (!wire_builder.IsDone()) fail("planar CAD boundary contains disconnected edges");
  }
  const TopoDS_Wire wire = wire_builder.Wire();
  if (!BRep_Tool::IsClosed(wire) || !BRepCheck_Analyzer(wire, Standard_True).IsValid())
    fail("planar CAD boundary cycle is not a valid closed wire");
  return wire;
}

std::vector<TopoDS_Face> make_planar_faces(const PlanarCadEntity& entity) {
  const gp_Pln plane(gp_Pnt(0.0, 0.0, entity.z_min), gp_Dir(0.0, 0.0, 1.0));
  struct CycleWire { TopoDS_Wire wire; double area; gp_Pnt witness; };
  std::vector<CycleWire> outers, holes;
  for (const auto& cycle : entity.cycles) {
    const double area = signed_cycle_area(cycle);
    if (!std::isfinite(area) || std::abs(area) <= 1e-14)
      fail("planar CAD cycle has zero or non-finite signed area");
    const TopoDS_Wire wire = make_planar_wire(cycle, entity.z_min);
    const PlanarCadPiece& first = cycle.front();
    const gp_Pnt witness(first.a, first.b, entity.z_min);
    (area > 0.0 ? outers : holes).push_back({wire, area, witness});
  }
  if (outers.empty()) fail("planar CAD entity has no positively oriented outer cycle");

  std::vector<std::vector<TopoDS_Wire>> assigned_holes(outers.size());
  for (const CycleWire& hole : holes) {
    std::size_t owner = outers.size();
    for (std::size_t index = 0; index < outers.size(); ++index) {
      BRepBuilderAPI_MakeFace outer_face(plane, outers[index].wire, Standard_True);
      if (!outer_face.IsDone()) fail("failed to build planar CAD outer face");
      BRepClass_FaceClassifier classifier(outer_face.Face(), hole.witness,
                                          1e-9);
      if (classifier.State() == TopAbs_IN) {
        if (owner != outers.size())
          fail("planar CAD hole is ambiguously contained in multiple outer cycles");
        owner = index;
      }
    }
    if (owner == outers.size())
      fail("planar CAD hole is not contained in an outer material cycle");
    assigned_holes[owner].push_back(hole.wire);
  }

  std::vector<TopoDS_Face> faces;
  faces.reserve(outers.size());
  for (std::size_t index = 0; index < outers.size(); ++index) {
    BRepBuilderAPI_MakeFace face_builder(plane, outers[index].wire, Standard_True);
    if (!face_builder.IsDone()) fail("failed to build planar CAD outer face");
    for (const TopoDS_Wire& hole : assigned_holes[index]) face_builder.Add(hole);
    if (!face_builder.IsDone()) fail("failed to add a planar CAD hole wire");
    const TopoDS_Face face = face_builder.Face();
    if (!BRepCheck_Analyzer(face, Standard_True).IsValid())
      fail("planar CAD face failed OCCT topology validation");
    faces.push_back(face);
  }
  return faces;
}

std::vector<PlanarCadEntity> read_planar_cad_input(const std::string& path) {
  std::ifstream input(path);
  if (!input) fail("cannot open layered-extrusion input");
  if (read_token(input, "schema") != "PLANAR_LAYERED_EXTRUSION_V2")
    fail("unsupported layered-extrusion input schema");
  const std::size_t entity_count = read_count(input, "entity count", 1, 10000);
  std::vector<PlanarCadEntity> entities;
  entities.reserve(entity_count);
  for (std::size_t entity_index = 0; entity_index < entity_count; ++entity_index) {
    if (read_token(input, "entity marker") != "ENTITY")
      fail("missing layered-extrusion entity marker");
    PlanarCadEntity entity;
    entity.role = read_token(input, "entity role");
    if ((entity_index == 0 && entity.role != "editable_design") ||
        (entity_index > 0 && entity.role != "fixed_material"))
      fail("layered-extrusion entity order/role is invalid");
    entity.identity = read_token(input, "entity geometry/density identity");
    if (entity.identity.size() != 64 ||
        entity.identity.find_first_not_of("0123456789abcdef") != std::string::npos)
      fail("layered-extrusion entity identity must be a lowercase SHA-256 token");
    entity.density = read_real(input, "entity density");
    entity.z_min = read_real(input, "entity lower z");
    entity.z_max = read_real(input, "entity upper z");
    if (!(entity.density > 0.0) || !(entity.z_max > entity.z_min))
      fail("layered-extrusion density or z interval is invalid");
    const std::size_t cycle_count = read_count(input, "cycle count", 1, 100000);
    entity.cycles.reserve(cycle_count);
    for (std::size_t cycle_index = 0; cycle_index < cycle_count; ++cycle_index) {
      const std::size_t piece_count = read_count(input, "piece count", 1, 1000000);
      std::vector<PlanarCadPiece> cycle;
      cycle.reserve(piece_count);
      for (std::size_t piece_index = 0; piece_index < piece_count; ++piece_index) {
        const std::string type = read_token(input, "boundary piece type");
        PlanarCadPiece piece;
        if (type == "L") {
          piece.kind = 'L';
          piece.a = read_real(input, "line start x");
          piece.b = read_real(input, "line start y");
          piece.c = read_real(input, "line end x");
          piece.d = read_real(input, "line end y");
        } else if (type == "A") {
          piece.kind = 'A';
          piece.a = read_real(input, "arc center x");
          piece.b = read_real(input, "arc center y");
          piece.c = read_real(input, "arc radius");
          piece.d = read_real(input, "arc start angle");
          piece.e = read_real(input, "arc sweep");
          if (!(piece.c > 0.0) || piece.e == 0.0 ||
              std::abs(piece.e) > 2.0 * std::acos(-1.0) + 1e-12)
            fail("circular CAD arc has invalid radius or sweep");
        } else {
          fail("unknown layered-extrusion boundary piece type");
        }
        cycle.push_back(piece);
      }
      entity.cycles.push_back(std::move(cycle));
    }
    entities.push_back(std::move(entity));
  }
  std::string trailing;
  if (input >> trailing) fail("layered-extrusion input has trailing tokens");
  return entities;
}

Moments weighted_moments(const Moments& geometric, double density) {
  Moments result = geometric;
  for (double& value : result.raw) value *= density;
  result.volume *= density;
  result.x *= density;
  result.y *= density;
  result.z *= density;
  return result;
}

int run_planar_layered_extrusion_roundtrip(int argc, char** argv) {
  if (argc != 5)
    fail("usage: occt_step_roundtrip --planar-layered-extrusion INPUT.txt MATERIAL.brep MATERIAL.step");
  const std::vector<PlanarCadEntity> entities = read_planar_cad_input(argv[2]);
  const std::string brep_path(argv[3]), step_path(argv[4]);
  if (brep_path == step_path) fail("B-rep and STEP output paths must differ");
  Message::DefaultMessenger()->ChangePrinters().Clear();
  BRep_Builder compound_builder;
  TopoDS_Compound material;
  compound_builder.MakeCompound(material);
  std::vector<Moments> source_solids;
  std::vector<double> solid_densities;
  std::vector<Moments> entity_mass_moments;
  std::vector<std::size_t> entity_solid_counts;
  for (const PlanarCadEntity& entity : entities) {
    const std::vector<TopoDS_Face> faces = make_planar_faces(entity);
    BRep_Builder entity_builder;
    TopoDS_Compound entity_shape;
    entity_builder.MakeCompound(entity_shape);
    std::size_t solid_count = 0;
    for (const TopoDS_Face& face : faces) {
      BRepPrimAPI_MakePrism prism(face, gp_Vec(0.0, 0.0, entity.z_max - entity.z_min),
                                  Standard_True);
      if (!prism.IsDone() || prism.Shape().IsNull())
        fail("exact planar profile extrusion failed");
      if (!BRepCheck_Analyzer(prism.Shape(), Standard_True).IsValid())
        fail("extruded planar material failed OCCT topology validation");
      bool found_solid = false;
      for (TopExp_Explorer it(prism.Shape(), TopAbs_SOLID); it.More(); it.Next()) {
        const TopoDS_Shape solid = it.Current();
        const Moments moments = volume_moments(solid);
        if (!(moments.volume > 0.0)) fail("extruded material solid has non-positive volume");
        compound_builder.Add(material, solid);
        entity_builder.Add(entity_shape, solid);
        source_solids.push_back(moments);
        solid_densities.push_back(entity.density);
        ++solid_count;
        found_solid = true;
      }
      if (!found_solid) fail("exact planar extrusion produced no solid");
    }
    entity_solid_counts.push_back(solid_count);
    entity_mass_moments.push_back(weighted_moments(volume_moments(entity_shape), entity.density));
  }
  if (source_solids.empty()) fail("layered planar result produced no material solids");
  validate_roundtrip_solid(material, "layered planar material compound");
  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write layered planar B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "layered planar B-rep reload");
  validate_roundtrip_solid(step_shape, "layered planar STEP reload");
  const auto solids_in = [](const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> result;
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next())
      result.push_back(it.Current());
    return result;
  };
  const std::vector<TopoDS_Shape> brep_solids = solids_in(brep_shape);
  const std::vector<TopoDS_Shape> step_solids = solids_in(step_shape);
  if (brep_solids.size() != source_solids.size() || step_solids.size() != source_solids.size())
    fail("B-rep/STEP round-trip changed layered material solid count");
  const auto integrate_reloaded = [&](const std::vector<TopoDS_Shape>& solids,
                                      const std::string& stage) {
    Moments total{};
    for (std::size_t index = 0; index < solids.size(); ++index) {
      const Moments observed = volume_moments(solids[index]);
      for (int moment = 0; moment < 10; ++moment) {
        const double expected = source_solids[index].raw[moment];
        const double scale = std::max(1.0, std::abs(expected));
        if (std::abs(observed.raw[moment] - expected) > 5e-8 * scale + 1e-9)
          fail(stage + " changed an individual material solid's geometric moments");
        total.raw[moment] += observed.raw[moment] * solid_densities[index];
      }
    }
    return total;
  };
  const Moments brep_mass = integrate_reloaded(brep_solids, "B-rep reload");
  const Moments step_mass = integrate_reloaded(step_solids, "STEP reload");
  std::vector<double> total_mass_moments(10, 0.0);
  for (const Moments& entity : entity_mass_moments)
    for (int moment = 0; moment < 10; ++moment)
      total_mass_moments[static_cast<std::size_t>(moment)] += entity.raw[moment];
  const auto emit_vector = [](const std::vector<double>& values) {
    std::cout << '[';
    for (std::size_t index = 0; index < values.size(); ++index) {
      if (index) std::cout << ',';
      std::cout << values[index];
    }
    std::cout << ']';
  };
  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_FASTER_PLANAR_CAD_METHOD\","
               "\"status\":\"PLANAR_LAYERED_STEP_BREP_ROUNDTRIP_VALIDATED\","
               "\"material_solid_count\":" << source_solids.size()
            << ",\"material_entity_ids\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '\"' << entities[index].identity << '\"';
  }
  std::cout << "],\"entity_roles\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '\"' << entities[index].role << '\"';
  }
  std::cout << "],\"entity_densities\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << entities[index].density;
  }
  std::cout << "],\"entity_solid_counts\":[";
  for (std::size_t index = 0; index < entity_solid_counts.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << entity_solid_counts[index];
  }
  std::cout << "],\"entity_mass_moments\":[";
  for (std::size_t index = 0; index < entity_mass_moments.size(); ++index) {
    if (index) std::cout << ',';
    emit_raw_moments(entity_mass_moments[index]);
  }
  std::cout << "],\"solver_mass_moments\":";
  emit_vector(total_mass_moments);
  std::cout << ",\"raw_mass_moments\":{\"brep\":";
  emit_vector(std::vector<double>(brep_mass.raw, brep_mass.raw + 10));
  std::cout << ",\"step\":";
  emit_vector(std::vector<double>(step_mass.raw, step_mass.raw + 10));
  std::cout << "}}\n";
  return 0;
}

struct Faster3dCadEntity {
  std::string role;
  std::string identity;
  double density = 0.0;
  std::string mesh_path;
};

TopoDS_Shape negative_plane_halfspace(
    const std::array<double, 10>& coefficients, const TopoDS_Shape& domain) {
  const gp_Vec linear(coefficients[1], coefficients[2], coefficients[3]);
  const double norm_squared = linear.SquareMagnitude();
  if (!std::isfinite(norm_squared) || !(norm_squared > 0.0))
    fail("plane CAD quadric has a zero or nonfinite normal");
  const gp_Vec normal = linear / std::sqrt(norm_squared);
  const gp_Pnt plane_point = gp_Pnt(0.0, 0.0, 0.0).Translated(
      linear * (-coefficients[0] / norm_squared));

  Bnd_Box bounds;
  BRepBndLib::AddOptimal(domain, bounds, Standard_False, Standard_True);
  if (bounds.IsVoid() || bounds.IsOpenXmin() || bounds.IsOpenXmax() ||
      bounds.IsOpenYmin() || bounds.IsOpenYmax() || bounds.IsOpenZmin() ||
      bounds.IsOpenZmax())
    fail("plane Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  const gp_Pnt center(0.5 * (xmin + xmax), 0.5 * (ymin + ymax),
                      0.5 * (zmin + zmax));
  const double diagonal = center.Distance(gp_Pnt(xmax, ymax, zmax)) * 2.0;
  const double extent = 4.0 * (diagonal + center.Distance(plane_point) + 1.0);
  if (!std::isfinite(extent) || !(extent > 0.0))
    fail("plane CAD halfspace face extent is not finite");

  const gp_Pln plane(plane_point, gp_Dir(normal));
  BRepBuilderAPI_MakeFace face(plane, -extent, extent, -extent, extent);
  if (!face.IsDone()) fail("could not build finite support face for plane halfspace");
  const gp_Pnt negative_side = plane_point.Translated(
      normal * -(diagonal + center.Distance(plane_point) + 1.0));
  BRepPrimAPI_MakeHalfSpace halfspace(face.Face(), negative_side);
  if (!halfspace.IsDone() || halfspace.Solid().IsNull())
    fail("OCCT could not construct the negative plane halfspace");
  return halfspace.Solid();
}

TopoDS_Shape finite_elliptic_paraboloid_carrier(
    const std::array<double, 10>& coefficients, const TopoDS_Shape& domain) {
  std::array<double, 10> canonical_coefficients = coefficients;
  const double quadratic_trace = coefficients[4] + coefficients[5] + coefficients[6];
  if (!std::isfinite(quadratic_trace) || std::abs(quadratic_trace) <= 1e-14)
    fail("paraboloid quadratic block has ambiguous definite orientation");
  const double orientation = quadratic_trace > 0.0 ? 1.0 : -1.0;
  if (orientation < 0.0)
    for (double& value : canonical_coefficients) value = -value;
  // The transverse eigenvalues stay in the canonical positive orientation,
  // while the axial linear term determines which direction is the negative
  // halfspace of the actual input polynomial. Keep both facts: the carrier
  // opens into q<0, and its curvature flips when the whole polynomial flips.
  const RankTwoQuadricFrame canonical_quadric =
      rank_two_quadric_from_solver_coefficients(canonical_coefficients);
  const RankTwoQuadricFrame quadric{
      canonical_quadric.axis_point, canonical_quadric.axis,
      canonical_quadric.radial_x, canonical_quadric.radial_y,
      canonical_quadric.lambda_x, canonical_quadric.lambda_y,
      orientation * canonical_quadric.center_level,
      orientation * canonical_quadric.axial_linear};
  const double linear_scale = std::max({
      std::abs(coefficients[1]), std::abs(coefficients[2]),
      std::abs(coefficients[3]), std::abs(coefficients[4]),
      std::abs(coefficients[5]), std::abs(coefficients[6]), 1e-300});
  if (std::abs(quadric.axial_linear) <= 1e-10 * linear_scale)
    fail("paraboloid quadric requires a nonzero linear component along its axis");
  const double transverse_scale_x =
      std::sqrt(std::abs(quadric.axial_linear) / quadric.lambda_x);
  const double transverse_scale_y =
      std::sqrt(std::abs(quadric.axial_linear) / quadric.lambda_y);
  if (!std::isfinite(transverse_scale_x) || !std::isfinite(transverse_scale_y) ||
      !(transverse_scale_x > 0.0) || !(transverse_scale_y > 0.0))
    fail("paraboloid transverse scaling is non-finite or degenerate");
  const gp_Dir canonical_axis = quadric.axial_linear < 0.0
      ? quadric.axis
      : gp_Dir(-quadric.axis.X(), -quadric.axis.Y(), -quadric.axis.Z());
  const gp_Pnt vertex = quadric.axis_point.Translated(
      gp_Vec(quadric.axis) * (quadric.center_level / quadric.axial_linear));

  Bnd_Box bounds;
  BRepBndLib::AddOptimal(domain, bounds, Standard_False, Standard_True);
  if (bounds.IsVoid() || bounds.IsOpenXmin() || bounds.IsOpenXmax() ||
      bounds.IsOpenYmin() || bounds.IsOpenYmax() || bounds.IsOpenZmin() ||
      bounds.IsOpenZmax())
    fail("paraboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  double maximum_radius_squared = 0.0;
  double maximum_axial = -std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const gp_Vec relative(vertex, corner);
        const double local_x = relative.Dot(gp_Vec(quadric.radial_x)) /
                               transverse_scale_x;
        const double local_y = relative.Dot(gp_Vec(quadric.radial_y)) /
                               transverse_scale_y;
        maximum_radius_squared = std::max(
            maximum_radius_squared, local_x * local_x + local_y * local_y);
        maximum_axial = std::max(maximum_axial,
                                 relative.Dot(gp_Vec(canonical_axis)));
      }
  const double diagonal = std::sqrt((xmax - xmin) * (xmax - xmin) +
                                    (ymax - ymin) * (ymax - ymin) +
                                    (zmax - zmin) * (zmax - zmin));
  const double margin = std::max(1e-7, 0.1 * diagonal);
  const double radial_limit = std::sqrt(maximum_radius_squared) +
                              std::max(1e-7, 0.1 * std::sqrt(maximum_radius_squared));
  const double cap_height = std::max(maximum_axial, radial_limit * radial_limit) + margin;
  if (!std::isfinite(radial_limit) || !std::isfinite(cap_height) ||
      !(radial_limit > 0.0) || !(cap_height > maximum_axial))
    fail("could not construct a finite paraboloid carrier outside the design domain");

  const gp_Pnt p0(0., 0., 0.), p1(0., 0., cap_height);
  const gp_Pnt p2(radial_limit, 0., cap_height);
  // The ruled cap closes the bounded carrier on the negative-input side. If
  // the input polynomial was globally flipped to get the canonical positive
  // quadratic frame, its zero-set curvature in this carrier frame is flipped
  // as well (the selected Boolean side remains q_input < 0).
  const gp_Pnt p3(radial_limit, 0., orientation * radial_limit * radial_limit);
  BRepBuilderAPI_MakeWire wire;
  wire.Add(BRepBuilderAPI_MakeEdge(p0, p1).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
  wire.Add(BRepBuilderAPI_MakeEdge(p2, p3).Edge());
  TColgp_Array1OfPnt poles(1, 3);
  poles.SetValue(1, p0);
      poles.SetValue(2, gp_Pnt(0.5 * radial_limit, 0., 0.));
  poles.SetValue(3, p3);
  Handle(Geom_BezierCurve) parabola = new Geom_BezierCurve(poles);
  TopoDS_Edge parabola_edge = BRepBuilderAPI_MakeEdge(parabola).Edge();
  parabola_edge.Reverse();
  wire.Add(parabola_edge);
  if (!wire.IsDone()) fail("failed to close the paraboloid meridian profile");
  const gp_Pln profile_plane(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.));
  BRepBuilderAPI_MakeFace profile_face(profile_plane, wire.Wire(), Standard_True);
  if (!profile_face.IsDone()) fail("failed to construct the paraboloid meridian face");
  const gp_Ax1 revolution_axis(gp_Pnt(0., 0., 0.), gp_Dir(0., 0., 1.));
  BRepPrimAPI_MakeRevol revolution(profile_face.Face(), revolution_axis,
                                    2.0 * std::acos(-1.0), Standard_True);
  if (!revolution.IsDone() || revolution.Shape().IsNull())
    fail("failed to construct the finite paraboloid carrier");
  const double map[3][3] = {
      {transverse_scale_x * quadric.radial_x.X(),
       transverse_scale_y * quadric.radial_y.X(), canonical_axis.X()},
      {transverse_scale_x * quadric.radial_x.Y(),
       transverse_scale_y * quadric.radial_y.Y(), canonical_axis.Y()},
      {transverse_scale_x * quadric.radial_x.Z(),
       transverse_scale_y * quadric.radial_y.Z(), canonical_axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_GTrsf transform(vectorial, gp_XYZ(vertex.X(), vertex.Y(), vertex.Z()));
  BRepBuilderAPI_GTransform affine(revolution.Shape(), transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not place the finite paraboloid carrier in solver coordinates");
  const TopoDS_Shape carrier = affine.Shape();
  if (!BRepCheck_Analyzer(carrier, Standard_True).IsValid())
    fail("elliptic-paraboloid carrier failed OCCT topology validation");
  return carrier;
}

TopoDS_Shape finite_hyperbolic_paraboloid_carrier(
    const std::array<double, 10>& negative_inside_coefficients,
    const TopoDS_Shape& domain) {
  const RankTwoQuadricFrame quadric =
      rank_two_quadric_from_solver_coefficients(
          negative_inside_coefficients, true);
  const double matrix_scale = std::max(std::abs(quadric.lambda_x),
                                       std::abs(quadric.lambda_y));
  if (std::abs(quadric.axial_linear) <= 1e-10 * matrix_scale)
    fail("hyperbolic paraboloid requires a nonzero linear component along its null axis");
  const double scale_x = std::sqrt(std::abs(quadric.axial_linear) /
                                   std::abs(quadric.lambda_x));
  const double scale_y = std::sqrt(std::abs(quadric.axial_linear) /
                                   std::abs(quadric.lambda_y));
  if (!std::isfinite(scale_x) || !std::isfinite(scale_y) ||
      !(scale_x > 0.0) || !(scale_y > 0.0))
    fail("hyperbolic-paraboloid transverse scales are invalid");

  // Canonical local coordinates satisfy q = u^2 - v^2 + w, so the
  // unbounded q<=0 halfspace is below the saddle w=-u^2+v^2.
  const gp_Dir canonical_axis = quadric.axial_linear < 0.0
      ? gp_Dir(-quadric.axis.X(), -quadric.axis.Y(), -quadric.axis.Z())
      : quadric.axis;
  const gp_Pnt vertex = quadric.axis_point.Translated(
      gp_Vec(quadric.axis) *
      (quadric.center_level / quadric.axial_linear));

  Bnd_Box bounds;
  BRepBndLib::AddOptimal(domain, bounds, Standard_False, Standard_True);
  if (bounds.IsVoid() || bounds.IsOpenXmin() || bounds.IsOpenXmax() ||
      bounds.IsOpenYmin() || bounds.IsOpenYmax() ||
      bounds.IsOpenZmin() || bounds.IsOpenZmax())
    fail("hyperbolic-paraboloid Boolean requires a finite design-domain bounding box");
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bounds.Get(xmin, ymin, zmin, xmax, ymax, zmax);
  double u_extent = 0.0;
  double v_extent = 0.0;
  double minimum_w = std::numeric_limits<double>::infinity();
  for (int x = 0; x < 2; ++x)
    for (int y = 0; y < 2; ++y)
      for (int z = 0; z < 2; ++z) {
        const gp_Pnt corner(x == 0 ? xmin : xmax, y == 0 ? ymin : ymax,
                            z == 0 ? zmin : zmax);
        const gp_Vec relative(vertex, corner);
        u_extent = std::max(u_extent,
            std::abs(relative.Dot(gp_Vec(quadric.radial_x)) / scale_x));
        v_extent = std::max(v_extent,
            std::abs(relative.Dot(gp_Vec(quadric.radial_y)) / scale_y));
        minimum_w = std::min(minimum_w,
            relative.Dot(gp_Vec(canonical_axis)));
      }
  const double diagonal = std::sqrt((xmax - xmin) * (xmax - xmin) +
                                    (ymax - ymin) * (ymax - ymin) +
                                    (zmax - zmin) * (zmax - zmin));
  const double margin = std::max(1e-7, 0.1 * diagonal);
  u_extent += std::max(1e-7, 0.1 * u_extent);
  v_extent += std::max(1e-7, 0.1 * v_extent);
  const double carrier_depth = v_extent * v_extent - minimum_w + margin;
  if (!std::isfinite(u_extent) || !std::isfinite(v_extent) ||
      !std::isfinite(carrier_depth) || !(u_extent > 0.0) ||
      !(v_extent > 0.0) || !(carrier_depth > 0.0))
    fail("could not construct a finite hyperbolic-paraboloid carrier outside the domain");

  TColgp_Array2OfPnt poles(1, 3, 1, 3);
  const double u_poles[3] = {-u_extent, 0.0, u_extent};
  const double u_square_poles[3] = {-u_extent * u_extent,
                                     u_extent * u_extent,
                                    -u_extent * u_extent};
  const double v_poles[3] = {-v_extent, 0.0, v_extent};
  const double v_square_poles[3] = {v_extent * v_extent,
                                    -v_extent * v_extent,
                                     v_extent * v_extent};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      poles.SetValue(i + 1, j + 1,
          gp_Pnt(u_poles[i], v_poles[j],
                 u_square_poles[i] + v_square_poles[j]));
  Handle(Geom_BezierSurface) saddle_surface = new Geom_BezierSurface(poles);
  BRepBuilderAPI_MakeFace face(saddle_surface, 1e-9);
  if (!face.IsDone())
    fail("failed to construct exact hyperbolic-paraboloid surface face");
  BRepPrimAPI_MakePrism prism(face.Face(),
      gp_Vec(0.0, 0.0, -carrier_depth), Standard_True, Standard_True);
  if (!prism.IsDone() || prism.Shape().IsNull())
    fail("failed to construct finite hyperbolic-paraboloid carrier");

  const double map[3][3] = {
      {scale_x * quadric.radial_x.X(), scale_y * quadric.radial_y.X(),
       canonical_axis.X()},
      {scale_x * quadric.radial_x.Y(), scale_y * quadric.radial_y.Y(),
       canonical_axis.Y()},
      {scale_x * quadric.radial_x.Z(), scale_y * quadric.radial_y.Z(),
       canonical_axis.Z()}};
  const gp_Mat vectorial(map[0][0], map[0][1], map[0][2],
                         map[1][0], map[1][1], map[1][2],
                         map[2][0], map[2][1], map[2][2]);
  const gp_GTrsf transform(vectorial,
      gp_XYZ(vertex.X(), vertex.Y(), vertex.Z()));
  BRepBuilderAPI_GTransform affine(prism.Shape(), transform, Standard_True);
  if (!affine.IsDone() || affine.Shape().IsNull())
    fail("OCCT could not place the hyperbolic-paraboloid carrier");
  const TopoDS_Shape carrier = affine.Shape();
  if (!BRepCheck_Analyzer(carrier, Standard_True).IsValid())
    fail("hyperbolic-paraboloid carrier failed OCCT topology validation");
  return carrier;
}

int run_faster_3d_quadric_material_roundtrip(int argc, char** argv) {
  if (argc != 5)
    fail("usage: occt_step_roundtrip --faster-3d-quadric-material INPUT.txt MATERIAL.brep MATERIAL.step");
  std::ifstream input(argv[2]);
  if (!input) fail("could not open Faster 3D material protocol");
  std::string token;
  if (!(input >> token) ||
      (token != "FASTER_3D_QUADRIC_MATERIAL_V1" &&
       token != "FASTER_3D_QUADRIC_MATERIAL_V2"))
    fail("unsupported Faster 3D material protocol header");
  const bool typed_protocol = token == "FASTER_3D_QUADRIC_MATERIAL_V2";
  int positive_region = -1;
  if (!(input >> positive_region) || (positive_region != 0 && positive_region != 1))
    fail("Faster quadric positive-region flag must be 0 or 1");
  std::array<double, 10> coefficients{};
  for (double& value : coefficients) {
    if (!(input >> token)) fail("truncated Faster quadric coefficient record");
    value = parse_finite_double(token, "Faster quadric coefficient");
  }
  std::string quadric_type = "ELLIPSOID";
  if (typed_protocol && !(input >> quadric_type))
    fail("Faster quadric type is missing from V2 material protocol");
  if (quadric_type != "PLANE" && quadric_type != "ELLIPSOID" &&
      quadric_type != "PARABOLOID" &&
      quadric_type != "HYPERBOLIC_PARABOLOID" &&
      quadric_type != "ELLIPTIC_CYLINDER" &&
      quadric_type != "ONE_SHEET_HYPERBOLOID" &&
      quadric_type != "TWO_SHEET_HYPERBOLOID")
    fail("unsupported Faster CAD quadric type: " + quadric_type);
  std::array<double, 10> negative_inside_coefficients = coefficients;
  if (positive_region == 1)
    for (double& value : negative_inside_coefficients) value = -value;
  std::array<double, 10> canonical_saddle_coefficients = coefficients;
  int saddle_coefficient_orientation = 1;
  if (quadric_type == "HYPERBOLIC_PARABOLOID") {
    bool found_nonzero_coefficient = false;
    for (const double value : coefficients) {
      if (value == 0.0) continue;
      saddle_coefficient_orientation = value < 0.0 ? -1 : 1;
      found_nonzero_coefficient = true;
      break;
    }
    if (!found_nonzero_coefficient)
      fail("hyperbolic-paraboloid coefficients are all zero");
    if (saddle_coefficient_orientation < 0)
      for (double& value : canonical_saddle_coefficients) value = -value;
  }
  if (quadric_type == "PLANE") {
    const double normal_squared = coefficients[1] * coefficients[1] +
                                  coefficients[2] * coefficients[2] +
                                  coefficients[3] * coefficients[3];
    if (!(normal_squared > 0.0) || !std::isfinite(normal_squared))
      fail("plane CAD quadric has a zero or nonfinite normal");
  } else if (quadric_type == "ELLIPSOID")
    (void)ellipsoid_from_solver_coefficients(negative_inside_coefficients);
  else if (quadric_type == "PARABOLOID") {
    // Validated against the bounded design domain during exact carrier construction.
  }
  else if (quadric_type == "HYPERBOLIC_PARABOLOID")
    (void)rank_two_quadric_from_solver_coefficients(
        negative_inside_coefficients, true);
  else if (quadric_type == "ELLIPTIC_CYLINDER")
    (void)rank_two_quadric_from_solver_coefficients(
        negative_inside_coefficients);
  else if (quadric_type == "ONE_SHEET_HYPERBOLOID")
    (void)hyperboloid_from_solver_coefficients(negative_inside_coefficients);
  else
    (void)hyperboloid_from_solver_coefficients(negative_inside_coefficients, true);
  std::size_t entity_count = 0;
  if (!(input >> entity_count) || entity_count == 0 || entity_count > 10000)
    fail("Faster material entity count must be in [1,10000]");
  std::vector<Faster3dCadEntity> entities;
  entities.reserve(entity_count);
  for (std::size_t index = 0; index < entity_count; ++index) {
    Faster3dCadEntity entity;
    if (!(input >> token) || token != "ENTITY" ||
        !(input >> entity.role >> entity.identity >> token >> std::quoted(entity.mesh_path)))
      fail("malformed Faster material entity record");
    entity.density = parse_finite_double(token, "Faster material density");
    if (!(entity.density > 0.0) || entity.identity.empty() || entity.mesh_path.empty())
      fail("Faster material entity identity, path, and density must be valid");
    if (index == 0 && entity.role != "design")
      fail("first Faster material entity must be the editable design region");
    if (index > 0 && entity.role != "fixed")
      fail("subsequent Faster material entities must be fixed regions");
    entities.push_back(std::move(entity));
  }
  if (input >> token) fail("unexpected trailing Faster 3D material protocol data");

  Message::DefaultMessenger()->ChangePrinters().Clear();
  BRep_Builder compound_builder;
  TopoDS_Compound material;
  compound_builder.MakeCompound(material);
  std::vector<TopoDS_Shape> source_solids;
  std::vector<double> solid_densities;
  std::vector<Moments> entity_mass_moments;
  std::vector<std::size_t> entity_solid_counts;
  for (std::size_t index = 0; index < entities.size(); ++index) {
    const Faster3dCadEntity& entity = entities[index];
    TopoDS_Shape source = mesh_to_solid(entity.mesh_path, 1.0);
    if (index == 0) {
      TopoDS_Shape quadric;
      if (quadric_type == "PLANE") {
        quadric = negative_plane_halfspace(negative_inside_coefficients, source);
      } else if (quadric_type == "ELLIPSOID") {
        const EllipsoidQuadric ellipsoid =
            ellipsoid_from_solver_coefficients(negative_inside_coefficients);
        const TopoDS_Shape unit_sphere =
            BRepPrimAPI_MakeSphere(gp_Pnt(0.0, 0.0, 0.0), 1.0).Shape();
        BRepBuilderAPI_GTransform affine(
            unit_sphere, ellipsoid.unit_sphere_map, Standard_True);
        if (!affine.IsDone() || affine.Shape().IsNull())
          fail("OCCT could not construct the exact Faster ellipsoid");
        quadric = affine.Shape();
      } else if (quadric_type == "ELLIPTIC_CYLINDER") {
        quadric = finite_elliptic_cylinder_carrier(
            negative_inside_coefficients, source);
      } else if (quadric_type == "PARABOLOID") {
        quadric = finite_elliptic_paraboloid_carrier(
            negative_inside_coefficients, source);
      } else if (quadric_type == "HYPERBOLIC_PARABOLOID") {
        quadric = finite_hyperbolic_paraboloid_carrier(
            canonical_saddle_coefficients, source);
      } else if (quadric_type == "ONE_SHEET_HYPERBOLOID") {
        quadric = finite_one_sheet_hyperboloid_carrier(
            negative_inside_coefficients, source);
      } else {
        quadric = finite_two_sheet_hyperboloid_carrier(
            negative_inside_coefficients, source);
      }
      // OCCT's half-space primitive is deliberately unbounded, so its carrier
      // is not a closed finite solid and is not subject to the closed-carrier
      // analyzer used for bounded quadrics. The finite design-domain Boolean
      // result is still fully validated below.
      if (quadric_type != "PLANE" &&
          !BRepCheck_Analyzer(quadric, Standard_True).IsValid())
        fail("exact Faster quadric failed OCCT topology validation");
      // Build one sign-canonical saddle B-rep for p and -p. OCCT's Boolean
      // builder can leave a narrow gap when complementary operations use
      // geometrically identical saddles with independently parameterized
      // surfaces. For p = s*C, where C is canonical, the retained p>=0 side
      // is C<=0 (Common) when s<0 and C>=0 (Cut) when s>0.
      const bool saddle_common =
          quadric_type == "HYPERBOLIC_PARABOLOID" &&
          saddle_coefficient_orientation < 0;
      if (quadric_type == "HYPERBOLIC_PARABOLOID" ? saddle_common
                                                     : positive_region == 1) {
        BRepAlgoAPI_Common common(source, quadric);
        common.SetRunParallel(Standard_False);
        common.Build();
        if (!common.IsDone() || common.Shape().IsNull())
          fail("exact Faster design/quadric intersection failed");
        source = common.Shape();
      } else {
        BRepAlgoAPI_Cut cut(source, quadric);
        cut.SetRunParallel(Standard_False);
        cut.Build();
        if (!cut.IsDone() || cut.Shape().IsNull())
          fail("exact Faster design-minus-quadric Boolean failed");
        source = cut.Shape();
      }
    }
    require_valid_topology(source,"Faster material entity "+entity.identity,
                          "Faster material entity failed OCCT topology validation");
    BRep_Builder entity_builder;
    TopoDS_Compound entity_shape;
    entity_builder.MakeCompound(entity_shape);
    std::size_t solid_count = 0;
    for (TopExp_Explorer it(source, TopAbs_SOLID); it.More(); it.Next()) {
      const TopoDS_Shape solid = it.Current();
      const Moments moments = volume_moments(solid);
      if (!(moments.volume > 0.0)) fail("Faster CAD entity has non-positive volume");
      compound_builder.Add(material, solid);
      entity_builder.Add(entity_shape, solid);
      source_solids.push_back(solid);
      solid_densities.push_back(entity.density);
      ++solid_count;
    }
    if (solid_count == 0) fail("Faster material Boolean produced an empty entity");
    entity_solid_counts.push_back(solid_count);
    entity_mass_moments.push_back(weighted_moments(volume_moments(entity_shape), entity.density));
  }
  validate_roundtrip_solid(material, "Faster 3D material compound");
  const std::string brep_path(argv[3]), step_path(argv[4]);
  if (brep_path == step_path) fail("B-rep and STEP output paths must differ");
  if (!BRepTools::Write(material, brep_path.c_str()))
    fail("failed to write Faster 3D material B-rep");
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(material, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);
  validate_roundtrip_solid(brep_shape, "Faster 3D B-rep reload");
  validate_roundtrip_solid(step_shape, "Faster 3D STEP reload");
  const auto solids_in = [](const TopoDS_Shape& shape) {
    std::vector<TopoDS_Shape> result;
    for (TopExp_Explorer it(shape, TopAbs_SOLID); it.More(); it.Next())
      result.push_back(it.Current());
    return result;
  };
  const std::vector<TopoDS_Shape> brep_solids = solids_in(brep_shape);
  const std::vector<TopoDS_Shape> step_solids = solids_in(step_shape);
  if (brep_solids.size() != source_solids.size() || step_solids.size() != source_solids.size())
    fail("Faster 3D B-rep/STEP round-trip changed material solid count");
  const auto integrate_reloaded = [&](const std::vector<TopoDS_Shape>& solids,
                                      const std::string& stage) {
    Moments total{};
    for (std::size_t index = 0; index < solids.size(); ++index) {
      const Moments observed = volume_moments(solids[index]);
      const Moments expected = volume_moments(source_solids[index]);
      for (int moment = 0; moment < 10; ++moment) {
        const double scale = std::max(1.0, std::abs(expected.raw[moment]));
        if (std::abs(observed.raw[moment] - expected.raw[moment]) > 5e-8 * scale + 1e-10)
          fail(stage + " changed a material solid's geometric moments");
        total.raw[moment] += observed.raw[moment] * solid_densities[index];
      }
    }
    return total;
  };
  const Moments brep_mass = integrate_reloaded(brep_solids, "Faster B-rep reload");
  const Moments step_mass = integrate_reloaded(step_solids, "Faster STEP reload");
  Moments solver_mass{};
  for (const Moments& entity : entity_mass_moments)
    for (int moment = 0; moment < 10; ++moment)
      solver_mass.raw[moment] += entity.raw[moment];
  const auto emit_vector = [](const double* values, std::size_t count) {
    std::cout << '[';
    for (std::size_t index = 0; index < count; ++index) {
      if (index) std::cout << ',';
      std::cout << values[index];
    }
    std::cout << ']';
  };
  std::cout << std::setprecision(17)
            << "{\"classification\":\"METHOD_VALIDATION_ONLY\","
               "\"status\":\"FASTER_3D_MATERIAL_STEP_BREP_VALIDATED\","
               "\"quadric_type\":\"" << quadric_type << "\","
               "\"boolean_backend\":\"OCCT_EXACT_ANALYTIC\","
               "\"occt_version\":\"" << OCC_VERSION_COMPLETE << "\","
               "\"material_solid_count\":" << source_solids.size() << ",\"material_entity_ids\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '\"' << entities[index].identity << '\"';
  }
  std::cout << "],\"entity_roles\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << '\"' << entities[index].role << '\"';
  }
  std::cout << "],\"entity_densities\":[";
  for (std::size_t index = 0; index < entities.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << entities[index].density;
  }
  std::cout << "],\"entity_solid_counts\":[";
  for (std::size_t index = 0; index < entity_solid_counts.size(); ++index) {
    if (index) std::cout << ',';
    std::cout << entity_solid_counts[index];
  }
  std::cout << "],\"entity_mass_moments\":[";
  for (std::size_t index = 0; index < entity_mass_moments.size(); ++index) {
    if (index) std::cout << ',';
    emit_raw_moments(entity_mass_moments[index]);
  }
  std::cout << "],\"cad_source_mass_moments\":";
  emit_vector(solver_mass.raw, 10);
  std::cout << ",\"raw_mass_moments\":{\"brep\":";
  emit_vector(brep_mass.raw, 10);
  std::cout << ",\"step\":";
  emit_vector(step_mass.raw, 10);
  std::cout << "}}\n";
  return 0;
}

int run_restore_material_frame(int argc, char** argv) {
  if (argc != 10)
    fail("usage: --restore-material-frame INPUT.brep L OX OY OZ UNIT OUTPUT.brep OUTPUT.step");
  Message::DefaultMessenger()->ChangePrinters().Clear();
  const double scale = parse_coordinate_scale(argv[3]);
  const gp_XYZ offset(parse_finite_double(argv[4], "origin x"),
                      parse_finite_double(argv[5], "origin y"),
                      parse_finite_double(argv[6], "origin z"));
  const std::string unit(argv[7]), brep_path(argv[8]), step_path(argv[9]);
  double mm_per_unit;
  DESTEP_Parameters params;
  if (unit == "m") { mm_per_unit=1000.; params.WriteUnit=UnitsMethods_LengthUnit_Meter; }
  else if (unit == "cm") { mm_per_unit=10.; params.WriteUnit=UnitsMethods_LengthUnit_Centimeter; }
  else if (unit == "mm") { mm_per_unit=1.; params.WriteUnit=UnitsMethods_LengthUnit_Millimeter; }
  else fail("unsupported original length unit (expected m, cm or mm)");
  if (std::filesystem::absolute(brep_path).lexically_normal() ==
      std::filesystem::absolute(step_path).lexically_normal()) fail("CAD output paths must differ");
  for (const auto& path : {brep_path,step_path})
    if (std::filesystem::exists(path) || std::filesystem::is_symlink(path))
      fail("refusing to overwrite CAD output");
  const TopoDS_Shape input=read_brep(argv[2]);
  validate_roundtrip_solid(input,"normalized CAD source");
  gp_GTrsf frame;
  frame.SetVectorialPart(gp_Mat(scale,0.,0.,0.,scale,0.,0.,0.,scale));
  frame.SetTranslationPart(offset);
  BRepBuilderAPI_GTransform transformed(input,frame,Standard_True);
  if (!transformed.IsDone()) fail("original coordinate transformation failed");
  const TopoDS_Shape material=transformed.Shape();
  validate_roundtrip_solid(material,"original-coordinate CAD");
  if (!BRepTools::Write(material,brep_path.c_str())) fail("failed to write original B-rep");
  STEPControl_Writer writer;
  writer.Model()->SetLocalLengthUnit(mm_per_unit);
  if (writer.Transfer(material,STEPControl_AsIs,params) != IFSelect_RetDone)
    fail("failed to transfer original-unit STEP");
  if (writer.Write(step_path.c_str()) != IFSelect_RetDone) fail("failed to write original-unit STEP");
  const auto step_in_unit = [&](double length_unit) {
    STEPControl_Reader reader;
    if (reader.ReadFile(step_path.c_str()) != IFSelect_RetDone) fail("cannot reload original STEP");
    reader.SetSystemLengthUnit(length_unit);
    if (reader.TransferRoots() <= 0) fail("cannot transfer original STEP");
    return reader.OneShape();
  };
  const std::vector<std::pair<std::string,TopoDS_Shape>> stages={
    {"source",input},{"transformed",material},{"brep",read_brep(brep_path)},
    {"step",step_in_unit(mm_per_unit)},{"step_mm",step_in_unit(1.)}};
  std::size_t expected_count=0;
  std::cout << std::setprecision(17) << "{\"schema\":\"original-cad-frame-v1\",\"length_unit\":\""
            << unit << "\",\"millimetres_per_unit\":" << mm_per_unit
            << ",\"solid_geometric_moments\":{";
  for (std::size_t s=0;s<stages.size();++s) {
    validate_roundtrip_solid(stages[s].second,stages[s].first);
    if (s) std::cout << ',';
    std::cout << '\"' << stages[s].first << "\":[";
    std::size_t count=0;
    for (TopExp_Explorer it(stages[s].second,TopAbs_SOLID);it.More();it.Next()) {
      validate_roundtrip_solid(it.Current(),stages[s].first+" solid");
      if (count++) std::cout << ',';
      emit_raw_moments(volume_moments(it.Current()));
    }
    if (s==0) expected_count=count;
    if (count==0 || count!=expected_count) fail("coordinate round trip changed solid count");
    std::cout << ']';
  }
  std::cout << "}}\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc >= 2 && std::string(argv[1]) == "--restore-material-frame")
    return run_restore_material_frame(argc,argv);
  if (argc >= 2 && std::string(argv[1]) == "--mesh") {
    return run_mesh_roundtrip(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-sphere") {
    return run_solver_sphere_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-ellipsoid") {
    return run_solver_ellipsoid_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-cylinder") {
    return run_solver_cylinder_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-paraboloid") {
    return run_solver_paraboloid_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-hyperbolic-paraboloid") {
    return run_solver_hyperbolic_paraboloid_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-hyperboloid-one-sheet") {
    return run_solver_one_sheet_hyperboloid_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--quadric-hyperboloid-two-sheet") {
    return run_solver_two_sheet_hyperboloid_quadric(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--planar-layered-extrusion") {
    return run_planar_layered_extrusion_roundtrip(argc, argv);
  }
  if (argc >= 2 && std::string(argv[1]) == "--faster-3d-quadric-material") {
    return run_faster_3d_quadric_material_roundtrip(argc, argv);
  }
  if (argc != 3) {
    fail("usage: occt_step_roundtrip MATERIAL.brep MATERIAL.step");
  }
  const std::string brep_path(argv[1]);
  const std::string step_path(argv[2]);

  // OCCT's default STEP transfer messenger writes informational progress to
  // stdout.  This single-purpose probe reserves stdout for its JSON evidence;
  // return statuses below remain mandatory and failures are sent to stderr.
  Message::DefaultMessenger()->ChangePrinters().Clear();

  const TopoDS_Shape outer = BRepPrimAPI_MakeBox(2.0, 2.0, 2.0).Shape();
  const TopoDS_Shape cavity = BRepPrimAPI_MakeSphere(gp_Pnt(1.0, 1.0, 1.0), 0.5).Shape();
  BRepAlgoAPI_Cut cut(outer, cavity);
  if (!cut.IsDone()) {
    fail("box-minus-sphere Boolean did not complete");
  }
  const TopoDS_Shape original = cut.Shape();

  if (!BRepTools::Write(original, brep_path.c_str())) {
    fail("failed to write B-rep material artifact");
  }
  const TopoDS_Shape brep_shape = read_brep(brep_path);
  write_step(original, step_path);
  const TopoDS_Shape step_shape = read_step(step_path);

  const double analytic_volume = 8.0 - std::acos(-1.0) / 6.0;
  const Moments original_moments = volume_moments(original);
  const Moments brep_moments = volume_moments(brep_shape);
  const Moments step_moments = volume_moments(step_shape);
  std::cout << std::setprecision(17)
            << "{\"classification\":\"SYNTHETIC_REGRESSION\","
               "\"status\":\"SYNTHETIC_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED\","
               "\"analytic_volume\":"
            << analytic_volume << ",\"moments\":{\"original\":";
  emit_moments(original_moments);
  std::cout << ",\"brep\":";
  emit_moments(brep_moments);
  std::cout << ",\"step\":";
  emit_moments(step_moments);
  std::cout << "}}\n";
  return 0;
}
