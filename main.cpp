#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <math.h>

#include <Eigen/Geometry>

#include <imgui/imgui.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/bounding_box.h>
#include <igl/grid.h>

#include "include/grid_util.h"
#include "include/inner_void_mesh.h"
#include "include/physics/MassProperties.h"
#include "include/stability/SupportRegion.h"
#include "include/viewer/ResultView.h"


#define PI 3.14159265

const bool DEBUG = false;

const std::string DEFAULT_MESH_FILE = "../data/bunny.off";

// Viewer data indices
int mesh_data_id;
int inner_data_mesh_id;
int design_data_mesh_id;
int material_data_mesh_id;
int plane_data_id;
int axes_data_id;
int gravity_data_id;
int carve_plane_data_id;

const Eigen::Vector3d DEFAULT_GRAVITY(0.0,-1.0,0.0);
const float DEFAULT_GRAVITY_YAW = 0.0;
const float DEFAULT_GRAVITY_PITCH = 0.0;
const float DEFAULT_GRAVITY_ROLL = 180.0;
const float DEFAULT_CONTACT_THRESHOLD = 0.0f;
const float DEFAULT_TOPPLING_ANGLE = 5.0f;
const float DEFAULT_MINIMUM_WALL_THICKNESS = 0.0f;

// Predefined colors
const Eigen::RowVector3d orange(1.0,0.7,0.2);
const Eigen::RowVector3d yellow(1.0,0.9,0.2);
const Eigen::RowVector3d blue(0.2,0.3,0.8);
const Eigen::RowVector3d green(0.2,0.6,0.3);
const Eigen::RowVector3d black(0.0,0.0,0.0);
const Eigen::RowVector4d resultShellColor(0.78, 0.84, 0.95, 0.24);
const Eigen::RowVector4d resultVoidColor(0.90, 0.10, 0.08, 0.94);
const Eigen::RowVector4d resultGroundColor(0.50, 0.66, 0.36, 0.70);
const Eigen::RowVector4d spinDesignColor(0.18, 0.48, 0.88, 0.16);
const Eigen::RowVector4d spinMaterialColor(0.95, 0.68, 0.12, 0.82);

// GUI objects
igl::opengl::glfw::Viewer viewer;
igl::opengl::glfw::imgui::ImGuiMenu menu;
bool resultViewMode = false;
bool spinResultViewMode = false;
std::string resultWindowTitle = "Make It Stand — sealed result";
bool showSpinOuter = true;
bool showSpinDesign = true;
bool showSpinMaterial = true;
bool showSpinCavity = true;
bool showSpinAxes = true;
bool showSpinTargetPlane = true;

// Undeformed objects
Eigen::MatrixXd origV;
Eigen::MatrixXi origF;

struct State
{
	// object current position/orientation
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	Eigen::MatrixXd innerV;
	Eigen::MatrixXi innerF;
	Eigen::MatrixXd designV;
	Eigen::MatrixXi designF;
	Eigen::MatrixXd materialV;
	Eigen::MatrixXi materialF;
	result_view::SpinResultViewMetadata spinMetadata;
	double spinAxisDisplayLength = 1.0;

	Eigen::MatrixXd planeV;
	Eigen::MatrixXi planeF;
	Eigen::RowVector3d planeCenter;

	int balancePointIdx;

	// gravity unit vector (relative to original mesh position)
	Eigen::Vector3d gravity;

	// object orientation 
	float yaw;
	float pitch;	
	float roll;
	float contactThreshold;
	float topplingAngle;
	float minimumWallThickness;
	std::string stabilityMessage;

	// user workflow state
	bool selectBalancePoint = true;
	bool selectOrientation = false;
	bool isCarving = false;

	// carving state
	InnerVoidMesh *innerMesh = NULL;
} state;

double sin_deg(double angle)
{
	return std::sin(angle * PI / 180.0);
}

double cos_deg(double angle)
{
	return std::cos(angle * PI / 180.0);
}

Eigen::RowVector3d getBalancePoint()
{
	return state.V.row(state.balancePointIdx);
}

// Update display of balance point
void updateBalancePoint()
{
	viewer.data(mesh_data_id).clear_points();
	viewer.data(mesh_data_id).set_points(getBalancePoint(), state.selectBalancePoint ? yellow : blue);
}


void updateGravityVector() {
	// Convert yaw, pitch, roll to vector
	// - source: https://stackoverflow.com/questions/1568568/how-to-convert-euler-angles-to-directional-vector
	state.gravity(0) = -cos_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) - sin_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(1) = -sin_deg(state.yaw) * sin_deg(state.pitch) * sin_deg(state.roll) + cos_deg(state.yaw) * cos_deg(state.roll);
	state.gravity(2) = cos_deg(state.pitch) * sin_deg(state.roll);
}

void updateGravity()
{
	updateGravityVector();

	if (DEBUG)
	{
		viewer.data(gravity_data_id).clear_points();
		viewer.data(gravity_data_id).add_points(Eigen::RowVector3d(0.0, 0.0, 0.0), black);
		viewer.data(gravity_data_id).add_points(state.gravity.transpose() / 10.0, black);
	}
}

void updateMesh()
{
	// Rotate mesh in the opposite direction from gravity
	// - source: https://stackoverflow.com/questions/1568568/how-to-convert-euler-angles-to-directional-vector
	Eigen::Matrix3d R, rollR, pitchR, yawR;
	rollR <<
		1, 0, 0,
		0, cos_deg((double)state.roll-DEFAULT_GRAVITY_ROLL), -sin_deg((double)state.roll-DEFAULT_GRAVITY_ROLL),
		0, sin_deg((double)state.roll-DEFAULT_GRAVITY_ROLL), cos_deg((double)state.roll-DEFAULT_GRAVITY_ROLL);
	pitchR <<
		cos_deg(state.pitch), 0, -sin_deg(state.pitch),
		0, 1, 0,
		sin_deg(state.pitch), 0, cos_deg(state.pitch);
	yawR <<
		cos_deg(state.yaw), -sin_deg(state.yaw), 0,
		sin_deg(state.yaw), cos_deg(state.yaw), 0,
		0, 0, 1;
	
	state.V = origV * rollR.transpose();
	state.V *= pitchR.transpose();
	state.V *= yawR.transpose();

	// Place the actual lowest geometry on the ground. The selected point is a
	// lateral posture anchor; the support target is derived from all contacts.
	Eigen::RowVector3d translate = state.planeCenter - getBalancePoint();
	translate(1) = state.planeCenter(1) - state.V.col(1).minCoeff();
	for (int i = 0; i < state.V.rows(); i++) 
	{
		state.V.row(i) += translate;
	}

	// Update viewer
	viewer.data(mesh_data_id).set_vertices(state.V);
	updateBalancePoint();
}

void updateInnerMesh()
{
	viewer.data(inner_data_mesh_id).clear();
	viewer.data(inner_data_mesh_id).set_mesh(state.innerV, state.innerF);
}

void clearInnerMesh()
{
	state.innerV.resize(0,3);
	state.innerF.resize(0,3);
	viewer.data(inner_data_mesh_id).clear();

	if (DEBUG)
	{
		viewer.data(carve_plane_data_id).clear_points();
	}
}

bool findLowestPointIdx(const Eigen::MatrixXd &V, int &lowestIdx)
{
	int heightCol = 1;
	double minZ = V.col(heightCol).minCoeff();
	lowestIdx = -1;
	for (int i = 0; i < V.rows(); i++) {
		if (V(i, heightCol) == minZ) {
			lowestIdx = i;
			break;
		}
	}

	if (lowestIdx == -1) {
		return false;
	}

	return true;
}

// Note: nothing related to the viewer yet
void initState(int argc, char *argv[])
{
	const result_view::ResultViewArguments resultArguments =
		result_view::parseResultViewArguments(argc, argv);
	if (resultArguments.enabled)
	{
		Eigen::MatrixXd outerV;
		Eigen::MatrixXi outerF;
		Eigen::MatrixXd innerV;
		Eigen::MatrixXi innerF;
		if (!igl::read_triangle_mesh(resultArguments.outerPath, outerV, outerF))
		{
			throw std::runtime_error(
				"Failed to load sealed outer mesh: " + resultArguments.outerPath
			);
		}
		if (!igl::read_triangle_mesh(resultArguments.innerPath, innerV, innerF))
		{
			throw std::runtime_error(
				"Failed to load sealed cavity mesh: " + resultArguments.innerPath
			);
		}
		if (resultArguments.spinEnabled)
		{
			Eigen::MatrixXd designV;
			Eigen::MatrixXi designF;
			Eigen::MatrixXd materialV;
			Eigen::MatrixXi materialF;
			if (!igl::read_triangle_mesh(
					resultArguments.designPath, designV, designF))
			{
				throw std::runtime_error(
					"Failed to load Spin-It design domain: " +
					resultArguments.designPath
				);
			}
			if (!igl::read_triangle_mesh(
					resultArguments.materialPath, materialV, materialF))
			{
				throw std::runtime_error(
					"Failed to load Spin-It material boundary: " +
					resultArguments.materialPath
				);
			}
			const result_view::SpinResultViewMetadata metadata =
				result_view::loadSpinResultViewMetadata(
					resultArguments.diagnosticsPath
				);
			const result_view::SpinResultViewGeometry geometry =
				result_view::prepareSpinResultViewGeometry(
					outerV, outerF, designV, designF,
					materialV, materialF, innerV, innerF, metadata
				);
			origV = geometry.outerVertices;
			origF = geometry.outerFaces;
			state.V = geometry.outerVertices;
			state.F = geometry.outerFaces;
			state.designV = geometry.designVertices;
			state.designF = geometry.designFaces;
			state.materialV = geometry.materialVertices;
			state.materialF = geometry.materialFaces;
			state.innerV = geometry.cavityVertices;
			state.innerF = geometry.cavityFaces;
			state.planeV = geometry.targetPlaneVertices;
			state.planeF = geometry.targetPlaneFaces;
			state.spinMetadata = geometry.metadata;
			state.spinAxisDisplayLength = geometry.axisDisplayLength;
			state.planeCenter = state.spinMetadata.targetOrigin.transpose();
			spinResultViewMode = true;
		}
		else
		{
			const result_view::ResultViewGeometry geometry =
				result_view::prepareResultViewGeometry(
					outerV, outerF, innerV, innerF
				);
		origV = geometry.outerVertices;
		origF = geometry.outerFaces;
		state.V = geometry.outerVertices;
		state.F = geometry.outerFaces;
		state.innerV = geometry.innerVertices;
		state.innerF = geometry.innerFaces;
		state.planeV = geometry.groundVertices;
		state.planeF = geometry.groundFaces;
		state.planeCenter = state.planeV.colwise().mean();
		}
		state.gravity = DEFAULT_GRAVITY;
		state.yaw = DEFAULT_GRAVITY_YAW;
		state.pitch = DEFAULT_GRAVITY_PITCH;
		state.roll = DEFAULT_GRAVITY_ROLL;
		state.contactThreshold = DEFAULT_CONTACT_THRESHOLD;
		state.topplingAngle = DEFAULT_TOPPLING_ANGLE;
		state.minimumWallThickness = DEFAULT_MINIMUM_WALL_THICKNESS;
		state.stabilityMessage = "Read-only sealed result";
		if (!findLowestPointIdx(state.V, state.balancePointIdx))
		{
			throw std::runtime_error("Sealed outer mesh has no lowest point");
		}
		resultViewMode = true;
		resultWindowTitle = resultArguments.title;
		return;
	}

	// Load object
  	igl::read_triangle_mesh(argc > 1 ? argv[1] : DEFAULT_MESH_FILE, origV, origF);

	// Align object to axis planes
	Eigen::MatrixXd unalignedV = origV;
	alignToAxis(unalignedV, origV);

	// Initialize display V and F
	state.V = origV;
	state.F = origF;

	// Create plane from bottom of the object's bounding box
	Eigen::AlignedBox3d boundingBox;
	createAlignedBox(origV, boundingBox);
	state.planeV.resize(4,3);
	state.planeV <<
		boundingBox.corner(boundingBox.BottomLeftCeil).transpose(),
		boundingBox.corner(boundingBox.BottomRightCeil).transpose(),
		boundingBox.corner(boundingBox.BottomLeftFloor).transpose(),
		boundingBox.corner(boundingBox.BottomRightFloor).transpose();
	state.planeF.resize(2,3);
	state.planeF <<
		0, 1, 2,
		2, 1, 3;

	state.planeCenter = state.planeV.colwise().sum() / 4.0;

	// Gravity points down by default
	state.yaw = DEFAULT_GRAVITY_YAW;
	state.pitch = DEFAULT_GRAVITY_PITCH;
	state.roll = DEFAULT_GRAVITY_ROLL;
	state.contactThreshold = DEFAULT_CONTACT_THRESHOLD;
	state.topplingAngle = DEFAULT_TOPPLING_ANGLE;
	state.minimumWallThickness = DEFAULT_MINIMUM_WALL_THICKNESS;
	state.stabilityMessage = "Support target not constructed";
	updateGravityVector();
	
	// set lowest point as balance point
	if (!findLowestPointIdx(origV, state.balancePointIdx)) {
		std::cerr << "Failed to find index of lowest mesh point";
		exit(-1);
	}
}


void update_spin_visibility()
{
	viewer.data(mesh_data_id).show_faces = showSpinOuter;
	viewer.data(mesh_data_id).show_lines = showSpinOuter;
	viewer.data(design_data_mesh_id).show_faces = showSpinDesign;
	viewer.data(design_data_mesh_id).show_lines = showSpinDesign;
	viewer.data(material_data_mesh_id).show_faces = showSpinMaterial;
	viewer.data(material_data_mesh_id).show_lines = false;
	viewer.data(inner_data_mesh_id).show_faces = showSpinCavity;
	viewer.data(inner_data_mesh_id).show_lines = false;
	viewer.data(axes_data_id).show_lines = showSpinAxes;
	viewer.data(axes_data_id).show_overlay = showSpinAxes;
	viewer.data(plane_data_id).show_faces = showSpinTargetPlane;
}


void draw_spin_result_window()
{
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(365, 310), ImGuiCond_FirstUseEver);
	ImGui::Begin("Spin-It Faster artifact", nullptr,
		ImGuiWindowFlags_NoSavedSettings);
	ImGui::TextWrapped("%s", resultWindowTitle.c_str());
	ImGui::Text("Mode: %s", state.spinMetadata.mode.c_str());
	ImGui::Separator();
	bool changed = false;
	changed |= ImGui::Checkbox("Fixed outer shell", &showSpinOuter);
	changed |= ImGui::Checkbox("Design domain", &showSpinDesign);
	changed |= ImGui::Checkbox("Retained material", &showSpinMaterial);
	changed |= ImGui::Checkbox("Cavity / signed quadric", &showSpinCavity);
	changed |= ImGui::Checkbox("CoM, target and principal axes", &showSpinAxes);
	changed |= ImGui::Checkbox("Target/support plane", &showSpinTargetPlane);
	if (changed) update_spin_visibility();
	ImGui::Separator();
	ImGui::TextColored(ImVec4(0.95f, 0.68f, 0.12f, 1.0f),
		"Gold: retained design material");
	ImGui::TextColored(ImVec4(0.90f, 0.10f, 0.08f, 1.0f),
		"Red: cavity / quadric negative region");
	ImGui::TextColored(ImVec4(0.18f, 0.48f, 0.88f, 1.0f),
		"Blue: target axis; RGB: principal axes");
	ImGui::TextWrapped("Read-only sealed result; outer coordinates are unchanged.");
	ImGui::End();
}

void draw_result_window()
{
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(330, 170), ImGuiCond_FirstUseEver);
	ImGui::Begin(
		"Sealed carving result",
		nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);
	ImGui::TextWrapped("%s", resultWindowTitle.c_str());
	ImGui::Separator();
	ImGui::TextColored(ImVec4(0.90f, 0.10f, 0.08f, 1.0f),
		"Red: removed internal material / cavity");
	ImGui::TextColored(ImVec4(0.78f, 0.84f, 0.95f, 1.0f),
		"Gray: retained outer shell");
	ImGui::TextColored(ImVec4(0.50f, 0.66f, 0.36f, 1.0f),
		"Green: support plane");
	ImGui::TextWrapped("Read-only result; experiment posture is preserved.");
	ImGui::End();
}

// == Callback Functions ==
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int x, int y) 
{
	// Select balance point
	// - source: from deformation assignment starter code
	Eigen::RowVector3f last_mouse = Eigen::RowVector3f(
		viewer.current_mouse_x, viewer.core().viewport(3) - viewer.current_mouse_y, 0);
	if (state.selectBalancePoint)
	{
		// Find closest point on mesh to mouse position
		int fid;
		Eigen::Vector3f bary;
		if (igl::unproject_onto_mesh(
			last_mouse.head(2), 
			viewer.core().view,
			viewer.core().proj, 
			viewer.core().viewport, 
			state.V, state.F, 
			fid, bary))
		{
			long c;
			bary.maxCoeff(&c);
			state.balancePointIdx = state.F(fid, c);

			// Snap to closest vertex on hit face
			updateBalancePoint();
			return true;
		}
	}
	return false;
}

void draw_viewer_menu()
{
	// Draw parent menu content
	menu.draw_viewer_menu();
}

void draw_workflow_control_window()
{
	// Define next window position + size
	ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(250, 360), ImGuiCond_FirstUseEver);
	ImGui::Begin(
		"Make It Stand", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);

	// Posture-anchor selection
	ImGui::Text("1. Posture Anchor");
	if (ImGui::Checkbox("Select a posture anchor", &(state.selectBalancePoint)))
	{
		if (state.selectBalancePoint) // switching to balance point selection
		{
			// Set other checkboxes to false
			state.selectOrientation = false;
			state.isCarving = false;

			// Update balance point color
			updateBalancePoint();

			// Undo carving
			clearInnerMesh();
			viewer.data(mesh_data_id).show_faces = true;
		}
	}
	if (state.selectBalancePoint)
	{
		if (ImGui::Button("Update", ImVec2(-1, 0)))
		{
			updateMesh();
		}
	}

	// Orientation selection
	ImGui::Separator();
	ImGui::Text("2. Orientation");
	if (ImGui::Checkbox("Select an orientation", &(state.selectOrientation)))
	{
		if (state.selectOrientation) // switching to orientation selection
		{
			// Set other checkboxes to false
			state.selectBalancePoint = false;
			state.isCarving = false;

			// Undo carving
			clearInnerMesh();
			viewer.data(mesh_data_id).show_faces = true;

			// Align object to its balance point
			updateMesh();
		}
	}
	if (state.selectOrientation)
	{
		ImGui::SliderFloat("Yaw", &(state.yaw), 0.0, 360.0);
		if (ImGui::IsItemActive())
		{
			updateGravity();
			updateMesh();
		}
		ImGui::SliderFloat("Pitch", &(state.pitch), 0.0, 360.0);
		if (ImGui::IsItemActive())
		{
			updateGravity();
			updateMesh();
		}
		ImGui::SliderFloat("Roll", &(state.roll), 0.0, 360.0);
		if (ImGui::IsItemActive())
		{
			updateGravity();
			updateMesh();
		}
	}

	// Carving selection
	ImGui::Separator();
	ImGui::Text("3. Inner Carving");
	ImGui::InputFloat(
		"Planned flatten band",
		&(state.contactThreshold),
		0.001f,
		0.01f,
		"%.4f"
	);
	ImGui::SliderFloat(
		"Toppling angle",
		&(state.topplingAngle),
		0.0f,
		45.0f,
		"%.1f deg"
	);
	ImGui::InputFloat(
		"Minimum wall thickness",
		&(state.minimumWallThickness),
		0.001f,
		0.01f,
		"%.4f"
	);
	if (ImGui::Checkbox("Begin carving", &(state.isCarving)))
	{
		if (state.isCarving) // switching to carving
		{
			// Set other checkboxes to false
			state.selectBalancePoint = false;
			state.selectOrientation = false;

			// Align object to its balance point
			updateMesh();

			// Display setup: voxelize inner mesh, make outer mesh transparent
			viewer.data(mesh_data_id).show_faces = false;
			if (state.innerMesh != NULL)
			{
				delete state.innerMesh;
				state.innerMesh = NULL;
			}
			try
			{
				const Eigen::Vector3d worldGravity(0.0, -1.0, 0.0);
				const stability::SupportRegion support =
					stability::SupportRegion::fromMesh(
						state.V,
						worldGravity,
						static_cast<double>(state.contactThreshold)
					);

				const MassProperties massProperties =
					computeSurfaceMassProperties(state.V, state.F);
				Eigen::Vector3d centerOfMass;
				if (!massProperties.tryCenterOfMass(centerOfMass))
				{
					throw std::runtime_error("Could not compute a finite center of mass");
				}

				const stability::StandingTarget standingTarget =
					stability::makeStandingTarget(
						support,
						centerOfMass,
						static_cast<double>(state.topplingAngle) * PI / 180.0
					);

				if (!standingTarget.requestedMarginAttainable())
				{
					throw std::runtime_error(
						"Requested toppling margin exceeds the support inradius"
					);
				}

				BalanceObjectives objectives;
				objectives.push_back(standingTarget.objective());
				state.innerMesh = new InnerVoidMesh(
					state.V,
					state.F,
					objectives,
					static_cast<double>(state.minimumWallThickness)
				);
				state.stabilityMessage =
					state.contactThreshold > 0.0f ?
					"Target uses a hypothetical flatten band; mesh is not flattened in R3" :
					"Standing target constructed from actual lowest contacts";
				state.innerMesh->convertToMesh(state.innerV, state.innerF);
				updateInnerMesh();
			}
			catch (const std::exception& exception)
			{
				state.isCarving = false;
				state.stabilityMessage = exception.what();
				viewer.data(mesh_data_id).show_faces = true;
			}
		}
	}
	ImGui::TextWrapped("%s", state.stabilityMessage.c_str());
	if (state.isCarving)
	{
		if (ImGui::Button("Carve", ImVec2(-1, 0)))
		{
			if (!state.innerMesh->isOptimized())
			{
				Eigen::MatrixXd carvePlaneV;
				Eigen::MatrixXi carvePlaneF;
				Eigen::Vector3d com;
				state.innerMesh->carveInnerMeshStep(carvePlaneV, carvePlaneF, com);
				state.innerMesh->convertToMesh(state.innerV, state.innerF);
				updateInnerMesh();

				if (DEBUG)
				{
					viewer.data(carve_plane_data_id).clear();
					viewer.data(carve_plane_data_id).set_mesh(carvePlaneV, carvePlaneF);
					viewer.data(carve_plane_data_id).add_points(com.transpose(), orange);
				}
			}
		}
		if (ImGui::Button("Finish Carving", ImVec2(-1, 0)))
		{
			state.innerMesh->carveInnerMesh();
			state.innerMesh->convertToMesh(state.innerV, state.innerF);
			updateInnerMesh();
		}
		if (state.innerMesh->isOptimized())
		{
			ImGui::Text("Completed Carving");
		}
	}

	ImGui::End();
}

int main(int argc, char *argv[])
{
	try
	{
		initState(argc, argv);
	}
	catch (const std::exception& exception)
	{
		std::cerr << "make-it-stand viewer: " << exception.what() << '\n';
		return EXIT_FAILURE;
	}

	if (resultViewMode)
	{
		viewer.plugins.push_back(&menu);
		menu.callback_draw_viewer_menu = draw_viewer_menu;
		menu.callback_draw_custom_window = spinResultViewMode ?
			draw_spin_result_window : draw_result_window;

		if (spinResultViewMode)
		{
			plane_data_id = viewer.data().id;
			viewer.data(plane_data_id).set_mesh(state.planeV, state.planeF);
			viewer.data(plane_data_id).set_colors(resultGroundColor);
			viewer.data(plane_data_id).show_lines = false;

			viewer.append_mesh();
			inner_data_mesh_id = viewer.data().id;
			viewer.data(inner_data_mesh_id).set_mesh(state.innerV, state.innerF);
			viewer.data(inner_data_mesh_id).set_colors(resultVoidColor);

			viewer.append_mesh();
			material_data_mesh_id = viewer.data().id;
			viewer.data(material_data_mesh_id).set_mesh(
				state.materialV, state.materialF);
			viewer.data(material_data_mesh_id).set_colors(spinMaterialColor);

			viewer.append_mesh();
			design_data_mesh_id = viewer.data().id;
			viewer.data(design_data_mesh_id).set_mesh(
				state.designV, state.designF);
			viewer.data(design_data_mesh_id).set_colors(spinDesignColor);

			viewer.append_mesh();
			mesh_data_id = viewer.data().id;
			viewer.data(mesh_data_id).set_mesh(state.V, state.F);
			viewer.data(mesh_data_id).set_colors(resultShellColor);

			viewer.append_mesh();
			axes_data_id = viewer.data().id;
			const double length = state.spinAxisDisplayLength;
			const Eigen::RowVector3d center =
				state.spinMetadata.centerOfMass.transpose();
			viewer.data(axes_data_id).add_points(center, orange);
			viewer.data(axes_data_id).add_points(
				state.spinMetadata.targetOrigin.transpose(), blue);
			viewer.data(axes_data_id).add_edges(
				state.spinMetadata.targetOrigin.transpose(),
				(state.spinMetadata.targetOrigin +
				 length * state.spinMetadata.targetAxes.col(2)).transpose(),
				blue);
			const Eigen::RowVector3d axisColors[3] = {
				Eigen::RowVector3d(0.90, 0.15, 0.12),
				Eigen::RowVector3d(0.15, 0.75, 0.25),
				Eigen::RowVector3d(0.15, 0.35, 0.95)
			};
			for (int axis = 0; axis < 3; ++axis)
			{
				viewer.data(axes_data_id).add_edges(
					center,
					(state.spinMetadata.centerOfMass + length *
					 state.spinMetadata.principalAxes.col(axis)).transpose(),
					axisColors[axis]);
			}
			update_spin_visibility();
			viewer.selected_data_index = viewer.mesh_index(mesh_data_id);
			viewer.core().align_camera_center(state.V, state.F);
			return viewer.launch(true, false, resultWindowTitle, 1280, 800);
		}

		// Draw back-to-front so the opaque cavity remains visible through the
		// translucent, coordinate-preserving sealed outer shell.
		plane_data_id = viewer.data().id;
		viewer.data(plane_data_id).set_mesh(state.planeV, state.planeF);
		viewer.data(plane_data_id).set_colors(resultGroundColor);
		viewer.data(plane_data_id).show_lines = false;

		viewer.append_mesh();
		inner_data_mesh_id = viewer.data().id;
		viewer.data(inner_data_mesh_id).set_mesh(state.innerV, state.innerF);
		viewer.data(inner_data_mesh_id).set_colors(resultVoidColor);
		viewer.data(inner_data_mesh_id).show_faces = true;
		viewer.data(inner_data_mesh_id).show_lines = false;

		viewer.append_mesh();
		mesh_data_id = viewer.data().id;
		viewer.data(mesh_data_id).set_mesh(state.V, state.F);
		viewer.data(mesh_data_id).set_colors(resultShellColor);
		viewer.data(mesh_data_id).show_faces = true;
		viewer.data(mesh_data_id).show_lines = true;

		viewer.selected_data_index = viewer.mesh_index(mesh_data_id);
		viewer.core().align_camera_center(state.V, state.F);
		return viewer.launch(true, false, resultWindowTitle, 1280, 800);
	}

	// == GUI ==
	// menu
	viewer.plugins.push_back(&menu);
	menu.callback_draw_viewer_menu = draw_viewer_menu;
	menu.callback_draw_custom_window = draw_workflow_control_window;
	
	// callbacks
	viewer.callback_mouse_down = mouse_down;

	// == Display Objects ==
	// display main mesh 
	mesh_data_id = viewer.data().id;
	viewer.data(mesh_data_id).set_mesh(state.V, state.F);
	updateMesh();

	// display inner voxelized mesh (representing the empty space)
	viewer.append_mesh();
	inner_data_mesh_id = viewer.data().id;

	// display plane (representing the ground)
	viewer.append_mesh();
	plane_data_id = viewer.data().id;
	viewer.data(plane_data_id).set_mesh(state.planeV, state.planeF);

	// gravity vector (for debugging)
	if (DEBUG)
	{
		viewer.append_mesh();
		gravity_data_id = viewer.data().id;
		updateGravity();

		viewer.append_mesh();
		carve_plane_data_id = viewer.data().id;
	}

	// currently selected mesh is the input object
	viewer.selected_data_index = viewer.mesh_index(mesh_data_id);
	viewer.core().align_camera_center(state.planeV, state.planeF);

	viewer.launch();
}
