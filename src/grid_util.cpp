#include "grid_util.h"

#include <igl/bounding_box.h>
#include <igl/voxel_grid.h>
#include <igl/grid.h>
#include <igl/signed_distance.h>
#include <igl/rotation_matrix_from_directions.h>

#include <cmath>
#include <limits>
#include <stdexcept>

void createAlignedBox(const Eigen::MatrixXd &V, Eigen::AlignedBox3d &box) {
	for (int i = 0; i < V.rows(); i++) {
		box.extend(V.row(i).transpose());
	}
}

void transformVertices(const Eigen::MatrixXd &V, const Eigen::AlignedBox3d &newOrientation, bool keepSize, Eigen::MatrixXd &Vout) {
	Vout.resize(V.rows(), 3);

	Eigen::Vector3d corner1 = newOrientation.corner(newOrientation.BottomLeftFloor);
	Eigen::Vector3d corner2 = newOrientation.corner(newOrientation.TopRightCeil);

	// Compute diagonals of input and output grids
	Eigen::AlignedBox<double, 3> inBox;
	createAlignedBox(V, inBox);
	Eigen::Vector3d inBLF = inBox.corner(inBox.BottomLeftFloor);
	Eigen::Vector3d inTRC = inBox.corner(inBox.TopRightCeil);
	Eigen::Vector3d inDiag = inBLF - inTRC;

	Eigen::Vector3d outDiag = corner1 - corner2;

	// Scale V to desired size
	Eigen::Matrix3d scale = Eigen::Matrix3d::Identity();
	if (!keepSize) {
		scale.diagonal() << outDiag.cwiseQuotient(inDiag).array();
	}
	Vout = V * scale;

	// Rotate
	inDiag = scale * inBLF - scale * inTRC;
	Eigen::Matrix3d R = igl::rotation_matrix_from_directions(inDiag.normalized(), outDiag.normalized());
	Vout *= R.transpose();

	// Move V to desired location
	Eigen::Vector3d translate = corner1 - R * (scale * inBLF);
	for (int i = 0; i < Vout.rows(); i++) {
		Vout.row(i) += translate.transpose();
	}
}

void transformVertices(const Eigen::MatrixXd& V, const Eigen::AlignedBox3d& newOrientation, Eigen::MatrixXd& Vout) {
	transformVertices(V, newOrientation, false, Vout);
}

void createVoxelGrid(
	const Eigen::MatrixXd &V, 
	Eigen::MatrixXd &centers, 
	Eigen::MatrixXd &corners, 
	Eigen::Vector3i &dimensions, 
	Eigen::Vector3d &voxelSize,
	int longestAxisResolution
) {
	if (longestAxisResolution < 3) {
		throw std::invalid_argument("voxel resolution must be at least three");
	}
	if (V.rows() == 0 || V.cols() != 3 || !V.allFinite()) {
		throw std::invalid_argument("voxel grid vertices must be a nonempty finite 3D array");
	}
	// Compute bounding box
	Eigen::MatrixXd BV, BF;
	igl::bounding_box(V, BV, BF);
	Eigen::AlignedBox3d boundBox;
	for (int i = 0; i < V.rows(); i++) {
		boundBox.extend(V.row(i).transpose());
	}

	// Create voxel grid (defined by center vertices of voxels)
	Eigen::RowVector3i dimOut; // column vector of dimentions
	igl::voxel_grid(boundBox, longestAxisResolution, 0, centers, dimOut);
	dimensions = dimOut.transpose();
	if ((dimensions.array() < 3).any()) {
		throw std::invalid_argument("voxel grid requires at least three cells on every axis");
	}
	Eigen::AlignedBox3d centersBox;
	createAlignedBox(centers, centersBox);

	Eigen::Vector3i one(1.0, 1.0, 1.0);
	voxelSize = centersBox.sizes().cwiseQuotient((dimensions - one).cast<double>());
	if (!voxelSize.allFinite() || (voxelSize.array() <= 0.0).any()) {
		throw std::invalid_argument("voxel grid produced invalid cell dimensions");
	}

	// Compute corner vertices of voxel grid
	Eigen::MatrixXd defaultGrid;
	Eigen::Vector3d dimCorners = (dimensions + one).cast<double>();
	igl::grid(dimCorners, defaultGrid);

	// scale default grid to voxel grid
	Eigen::Vector3d outBLF = centersBox.corner(centersBox.BottomLeftFloor) - (voxelSize / 2.0);
	Eigen::Vector3d outTRC = centersBox.corner(centersBox.TopRightCeil) + (voxelSize / 2.0);
	Eigen::AlignedBox3d outBox(outBLF, outTRC);
	transformVertices(defaultGrid, outBox, corners);
}

void createUnitCubeVoxelGrid(
	Eigen::MatrixXd& centers,
	Eigen::MatrixXd& corners,
	Eigen::Vector3i& dimensions,
	Eigen::Vector3d& voxelSize,
	int resolution
) {
	if (resolution < 3) {
		throw std::invalid_argument("voxel resolution must be at least three");
	}

	const long long cornerResolution =
		static_cast<long long>(resolution) + 1LL;
	const long long maximumCount =
		static_cast<long long>(std::numeric_limits<int>::max());
	auto checkedCubeCount = [maximumCount](long long side) {
		if (
			side <= 0 ||
			side > maximumCount / side ||
			side * side > maximumCount / side
		) {
			throw std::overflow_error(
				"unit-cube voxel grid exceeds integer indexing"
			);
		}
		return side * side * side;
	};
	const long long cellCount = checkedCubeCount(resolution);
	const long long cornerCount = checkedCubeCount(cornerResolution);
	if (cornerResolution > static_cast<long long>(std::numeric_limits<int>::max())) {
		throw std::overflow_error("unit-cube voxel grid exceeds integer indexing");
	}

	dimensions = Eigen::Vector3i::Constant(resolution);
	voxelSize = Eigen::Vector3d::Constant(
		1.0 / static_cast<double>(resolution)
	);
	centers.resize(static_cast<Eigen::Index>(cellCount), 3);
	for (int z = 0; z < resolution; ++z) {
		for (int y = 0; y < resolution; ++y) {
			for (int x = 0; x < resolution; ++x) {
				const Eigen::Index row = static_cast<Eigen::Index>(
					x + resolution * (y + resolution * z)
				);
				centers.row(row) <<
					(static_cast<double>(x) + 0.5) * voxelSize(0),
					(static_cast<double>(y) + 0.5) * voxelSize(1),
					(static_cast<double>(z) + 0.5) * voxelSize(2);
			}
		}
	}

	const int nodes = resolution + 1;
	corners.resize(static_cast<Eigen::Index>(cornerCount), 3);
	for (int z = 0; z < nodes; ++z) {
		for (int y = 0; y < nodes; ++y) {
			for (int x = 0; x < nodes; ++x) {
				const Eigen::Index row = static_cast<Eigen::Index>(
					x + nodes * (y + nodes * z)
				);
				corners.row(row) <<
					(x == resolution
						? 1.0
						: static_cast<double>(x) * voxelSize(0)),
					(y == resolution
						? 1.0
						: static_cast<double>(y) * voxelSize(1)),
					(z == resolution
						? 1.0
						: static_cast<double>(z) * voxelSize(2));
			}
		}
	}
}

void alignToAxis(const Eigen::MatrixXd& V, Eigen::MatrixXd& alignedV) {
	Eigen::AlignedBox3d inputBox;
	createAlignedBox(V, inputBox);

	Eigen::AlignedBox3d alignedBox;
	Eigen::MatrixXd alignedGrid;
	Eigen::Vector3d dim = Eigen::Vector3d::Constant(10.0);
	igl::grid(dim, alignedGrid);
	createAlignedBox(alignedGrid, alignedBox);

	transformVertices(V, alignedBox, true, alignedV);
}
