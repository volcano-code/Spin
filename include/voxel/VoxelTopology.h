#pragma once


class VoxelGrid;


/** Six-neighbor component counts for discrete material and carved void. */
struct VoxelTopologyReport
{
    int materialComponentCount;
    int carvedVoidComponentCount;
    int materialVoxelCount;
    int carvedVoidVoxelCount;

    VoxelTopologyReport();
};


/**
 * Count components on the interior object-voxel domain.
 *
 * Material means an inside voxel whose fill remains true. Carved void means
 * an eligible interior voxel whose fill is false. Empty exterior grid cells
 * do not participate.
 */
VoxelTopologyReport analyzeVoxelTopology(
    const VoxelGrid& voxelGrid
);
