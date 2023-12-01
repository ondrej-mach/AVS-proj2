/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  Ondrej Mach <xmacho12@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    29. 11. 2023
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}


void TreeMeshBuilder::computeDivision(
    const ParametricScalarField &field,
    const Vec3_t<float> origin,
    const unsigned size,
    unsigned &triangleCount)
{
    if (size > divisionCutoff) {
        unsigned childSize = size / 2; // We assume a power of 2

        for (auto normOffset: sc_vertexNormPos) {
            Vec3_t<float> childOrigin(
                origin.x + normOffset.x * childSize,
                origin.y + normOffset.y * childSize,
                origin.z + normOffset.z * childSize
            );


            if (!nodeEmpty(field, childOrigin, childSize)) {
                #pragma omp task
                computeDivision(field, childOrigin, childSize, triangleCount);
            }
        }

    } else {
        #pragma omp atomic
        triangleCount += buildCube(origin, field);
    }
}


bool TreeMeshBuilder::nodeEmpty(const ParametricScalarField &field, const Vec3_t<float> &origin, const float size)
{
	Vec3_t<float> center(
                (origin.x + size/2.f) * mGridResolution,
                (origin.y + size/2.f) * mGridResolution,
                (origin.z + size/2.f) * mGridResolution
            );

    float fieldValue = evaluateFieldAt(center, field);

	return fieldValue > mIsoLevel + sqrt(3)/2 * size * mGridResolution;
}


unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.

    unsigned totalTriangles;
    Vec3_t<float> origin(0, 0, 0);

    #pragma omp parallel
    #pragma omp single
    computeDivision(field, origin, mGridSize, totalTriangles);
    #pragma omp taskwait

    return totalTriangles;
}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // NOTE: This method is called from "buildCube(...)"!

    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    #pragma omp simd reduction(min:value)
    for (unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    // NOTE: This method is called from "buildCube(...)"!

    // Store generated triangle into vector (array) of generated triangles.
    // The pointer to data in this array is return by "getTrianglesArray(...)" call
    // after "marchCubes(...)" call ends.
    #pragma omp critical(mTriangles)
    mTriangles.push_back(triangle);
}
