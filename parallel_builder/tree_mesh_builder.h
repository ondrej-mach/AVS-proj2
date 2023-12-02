/**
 * @file    tree_mesh_builder.h
 *
 * @author  Ondrej Mach <xmacho12@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    29. 11. 2023
 **/

#ifndef TREE_MESH_BUILDER_H
#define TREE_MESH_BUILDER_H

#include <vector>
#include <omp.h>
#include "base_mesh_builder.h"

class TreeMeshBuilder : public BaseMeshBuilder
{
public:
    TreeMeshBuilder(unsigned gridEdgeSize);

protected:
    unsigned marchCubes(const ParametricScalarField &field);
    bool nodeEmpty(const ParametricScalarField &field, const Vec3_t<float> &origin, const float size);
    void computeDivision(const ParametricScalarField &field,
                         const Vec3_t<float> origin,
                         const unsigned size,
                         unsigned &triangleCount
                         );
    float evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field);
    void emitTriangle(const Triangle_t &triangle);
    const Triangle_t *getTrianglesArray() const { return mTriangles.data(); }

    const unsigned divisionCutoff = 2;
    std::vector<Triangle_t> mTriangles; ///< Temporary array of triangles
};

#endif // TREE_MESH_BUILDER_H
