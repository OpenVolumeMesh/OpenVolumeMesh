#include <iostream>

#include <gtest/gtest.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include "unittests_common.hh"

using namespace OpenVolumeMesh;


TEST_F(TetrahedralMeshBase, SwapVerticesWithCommonEdge)
{
    generateTetrahedralMesh(mesh_);
    auto eh = HEH(0);
    auto v0 = mesh_.from_vertex_handle(eh);
    auto v1 = mesh_.to_vertex_handle(eh);
    mesh_.swap_vertex_indices(v0, v1);
    ASSERT_EQ(mesh_.from_vertex_handle(eh), v1);
    ASSERT_EQ(mesh_.to_vertex_handle(eh), v0);
}

