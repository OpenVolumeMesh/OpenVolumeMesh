#include <iostream>

#include <gtest/gtest.h>

#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include "unittests_common.hh"

using namespace OpenVolumeMesh;

/// Assert that the vertices 0, 1, 2, 3 are all present and correctly oriented
void ASSERT_VERTICES_0123(TetrahedralMesh& mesh, CellHandle ch, const std::vector<VertexHandle>& vhs)
{
    ASSERT_EQ(vhs.size(), 4);
    ASSERT_NE(std::find(vhs.begin(),vhs.end(),VH(0)), vhs.end());
    ASSERT_NE(std::find(vhs.begin(),vhs.end(),VH(1)), vhs.end());
    ASSERT_NE(std::find(vhs.begin(),vhs.end(),VH(2)), vhs.end());
    ASSERT_NE(std::find(vhs.begin(),vhs.end(),VH(3)), vhs.end());

    HalfFaceHandle hfh012 = mesh.find_halfface({vhs[0], vhs[1], vhs[2]});
    HalfFaceHandle hfh023 = mesh.find_halfface({vhs[0], vhs[2], vhs[3]});
    HalfFaceHandle hfh031 = mesh.find_halfface({vhs[0], vhs[3], vhs[1]});
    HalfFaceHandle hfh132 = mesh.find_halfface({vhs[1], vhs[3], vhs[2]});

    ASSERT_EQ(mesh.incident_cell(hfh012), ch);
    ASSERT_EQ(mesh.incident_cell(hfh023), ch);
    ASSERT_EQ(mesh.incident_cell(hfh031), ch);
    ASSERT_EQ(mesh.incident_cell(hfh132), ch);
}

TEST_F(TetrahedralMeshBase, GetCellVerticesHfh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);
    HalfFaceHandle hfh023 = mesh_.add_halfface(vh0, vh2, vh3);
    HalfFaceHandle hfh031 = mesh_.add_halfface(vh0, vh3, vh1);
    HalfFaceHandle hfh132 = mesh_.add_halfface(vh1, vh3, vh2);

    CellHandle ch = mesh_.add_cell({hfh012, hfh023, hfh031, hfh132});

    auto vhs = mesh_.get_cell_vertices(hfh012);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh3);

    vhs = mesh_.get_cell_vertices(hfh023);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh3);
    ASSERT_EQ(vhs[3], vh1);

    vhs = mesh_.get_cell_vertices(hfh031);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh3);
    ASSERT_EQ(vhs[2], vh1);
    ASSERT_EQ(vhs[3], vh2);

    vhs = mesh_.get_cell_vertices(hfh132);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh3);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh0);
}

TEST_F(TetrahedralMeshBase, GetCellVerticesCh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    CellHandle ch = mesh_.add_cell(vh0, vh1, vh2, vh3);

    auto vhs = mesh_.get_cell_vertices(ch);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh3);
}

TEST_F(TetrahedralMeshBase, GetCellVerticesChVh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    CellHandle ch = mesh_.add_cell(vh0, vh1, vh2, vh3);

    auto vhs = mesh_.get_cell_vertices(ch, vh0);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);

    vhs = mesh_.get_cell_vertices(ch, vh1);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh1);

    vhs = mesh_.get_cell_vertices(ch, vh2);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh2);

    vhs = mesh_.get_cell_vertices(ch, vh3);
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh3);
}

TEST_F(TetrahedralMeshBase, GetCellVerticesHfhHeh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);
    HalfFaceHandle hfh023 = mesh_.add_halfface(vh0, vh2, vh3);
    HalfFaceHandle hfh031 = mesh_.add_halfface(vh0, vh3, vh1);
    HalfFaceHandle hfh132 = mesh_.add_halfface(vh1, vh3, vh2);

    CellHandle ch = mesh_.add_cell({hfh012, hfh023, hfh031, hfh132});

    auto vhs = mesh_.get_cell_vertices(hfh012, mesh_.find_halfedge(vh0, vh1));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh3);

    vhs = mesh_.get_cell_vertices(hfh012, mesh_.find_halfedge(vh1, vh2));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh0);
    ASSERT_EQ(vhs[3], vh3);

    vhs = mesh_.get_cell_vertices(hfh012, mesh_.find_halfedge(vh2, vh0));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh2);
    ASSERT_EQ(vhs[1], vh0);
    ASSERT_EQ(vhs[2], vh1);
    ASSERT_EQ(vhs[3], vh3);

    vhs = mesh_.get_cell_vertices(hfh023, mesh_.find_halfedge(vh0, vh2));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh3);
    ASSERT_EQ(vhs[3], vh1);

    vhs = mesh_.get_cell_vertices(hfh023, mesh_.find_halfedge(vh2, vh3));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh2);
    ASSERT_EQ(vhs[1], vh3);
    ASSERT_EQ(vhs[2], vh0);
    ASSERT_EQ(vhs[3], vh1);

    vhs = mesh_.get_cell_vertices(hfh023, mesh_.find_halfedge(vh3, vh0));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh3);
    ASSERT_EQ(vhs[1], vh0);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh1);

    vhs = mesh_.get_cell_vertices(hfh031, mesh_.find_halfedge(vh0, vh3));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh3);
    ASSERT_EQ(vhs[2], vh1);
    ASSERT_EQ(vhs[3], vh2);

    vhs = mesh_.get_cell_vertices(hfh031, mesh_.find_halfedge(vh3, vh1));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh3);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh0);
    ASSERT_EQ(vhs[3], vh2);

    vhs = mesh_.get_cell_vertices(hfh031, mesh_.find_halfedge(vh1, vh0));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh0);
    ASSERT_EQ(vhs[2], vh3);
    ASSERT_EQ(vhs[3], vh2);

    vhs = mesh_.get_cell_vertices(hfh132, mesh_.find_halfedge(vh1, vh3));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh3);
    ASSERT_EQ(vhs[2], vh2);
    ASSERT_EQ(vhs[3], vh0);

    vhs = mesh_.get_cell_vertices(hfh132, mesh_.find_halfedge(vh3, vh2));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh3);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh1);
    ASSERT_EQ(vhs[3], vh0);

    vhs = mesh_.get_cell_vertices(hfh132, mesh_.find_halfedge(vh2, vh1));
    ASSERT_VERTICES_0123(mesh_, ch, vhs);
    ASSERT_EQ(vhs[0], vh2);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh3);
    ASSERT_EQ(vhs[3], vh0);
}

TEST_F(TetrahedralMeshBase, HalffaceOppositeVertex)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);
    HalfFaceHandle hfh023 = mesh_.add_halfface(vh0, vh2, vh3);
    HalfFaceHandle hfh031 = mesh_.add_halfface(vh0, vh3, vh1);
    HalfFaceHandle hfh132 = mesh_.add_halfface(vh1, vh3, vh2);

    mesh_.add_cell({hfh012, hfh023, hfh031, hfh132});

    ASSERT_EQ(mesh_.halfface_opposite_vertex(hfh012), vh3);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(hfh023), vh1);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(hfh031), vh2);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(hfh132), vh0);
}

TEST_F(TetrahedralMeshBase, VertexOppositeHalfface)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);
    HalfFaceHandle hfh023 = mesh_.add_halfface(vh0, vh2, vh3);
    HalfFaceHandle hfh031 = mesh_.add_halfface(vh0, vh3, vh1);
    HalfFaceHandle hfh132 = mesh_.add_halfface(vh1, vh3, vh2);

    CellHandle ch = mesh_.add_cell({hfh012, hfh023, hfh031, hfh132});

    ASSERT_EQ(mesh_.vertex_opposite_halfface(ch, vh0), hfh132);
    ASSERT_EQ(mesh_.vertex_opposite_halfface(ch, vh1), hfh023);
    ASSERT_EQ(mesh_.vertex_opposite_halfface(ch, vh2), hfh031);
    ASSERT_EQ(mesh_.vertex_opposite_halfface(ch, vh3), hfh012);
}

TEST_F(TetrahedralMeshBase, HalffaceOppositeVertexOppositeHalfface)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    CellHandle ch = mesh_.add_cell(vh0, vh1, vh2, vh3);

    ASSERT_EQ(mesh_.halfface_opposite_vertex(mesh_.vertex_opposite_halfface(ch, vh0)), vh0);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(mesh_.vertex_opposite_halfface(ch, vh1)), vh1);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(mesh_.vertex_opposite_halfface(ch, vh2)), vh2);
    ASSERT_EQ(mesh_.halfface_opposite_vertex(mesh_.vertex_opposite_halfface(ch, vh3)), vh3);
}
