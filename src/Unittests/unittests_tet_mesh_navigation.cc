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

TEST_F(TetrahedralMeshBase, GetCellVerticesConsistencyV33)
{
    auto get_cell_vertices_hfh_heh_v33 = [&](HalfFaceHandle hfh, HalfEdgeHandle heh) -> std::vector<VertexHandle>
    {
        std::vector<VertexHandle> vertices;

        // add vertices of halfface
        for (unsigned int i = 0; i < 3; ++i)
        {
            auto e = mesh_.halfedge(heh);
            vertices.push_back(e.from_vertex());
            heh = mesh_.next_halfedge_in_halfface(heh, hfh);
        }

        auto c = mesh_.cell(mesh_.incident_cell(hfh));
        HalfFaceHandle otherHfh = c.halffaces()[0];
        if (otherHfh == hfh)
            otherHfh = c.halffaces()[1];

        auto otherF = mesh_.halfface(otherHfh);

        for (unsigned int i = 0; i < otherF.halfedges().size(); ++i)
        {
            HalfEdgeHandle he = otherF.halfedges()[i];
            auto e = mesh_.halfedge(he);
            if (std::find(vertices.begin(), vertices.end(), e.to_vertex()) == vertices.end())
            {
                vertices.push_back(e.to_vertex());
                return vertices;
            }
        }

        return vertices;
    };

    auto get_cell_vertices_ch_vh_v33 = [&](CellHandle ch, VertexHandle vh) -> std::vector<VertexHandle>
    {
        HalfFaceHandle hfh = mesh_.cell(ch).halffaces()[0];
        auto f = mesh_.halfface(hfh);
        HalfEdgeHandle heh;
        for (unsigned int i = 0; i < 3; ++i)
        {
            auto e = mesh_.halfedge(f.halfedges()[i]);
            if (e.from_vertex() == vh)
            {
                heh = f.halfedges()[i];
                break;
            }
        }
        if (!heh.is_valid())
        {
            hfh = mesh_.adjacent_halfface_in_cell(hfh, f.halfedges()[0]);
            heh = mesh_.prev_halfedge_in_halfface(mesh_.opposite_halfedge_handle(f.halfedges()[0]), hfh);
        }

        return get_cell_vertices_hfh_heh_v33(hfh,heh);
    };

    auto get_cell_vertices_hfh_v33 = [&](HalfFaceHandle hfh) -> std::vector<VertexHandle>
    {
        return get_cell_vertices_hfh_heh_v33(hfh, mesh_.halfface(hfh).halfedges().front());
    };

    auto get_cell_vertices_ch_v33 = [&](CellHandle ch) -> std::vector<VertexHandle>
    {
        return get_cell_vertices_hfh_v33(mesh_.cell(ch).halffaces().front());
    };

    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));
    VertexHandle vh3 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));

    std::vector<VertexHandle> dummy1 = {vh0, vh1, vh2, vh3};
    std::vector<VertexHandle> dummy2 = {vh0, vh1, vh2, vh3};
    std::vector<VertexHandle> dummy3 = {vh0, vh1, vh3, vh2};
    ASSERT_EQ(dummy1, dummy2);
    ASSERT_NE(dummy1, dummy3);

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);
    HalfFaceHandle hfh023 = mesh_.add_halfface(vh0, vh2, vh3);
    HalfFaceHandle hfh031 = mesh_.add_halfface(vh0, vh3, vh1);
    HalfFaceHandle hfh132 = mesh_.add_halfface(vh1, vh3, vh2);

    CellHandle ch = mesh_.add_cell({hfh012, hfh023, hfh031, hfh132});

    ASSERT_EQ(get_cell_vertices_ch_v33(ch), mesh_.get_cell_vertices(ch));

    ASSERT_EQ(get_cell_vertices_hfh_v33(hfh012), mesh_.get_cell_vertices(hfh012));
    ASSERT_EQ(get_cell_vertices_hfh_v33(hfh023), mesh_.get_cell_vertices(hfh023));
    ASSERT_EQ(get_cell_vertices_hfh_v33(hfh031), mesh_.get_cell_vertices(hfh031));
    ASSERT_EQ(get_cell_vertices_hfh_v33(hfh132), mesh_.get_cell_vertices(hfh132));

    ASSERT_EQ(get_cell_vertices_ch_vh_v33(ch, vh0), mesh_.get_cell_vertices(ch, vh0));
    ASSERT_EQ(get_cell_vertices_ch_vh_v33(ch, vh1), mesh_.get_cell_vertices(ch, vh1));
    ASSERT_EQ(get_cell_vertices_ch_vh_v33(ch, vh2), mesh_.get_cell_vertices(ch, vh2));
    ASSERT_EQ(get_cell_vertices_ch_vh_v33(ch, vh3), mesh_.get_cell_vertices(ch, vh3));

    HalfEdgeHandle heh01 = mesh_.find_halfedge(vh0, vh1);
    HalfEdgeHandle heh12 = mesh_.find_halfedge(vh1, vh2);
    HalfEdgeHandle heh20 = mesh_.find_halfedge(vh2, vh0);

    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh012, heh01), mesh_.get_cell_vertices(hfh012, heh01));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh012, heh12), mesh_.get_cell_vertices(hfh012, heh12));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh012, heh20), mesh_.get_cell_vertices(hfh012, heh20));

    HalfEdgeHandle heh02 = mesh_.find_halfedge(vh0, vh2);
    HalfEdgeHandle heh23 = mesh_.find_halfedge(vh2, vh3);
    HalfEdgeHandle heh30 = mesh_.find_halfedge(vh3, vh0);

    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh023, heh02), mesh_.get_cell_vertices(hfh023, heh02));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh023, heh23), mesh_.get_cell_vertices(hfh023, heh23));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh023, heh30), mesh_.get_cell_vertices(hfh023, heh30));

    HalfEdgeHandle heh03 = heh30.opposite_handle();
    HalfEdgeHandle heh31 = mesh_.find_halfedge(vh3, vh1);
    HalfEdgeHandle heh10 = heh01.opposite_handle();

    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh031, heh03), mesh_.get_cell_vertices(hfh031, heh03));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh031, heh31), mesh_.get_cell_vertices(hfh031, heh31));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh031, heh10), mesh_.get_cell_vertices(hfh031, heh10));

    HalfEdgeHandle heh13 = heh31.opposite_handle();
    HalfEdgeHandle heh32 = heh23.opposite_handle();
    HalfEdgeHandle heh21 = heh12.opposite_handle();

    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh132, heh13), mesh_.get_cell_vertices(hfh132, heh13));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh132, heh32), mesh_.get_cell_vertices(hfh132, heh32));
    ASSERT_EQ(get_cell_vertices_hfh_heh_v33(hfh132, heh21), mesh_.get_cell_vertices(hfh132, heh21));
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
