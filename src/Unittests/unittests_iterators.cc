#include "unittests_common.hh"


using namespace OpenVolumeMesh;

TEST_F(HexahedralMeshBase, HexVertexIterTest) {

    generateHexahedralMesh(mesh_);

    HexVertexIter hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_TRUE(hv_it.valid());

    EXPECT_HANDLE_EQ(VertexHandle(0), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(1), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(2), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(3), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(4), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(7), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(6), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(5), *hv_it);
}

TEST_F(TetrahedralMeshBase, VertexVertexIteratorTest) {

    generateTetrahedralMesh(mesh_);

    {

        VertexVertexIter vv_it = mesh_.vv_iter(VertexHandle(0));

        EXPECT_TRUE(vv_it.valid());

        std::set<VertexHandle> onering;
        int valence = 0;

        while (vv_it.valid())
        {
          ++valence;
          onering.insert(*vv_it);
          ++vv_it;
        }

        // check that there have been three adjacent vertices
        EXPECT_EQ(3, valence);

        // check that no vertex was visited twice
        EXPECT_EQ(3u, onering.size());

        // check that no invalid vertex was adjacent
        EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }

    {

      std::set<VertexHandle> onering;
      int valence = 0;

      for (auto vh : mesh_.vertex_vertices(VertexHandle(0)))
      {
        ++valence;
        onering.insert(vh);
      }

      // check that there have been three adjacent vertices
      EXPECT_EQ(3, valence);

      // check that no vertex was visited twice
      EXPECT_EQ(3u, onering.size());

      // check that no invalid vertex was adjacent
      EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }

}

TEST_F(TetrahedralMeshBase, VertexFaceIteratorTest) {

    generateTetrahedralMesh(mesh_);

    {

        VertexFaceIter vf_it = mesh_.vf_iter(VertexHandle(0));

        EXPECT_TRUE(vf_it.valid());

        std::set<FaceHandle> incident_faces;
        int valence = 0;

        while (vf_it.valid())
        {
          ++valence;
          incident_faces.insert(*vf_it);
          ++vf_it;
        }

        // check that there have been three adjacent vertices
        EXPECT_EQ(3, valence);

        // check that no vertex was visited twice
        EXPECT_EQ(3u, incident_faces.size());

        // check that no invalid vertex was adjacent
        EXPECT_EQ(incident_faces.end(), std::find(incident_faces.begin(), incident_faces.end(), FaceHandle(-1)));

    }

    {

      std::set<VertexHandle> onering;
      int valence = 0;

      for (auto vh : mesh_.vertex_vertices(VertexHandle(0)))
      {
        ++valence;
        onering.insert(vh);
      }

      // check that there have been three adjacent vertices
      EXPECT_EQ(3, valence);

      // check that no vertex was visited twice
      EXPECT_EQ(3u, onering.size());

      // check that no invalid vertex was adjacent
      EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }

}

TEST_F(TetrahedralMeshBase, HalfFaceHalfEdgeIter)
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

    auto hfhe_it = mesh_.hfhe_iter(hfh012);
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh0, vh1));
    ++hfhe_it;
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh1, vh2));
    ++hfhe_it;
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh2, vh0));

    hfhe_it = mesh_.hfhe_iter(hfh012.opposite_handle());
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh0, vh2));
    ++hfhe_it;
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh2, vh1));
    ++hfhe_it;
    ASSERT_EQ(hfhe_it.cur_handle(), mesh_.find_halfedge(vh1, vh0));
}

TEST_F(TetrahedralMeshBase, HalfFaceVertexIter)
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

    auto hfv_it = mesh_.hfv_iter(hfh012);
    ASSERT_EQ(hfv_it.cur_handle(), vh0);
    ++hfv_it;
    ASSERT_EQ(hfv_it.cur_handle(), vh1);
    ++hfv_it;
    ASSERT_EQ(hfv_it.cur_handle(), vh2);

    hfv_it = mesh_.hfv_iter(hfh012.opposite_handle());
    ASSERT_EQ(hfv_it.cur_handle(), vh0);
    ++hfv_it;
    ASSERT_EQ(hfv_it.cur_handle(), vh2);
    ++hfv_it;
    ASSERT_EQ(hfv_it.cur_handle(), vh1);
}

TEST_F(TetrahedralMeshBase, GetHalfFaceVerticesHfh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);

    auto vhs = mesh_.get_halfface_vertices(hfh012);
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh2);

    vhs = mesh_.get_halfface_vertices(hfh012.opposite_handle());
    ASSERT_EQ(vhs[0], vh0);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh1);
}

TEST_F(TetrahedralMeshBase, GetHalfFaceVerticesHfhHeh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);

    auto vhs = mesh_.get_halfface_vertices(hfh012, mesh_.find_halfedge(vh1, vh2));
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh2);
    ASSERT_EQ(vhs[2], vh0);

    vhs = mesh_.get_halfface_vertices(hfh012.opposite_handle(), mesh_.find_halfedge(vh2, vh1));
    ASSERT_EQ(vhs[0], vh2);
    ASSERT_EQ(vhs[1], vh1);
    ASSERT_EQ(vhs[2], vh0);
}

TEST_F(TetrahedralMeshBase, GetHalfFaceVerticesHfhVh)
{
    VertexHandle vh0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VertexHandle vh1 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));
    VertexHandle vh2 = mesh_.add_vertex(Vec3d(1.0, 1.0, 0.0));

    HalfFaceHandle hfh012 = mesh_.add_halfface(vh0, vh1, vh2);

    auto vhs = mesh_.get_halfface_vertices(hfh012, vh2);
    ASSERT_EQ(vhs[0], vh2);
    ASSERT_EQ(vhs[1], vh0);
    ASSERT_EQ(vhs[2], vh1);

    vhs = mesh_.get_halfface_vertices(hfh012.opposite_handle(), vh1);
    ASSERT_EQ(vhs[0], vh1);
    ASSERT_EQ(vhs[1], vh0);
    ASSERT_EQ(vhs[2], vh2);
}

TEST_F(HexahedralMeshBase, RangeForTest) {
    // no EXPECTs here, if it compiles, it'll work.
    generateHexahedralMesh(mesh_);
    VertexHandle _dummy; // use vh to avoid compiler warnings
    for (const auto& vh: mesh_.vertices()) { _dummy = vh;}
    const auto& constref = mesh_;
    for (const auto& vh: constref.vertices()) { _dummy = vh;}
}
