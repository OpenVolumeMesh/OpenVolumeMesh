#include "unittests_common.hh"

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

using namespace OpenVolumeMesh;
using OpenVolumeMesh::Geometry::Vec3d;

class GarbageCollectorTest : public testing::Test {
public:
  using MeshT =  OpenVolumeMesh::TopologicTetrahedralMesh;
  MeshT mesh;

  void SetUp() override
  {
      VertexHandle v0 = mesh.add_vertex();
      VertexHandle v1 = mesh.add_vertex();
      VertexHandle v2 = mesh.add_vertex();
      VertexHandle v3 = mesh.add_vertex();

      FaceHandle f0 = mesh.add_face({v0, v1, v2});
      FaceHandle f1 = mesh.add_face({v0, v1, v3});
      FaceHandle f2 = mesh.add_face({v0, v2, v3});
      mesh.add_cell(v0, v1, v2, v3);

      mesh.enable_deferred_deletion(true);
  }

};

TEST_F(GarbageCollectorTest, n_entities)
{

  EXPECT_EQ(mesh.n_vertices(), 4);
  EXPECT_EQ(mesh.n_edges(), 6);
  EXPECT_EQ(mesh.n_halfedges(), 12);
  EXPECT_EQ(mesh.n_faces(), 4);
  EXPECT_EQ(mesh.n_halffaces(), 8);
  EXPECT_EQ(mesh.n_cells(), 1);

  EXPECT_EQ(mesh.n_logical_vertices(), 4);
  EXPECT_EQ(mesh.n_logical_edges(), 6);
  EXPECT_EQ(mesh.n_logical_halfedges(), 12);
  EXPECT_EQ(mesh.n_logical_faces(), 4);
  EXPECT_EQ(mesh.n_logical_halffaces(), 8);
  EXPECT_EQ(mesh.n_logical_cells(), 1);

  mesh.delete_vertex(VH(0));

  EXPECT_EQ(mesh.n_vertices(), 4);
  EXPECT_EQ(mesh.n_edges(), 6);
  EXPECT_EQ(mesh.n_halfedges(), 12);
  EXPECT_EQ(mesh.n_faces(), 4);
  EXPECT_EQ(mesh.n_halffaces(), 8);
  EXPECT_EQ(mesh.n_cells(), 1);

  EXPECT_EQ(mesh.n_logical_vertices(), 3);
  EXPECT_EQ(mesh.n_logical_edges(), 3);
  EXPECT_EQ(mesh.n_logical_halfedges(), 6);
  EXPECT_EQ(mesh.n_logical_faces(), 1);
  EXPECT_EQ(mesh.n_logical_halffaces(), 2);
  EXPECT_EQ(mesh.n_logical_cells(), 0);

  mesh.collect_garbage();

  EXPECT_EQ(mesh.n_vertices(), 3);
  EXPECT_EQ(mesh.n_edges(), 3);
  EXPECT_EQ(mesh.n_halfedges(), 6);
  EXPECT_EQ(mesh.n_faces(), 1);
  EXPECT_EQ(mesh.n_halffaces(), 2);
  EXPECT_EQ(mesh.n_cells(), 0);

  EXPECT_EQ(mesh.n_logical_vertices(), 3);
  EXPECT_EQ(mesh.n_logical_edges(), 3);
  EXPECT_EQ(mesh.n_logical_halfedges(), 6);
  EXPECT_EQ(mesh.n_logical_faces(), 1);
  EXPECT_EQ(mesh.n_logical_halffaces(), 2);
  EXPECT_EQ(mesh.n_logical_cells(), 0);

}

TEST_F(GarbageCollectorTest, n_edges)
{
}
