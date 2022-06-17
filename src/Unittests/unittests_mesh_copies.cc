#include "unittests_common.hh"

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <OpenVolumeMesh/Attribs/ColorAttrib.hh>

using namespace OpenVolumeMesh;
using OpenVolumeMesh::Geometry::Vec3d;

class MeshCopyTest : public testing::Test {
};

TEST_F(MeshCopyTest, CopyConstructor)
{
  using MeshT =  OpenVolumeMesh::GeometricPolyhedralMeshV3d;
  MeshT mesh;
  auto pos0 = Vec3d(-23.0, -42.0, -1.0);
  auto pos1 = Vec3d(99.0, -14.0, -7.0);
  auto pos2 = Vec3d(123, 456, 789);
  VertexHandle v0 = mesh.add_vertex(pos0);
  VertexHandle v1 = mesh.add_vertex(pos1);
  VertexHandle v2 = mesh.add_vertex(pos2);
  EdgeHandle e0 = mesh.add_edge(v0, v2);

  mesh.enable_deferred_deletion(true);
  mesh.delete_vertex(v1);
  PolyhedralMesh copy{mesh};

  ASSERT_EQ(copy.n_vertices(), mesh.n_vertices());
  ASSERT_EQ(mesh.is_deleted(v1), true);
  ASSERT_EQ(copy.is_deleted(v1), true);
  ASSERT_EQ(copy.n_logical_vertices(), mesh.n_logical_vertices());
  ASSERT_EQ(copy.n_edges(), mesh.n_edges());
  ASSERT_EQ(copy.n_faces(), mesh.n_faces());
  ASSERT_EQ(copy.n_cells(), mesh.n_cells());

  EXPECT_EQ(copy.vertex(v0), mesh.vertex(v0));
  EXPECT_EQ(copy.vertex(v1), mesh.vertex(v1));
  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(0)), mesh.from_vertex_handle(e0.halfedge_handle(0)));
  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(1)), mesh.from_vertex_handle(e0.halfedge_handle(1)));
  // TODO: test face topology
  // TODO: test cell topology

  VertexHandle orig_v3 = mesh.add_vertex(pos0+pos1);
  ASSERT_EQ(copy.n_vertices(), mesh.n_vertices()-1);
  mesh.set_vertex(v0, Vec3d(41));
  mesh.set_vertex(v1, Vec3d(55));
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);
  mesh.clear();
  ASSERT_EQ(copy.n_logical_vertices(), 2);
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);

  // TODO: test persistent properties
  // TODO: test tet and hex meshes
  // TODO: a deep comparison operator would help here. or save to (in-memory) file to compare binary representation?
}
