#include "unittests_common.hh"

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <OpenVolumeMesh/Attribs/ColorAttrib.hh>

using namespace OpenVolumeMesh;
using OpenVolumeMesh::Geometry::Vec3d;

TEST_F(PolyhedralMeshBase, CopyConstructor)
{
  auto pos0 = Vec3d(-23.0, -42.0, -1.0);
  auto pos1 = Vec3d(99.0, -14.0, -7.0);
  VertexHandle v0 = mesh_.add_vertex(pos0);
  VertexHandle v1 = mesh_.add_vertex(pos1);
  EdgeHandle e0 = mesh_.add_edge(v0, v1);

  PolyhedralMesh copy{mesh_};

  ASSERT_EQ(copy.n_vertices(), mesh_.n_vertices());
  ASSERT_EQ(copy.n_edges(), mesh_.n_edges());
  ASSERT_EQ(copy.n_faces(), mesh_.n_faces());
  ASSERT_EQ(copy.n_cells(), mesh_.n_cells());

  EXPECT_EQ(copy.vertex(v0), mesh_.vertex(v0));
  EXPECT_EQ(copy.vertex(v1), mesh_.vertex(v1));
  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(0)), mesh_.from_vertex_handle(e0.halfedge_handle(0)));
  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(1)), mesh_.from_vertex_handle(e0.halfedge_handle(1)));
  // TODO: test face topology
  // TODO: test cell topology

  VertexHandle orig_v2 = mesh_.add_vertex(pos0+pos1);
  ASSERT_EQ(copy.n_vertices(), mesh_.n_vertices()-1);
  mesh_.set_vertex(v0, Vec3d(41));
  mesh_.set_vertex(v1, Vec3d(55));
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);
  mesh_.clear();
  ASSERT_EQ(copy.n_vertices(), 2);
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);

  // TODO: test persistent properties
  // TODO: test tet and hex meshes
  // TODO: a deep comparison operator would help here. or save to (in-memory) file to compare binary representation?
}
