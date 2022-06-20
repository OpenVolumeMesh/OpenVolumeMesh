#include "unittests_common.hh"

#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

// TODO: parameterize tests on mesh type

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

  const std::string propname = "somename";
  auto mesh_prop = *mesh.create_persistent_property<int, Entity::Vertex>(propname);
  mesh_prop[v0] = 92;
  mesh_prop[v2] = 94;
  mesh.enable_deferred_deletion(true);
  mesh.delete_vertex(v1);


  PolyhedralMesh copy{mesh};

  PolyhedralMesh assigned;;
  assigned.add_vertex(Vec3d(-1,-2,-3));
  auto oldprop_pers = *assigned.create_persistent_property<int, Entity::Vertex>("pers");
  auto oldprop_pers2 = *assigned.create_persistent_property<int, Entity::Vertex>("pers2");
  auto oldprop_pers3 = *assigned.create_persistent_property<int, Entity::Vertex>("pers3");
  auto oldprop_shared = *assigned.create_shared_property<int, Entity::Vertex>("shared");
  auto oldprop_private = assigned.create_private_property<int, Entity::Vertex>("private");
  EXPECT_EQ(assigned.n_persistent_props<Entity::Vertex>(), 3);
  EXPECT_EQ(assigned.n_props<Entity::Vertex>(), 6); // count internal geometry prop
  EXPECT_EQ(oldprop_pers.size(), 1);
  EXPECT_EQ(oldprop_shared.size(), 1);
  EXPECT_EQ(oldprop_private.size(), 1);

  assigned = mesh;

  EXPECT_EQ(mesh.n_persistent_props<Entity::Vertex>(), 1);
  EXPECT_EQ(assigned.n_persistent_props<Entity::Vertex>(), 1);
  // existing, anonymized props should have a matching size for the new assigned contents:
  EXPECT_EQ(oldprop_pers.size(), 3);
  EXPECT_EQ(oldprop_shared.size(), 3);
  EXPECT_EQ(oldprop_private.size(), 3);

  EXPECT_EQ(mesh.n_vertices(), 3);
  EXPECT_EQ(copy.n_vertices(), 3);
  EXPECT_EQ(assigned.n_vertices(), 3);

  EXPECT_EQ(mesh.n_logical_vertices(), 2);
  EXPECT_EQ(copy.n_logical_vertices(), 2);
  EXPECT_EQ(assigned.n_logical_vertices(), 2);

  EXPECT_EQ(mesh.is_deleted(v1), true);
  EXPECT_EQ(copy.is_deleted(v1), true);
  EXPECT_EQ(assigned.is_deleted(v1), true);

  EXPECT_EQ(copy.n_edges(), mesh.n_edges());
  EXPECT_EQ(copy.n_faces(), mesh.n_faces());
  EXPECT_EQ(copy.n_cells(), mesh.n_cells());

  EXPECT_EQ(assigned.n_edges(), mesh.n_edges());
  EXPECT_EQ(assigned.n_faces(), mesh.n_faces());
  EXPECT_EQ(assigned.n_cells(), mesh.n_cells());

  EXPECT_EQ(copy.vertex(v0), mesh.vertex(v0));
  EXPECT_EQ(copy.vertex(v1), mesh.vertex(v1));
  EXPECT_EQ(assigned.vertex(v0), mesh.vertex(v0));
  EXPECT_EQ(assigned.vertex(v1), mesh.vertex(v1));

  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(0)), mesh.from_vertex_handle(e0.halfedge_handle(0)));
  EXPECT_EQ(copy.from_vertex_handle(e0.halfedge_handle(1)), mesh.from_vertex_handle(e0.halfedge_handle(1)));
  EXPECT_EQ(assigned.from_vertex_handle(e0.halfedge_handle(0)), mesh.from_vertex_handle(e0.halfedge_handle(0)));
  EXPECT_EQ(assigned.from_vertex_handle(e0.halfedge_handle(1)), mesh.from_vertex_handle(e0.halfedge_handle(1)));
  // TODO: test face topology
  // TODO: test cell topology

  VertexHandle orig_v3 = mesh.add_vertex(pos0+pos1);

  EXPECT_EQ(copy.n_vertices(), 3);
  EXPECT_EQ(assigned.n_vertices(), 3);
  EXPECT_EQ(copy.n_logical_vertices(), 2);
  EXPECT_EQ(assigned.n_logical_vertices(), 2);

  mesh.set_vertex(v0, Vec3d(41));
  mesh.set_vertex(v1, Vec3d(55));
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);
  EXPECT_EQ(assigned.vertex(v0), pos0);
  EXPECT_EQ(assigned.vertex(v1), pos1);
  mesh.clear(false);
  //mesh.clear(true);
  EXPECT_EQ(copy.n_logical_vertices(), 2);
  EXPECT_EQ(copy.vertex(v0), pos0);
  EXPECT_EQ(copy.vertex(v1), pos1);
  EXPECT_EQ(assigned.n_logical_vertices(), 2);
  EXPECT_EQ(assigned.vertex(v0), pos0);
  EXPECT_EQ(assigned.vertex(v1), pos1);


  auto copy_prop = *copy.get_property<int, Entity::Vertex>(propname);
  auto assigned_prop = *copy.get_property<int, Entity::Vertex>(propname);
  EXPECT_EQ(mesh_prop.size(), mesh.n_vertices());
  EXPECT_EQ(copy_prop.size(), copy.n_vertices());
  EXPECT_EQ(assigned_prop.size(), copy.n_vertices());

  EXPECT_EQ(copy_prop[v0], 92);
  EXPECT_EQ(copy_prop[v2], 94);
  EXPECT_EQ(assigned_prop[v0], 92);
  EXPECT_EQ(assigned_prop[v2], 94);
  // TODO: test tet and hex meshes
  // TODO: a deep comparison operator would help here. or save to (in-memory) file to compare binary representation?
}

TEST_F(MeshCopyTest, AssignMixedTopologies)
{
    PolyhedralMesh poly;
    TetrahedralMesh tet;
    HexahedralMesh hex;
    tet = poly;
    hex = poly;
    poly = tet;
    poly = hex;
    // TODO: add some content and ensure it transfers correctly.
}
