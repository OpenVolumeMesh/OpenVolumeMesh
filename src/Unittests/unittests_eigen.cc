#include <gtest/gtest.h>

#include <Eigen/Core>
#include <OpenVolumeMesh/IO/PropertyCodecsEigen.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Geometry/EigenTraits.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>

using namespace OpenVolumeMesh;
using namespace Geometry;

class EigenTest : public testing::Test
{
};

TEST_F(EigenTest, LoadSaveProps)
{
  using Mesh = GeometryKernel<Eigen::Vector3d, TetrahedralMeshTopologyKernel>;
  Mesh mesh;
  auto vh0 = mesh.add_vertex(Eigen::Vector3d{1, 2, 3});
  auto vh1 = mesh.add_vertex(Eigen::Vector3d{-7, 99, 87});
  auto vh2 = mesh.add_vertex(Eigen::Vector3d{23, 42, 67});


  auto prop_v3d = mesh.request_vertex_property<Eigen::Vector2f>("vec2f");

  EXPECT_EQ(0u, mesh.n_persistent_props<Entity::Vertex>());
  prop_v3d[vh0] = {23,42};
  mesh.set_persistent(prop_v3d);
  EXPECT_EQ(1u, mesh.n_persistent_props<Entity::Vertex>());

  OpenVolumeMesh::IO::PropertyCodecs codecs;
  OpenVolumeMesh::IO::register_eigen_codecs(codecs);

  // Write file
  auto write_opts = OpenVolumeMesh::IO::WriteOptions{};
  auto write_result = OpenVolumeMesh::IO::ovmb_write("eigen.ovmb", mesh, write_opts, codecs);
  ASSERT_EQ(write_result, OpenVolumeMesh::IO::WriteResult::Ok);

  Mesh loaded;
  auto read_opts = OpenVolumeMesh::IO::ReadOptions{};
  auto read_result = OpenVolumeMesh::IO::ovmb_read("eigen.ovmb", loaded, read_opts, codecs);
  ASSERT_EQ(read_result, OpenVolumeMesh::IO::ReadResult::Ok);

  EXPECT_EQ(1u, loaded.n_persistent_props<Entity::Vertex>());

  auto maybe_prop_v2d= loaded.get_vertex_property<Eigen::Vector2f>("vec2f");
  ASSERT_TRUE(maybe_prop_v2d.has_value());

  ASSERT_EQ((*maybe_prop_v2d)[vh0], (Eigen::Vector2f{23,42}));
  ASSERT_EQ(loaded.vertex(vh0), (Eigen::Vector3d{1, 2, 3}));
  ASSERT_EQ(loaded.vertex(vh1), (Eigen::Vector3d{-7, 99, 87}));
  ASSERT_EQ(loaded.vertex(vh2), (Eigen::Vector3d{23, 42, 67}));
}

