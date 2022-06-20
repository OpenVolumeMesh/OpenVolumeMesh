#include "unittests_common.hh"

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>
#include <OpenVolumeMesh/Attribs/InterfaceAttrib.hh>
#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include <OpenVolumeMesh/Attribs/TexCoordAttrib.hh>


using namespace OpenVolumeMesh;
using OpenVolumeMesh::Geometry::Vec3d;

class AttribTest : public testing::Test
{
public:
  virtual void SetUp()
  {
      v0_ = mesh_.add_vertex();
      v1_ = mesh_.add_vertex();
  }
  TopologyKernel mesh_;
  VertexHandle v0_, v1_;
};

TEST_F(AttribTest, StatusAttrib)
{
  auto status = StatusAttrib(mesh_);
  EXPECT_EQ(status[v0_].selected(), false);
  EXPECT_EQ(status[v0_].tagged(), false);
  EXPECT_EQ(status[v0_].deleted(), false);
  EXPECT_EQ(status[v0_].hidden(), false);
  status[v0_].set_selected(true);
  EXPECT_EQ(status[v0_].selected(), true);
  status[v0_].set_tagged(true);
  EXPECT_EQ(status[v0_].tagged(), true);
  status[v0_].set_hidden(true);
  EXPECT_EQ(status[v0_].hidden(), true);
}

TEST_F(AttribTest, InterfaceAttrib)
{
  auto iface = InterfaceAttrib(mesh_);
  EXPECT_EQ(iface[v0_], false);
  iface[v0_] = true;
  EXPECT_EQ(iface[v0_], true);
  // TODO: test other entitiesr (edges, faces)
}

// TOOD: write tests for remaining attribs
// TOOD: test assignment (and ensure consistent behaviour)
