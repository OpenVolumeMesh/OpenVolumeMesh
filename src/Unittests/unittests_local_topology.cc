#include <gtest/gtest.h>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Unstable/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Unstable/Topology/TriangleTopology.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>


using namespace OpenVolumeMesh;

struct TetTopologyTest : public ::testing::Test
{
  TetrahedralMeshTopologyKernel mesh;
  //TopologyKernel mesh;
  VH a = mesh.add_vertex();
  VH b = mesh.add_vertex();
  VH c = mesh.add_vertex();
  VH d = mesh.add_vertex();
  CH ch = mesh.add_cell({a,b,c,d}, true);
};

TEST_F(TetTopologyTest, tetTopo_vertices)
{
  ASSERT_TRUE(ch.is_valid());
  auto tetTopo = TetTopology(mesh, ch, a);
  EXPECT_EQ(tetTopo.a(), a);
  // The following are an implementation detail, only orientation has
  // to be correcet
  EXPECT_EQ(tetTopo.b(), b);
  EXPECT_EQ(tetTopo.c(), c);
  EXPECT_EQ(tetTopo.d(), d);
}

TEST_F(TetTopologyTest, tetTopo_halfedges)
{
  auto tetTopo = TetTopology(mesh, ch, a);

  EXPECT_EQ(mesh.halfedge(tetTopo.ab()).from_vertex(), a);
  EXPECT_EQ(mesh.halfedge(tetTopo.ab()).to_vertex(),   b);
  EXPECT_EQ(mesh.halfedge(tetTopo.ba()).from_vertex(), b);
  EXPECT_EQ(mesh.halfedge(tetTopo.ba()).to_vertex(),   a);
  EXPECT_EQ(mesh.halfedge(tetTopo.ac()).from_vertex(), a);
  EXPECT_EQ(mesh.halfedge(tetTopo.ac()).to_vertex(),   c);
  EXPECT_EQ(mesh.halfedge(tetTopo.ca()).from_vertex(), c);
  EXPECT_EQ(mesh.halfedge(tetTopo.ca()).to_vertex(),   a);
  EXPECT_EQ(mesh.halfedge(tetTopo.ad()).from_vertex(), a);
  EXPECT_EQ(mesh.halfedge(tetTopo.ad()).to_vertex(),   d);
  EXPECT_EQ(mesh.halfedge(tetTopo.da()).from_vertex(), d);
  EXPECT_EQ(mesh.halfedge(tetTopo.da()).to_vertex(),   a);
  EXPECT_EQ(mesh.halfedge(tetTopo.bc()).from_vertex(), b);
  EXPECT_EQ(mesh.halfedge(tetTopo.bc()).to_vertex(),   c);
  EXPECT_EQ(mesh.halfedge(tetTopo.cb()).from_vertex(), c);
  EXPECT_EQ(mesh.halfedge(tetTopo.cb()).to_vertex(),   b);
  EXPECT_EQ(mesh.halfedge(tetTopo.bd()).from_vertex(), b);
  EXPECT_EQ(mesh.halfedge(tetTopo.bd()).to_vertex(),   d);
  EXPECT_EQ(mesh.halfedge(tetTopo.db()).from_vertex(), d);
  EXPECT_EQ(mesh.halfedge(tetTopo.db()).to_vertex(),   b);
  EXPECT_EQ(mesh.halfedge(tetTopo.cd()).from_vertex(), c);
  EXPECT_EQ(mesh.halfedge(tetTopo.cd()).to_vertex(),   d);
  EXPECT_EQ(mesh.halfedge(tetTopo.dc()).from_vertex(), d);
  EXPECT_EQ(mesh.halfedge(tetTopo.dc()).to_vertex(),   c);
}

TEST_F(TetTopologyTest, tetTopo_halffaces)
{
  auto tetTopo = TetTopology(mesh, ch, a);

  EXPECT_EQ(tetTopo.bdc(), mesh.find_halfface({b,d,c}));
  EXPECT_EQ(tetTopo.dcb(), mesh.find_halfface({d,c,b}));
  EXPECT_EQ(tetTopo.cbd(), mesh.find_halfface({c,b,d}));

  EXPECT_EQ(tetTopo.acd(), mesh.find_halfface({a,c,d}));
  EXPECT_EQ(tetTopo.cda(), mesh.find_halfface({c,d,a}));
  EXPECT_EQ(tetTopo.dac(), mesh.find_halfface({d,a,c}));

  EXPECT_EQ(tetTopo.bad(), mesh.find_halfface({b,a,d}));
  EXPECT_EQ(tetTopo.adb(), mesh.find_halfface({a,d,b}));
  EXPECT_EQ(tetTopo.dba(), mesh.find_halfface({d,b,a}));

  EXPECT_EQ(tetTopo.abc(), mesh.find_halfface({a,b,c}));
  EXPECT_EQ(tetTopo.bca(), mesh.find_halfface({b,c,a}));
  EXPECT_EQ(tetTopo.cab(), mesh.find_halfface({c,a,b}));

}

TEST_F(TetTopologyTest, tetTopo_constructors)
{
  HFH abc = mesh.find_halfface({a,b,c});

  auto tetTopo1 = TetTopology(mesh, ch, abc, a);
  auto tetTopo2 = TetTopology(mesh, ch, abc);
  auto tetTopo3 = TetTopology(mesh, abc, a);
  auto tetTopo4 = TetTopology(mesh, abc);
  auto tetTopo5 = TetTopology(mesh, ch, a);
  auto tetTopo6 = TetTopology(mesh, ch);

  EXPECT_EQ(tetTopo1, tetTopo2);
  EXPECT_EQ(tetTopo1, tetTopo3);
  EXPECT_EQ(tetTopo1, tetTopo4);
  EXPECT_EQ(tetTopo1, tetTopo5);
  EXPECT_EQ(tetTopo1, tetTopo6);
}

TEST_F(TetTopologyTest, tetTopo_triTopo)
{
  HFH abc = mesh.find_halfface({a,b,c});
  HFH bcd = mesh.find_halfface({b,c,d});
  auto tt = TetTopology(mesh, ch, abc, a);

  auto tri_abc = tt.triangle_topology<TetTopology::ABC>();
  EXPECT_EQ(tri_abc, TriangleTopology(mesh, abc, a));
  auto tri_bca = tt.triangle_topology<TetTopology::BCA>();
  EXPECT_EQ(tri_bca, TriangleTopology(mesh, abc, b));
  auto tri_cab = tt.triangle_topology<TetTopology::CAB>();
  EXPECT_EQ(tri_cab, TriangleTopology(mesh, abc, c));

  auto tri_acb = tt.triangle_topology<TetTopology::ACB>();
  EXPECT_EQ(tri_acb, TriangleTopology(mesh, abc.opposite_handle(), a));

  auto tri_dcb = tt.triangle_topology<TetTopology::DCB>();
  EXPECT_EQ(tri_dcb, TriangleTopology(mesh, bcd.opposite_handle(), d));

  auto tri_bcd = tt.triangle_topology<TetTopology::BCD>();
  EXPECT_EQ(tri_bcd, TriangleTopology(mesh, bcd, b));

  // TODO: add more tests!

}

TEST_F(TetTopologyTest, tetTopo_triTopo_dynamic)
{
  HFH abc = mesh.find_halfface({a,b,c});
  HFH bcd = mesh.find_halfface({b,c,d});
  auto tt = TetTopology(mesh, ch, abc, a);

  auto tri_abc = tt.triangle_topology(TetTopology::ABC);
  EXPECT_EQ(tri_abc, TriangleTopology(mesh, abc, a));
  auto tri_dcb = tt.triangle_topology(TetTopology::DCB);
  EXPECT_EQ(tri_dcb, TriangleTopology(mesh, bcd.opposite_handle(), d));
}
TEST_F(TetTopologyTest, tetTopo_getlabel)
{
  HFH abc = mesh.find_halfface({a,b,c});
  auto tt = TetTopology(mesh, ch, abc, a);
  EXPECT_EQ(tt.get_label(abc, a), TetTopology::ABC);
  EXPECT_EQ(tt.get_label(abc, b), TetTopology::BCA);
  EXPECT_EQ(tt.get_label(abc, c), TetTopology::CAB);
  EXPECT_EQ(tt.get_label(abc.opposite_handle(), a), TetTopology::ACB);
}
