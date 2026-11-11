
#include "unittests_common.hh"

#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>

using namespace OpenVolumeMesh;
using namespace OpenVolumeMesh::Geometry;

// Test fixture specifically for TetrahedralGeometryKernel
class TetrahedralGeometryKernelBase : public ::testing::Test
{
protected:

    using Mesh = TetrahedralGeometryKernel<
        Vec3d,
        TetrahedralMeshTopologyKernel
    >;

    using VH = VertexHandle;
    using CH = CellHandle;

    Mesh mesh_;

    void SetUp() override
    {
        mesh_.enable_deferred_deletion(false);
        mesh_.enable_fast_deletion(false);
    }
};

// Test signed volume
TEST_F(TetrahedralGeometryKernelBase, SignedVolume)
{
    VH v0 = mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));
    VH v1 = mesh_.add_vertex(Vec3d(0.0, 1.0, 0.0));
    VH v2 = mesh_.add_vertex(Vec3d(0.0, 0.0, 1.0));
    VH v3 = mesh_.add_vertex(Vec3d(1.0, 0.0, 0.0));

    CH ch = mesh_.add_cell(v0, v1, v2, v3);

    EXPECT_DOUBLE_EQ(mesh_.signed_volume(ch), 1.0 / 6.0);

    // invert cell, should have negative volume after
    mesh_.set_vertex(v3,{-1.0, 0.0, 0.0});
    EXPECT_DOUBLE_EQ(mesh_.signed_volume(ch), -1.0 / 6.0);
}
