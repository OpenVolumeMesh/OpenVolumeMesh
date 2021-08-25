#include "unittests_common.hh"
#include <OpenVolumeMesh/Util/SmartTagger.hh>


using namespace OpenVolumeMesh;

TEST_F(TetrahedralMeshBase, SmartTagger) {

    generateTetrahedralMesh(mesh_);

    auto tagger = SmartTagger<Entity::Vertex>(mesh_);
    for (auto vh: mesh_.vertices()) {
      EXPECT_FALSE(tagger[vh]);
    }
    tagger.tag(VertexHandle(1));
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 1);
    }
    tagger.tag(VertexHandle(2));
    tagger.untag(VertexHandle(1));
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 2);
    }
    tagger.reset();
    for (auto vh: mesh_.vertices()) {
      EXPECT_FALSE(tagger[vh]);
    }
}

// TODO: test SmartTagger overflow
