#include "unittests_common.hh"
#include <OpenVolumeMesh/Util/SmartTagger.hh>


using namespace OpenVolumeMesh;

TEST_F(TetrahedralMeshBase, SmartTaggerBool) {

    generateTetrahedralMesh(mesh_);

    auto tagger = SmartTaggerBool<Entity::Vertex>(mesh_);
    for (auto vh: mesh_.vertices()) {
      EXPECT_FALSE(tagger[vh]);
    }
    tagger.set(VertexHandle(1), true);
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 1);
    }
    tagger.set(VertexHandle(2), true);
    tagger.set(VertexHandle(1), false);
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 2);
    }
    tagger.reset();
    for (auto vh: mesh_.vertices()) {
      EXPECT_FALSE(tagger[vh]);
    }
}
TEST_F(TetrahedralMeshBase, SmartTaggerEnum) {

  enum class TestEnum { Zero, One, Two };

    generateTetrahedralMesh(mesh_);

    auto tagger = SmartTagger<Entity::Vertex, uint8_t, TestEnum>(mesh_, 3);
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], TestEnum::Zero);
    }
    tagger.set(VertexHandle(1), TestEnum::One);
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 1 ? TestEnum::One : TestEnum::Zero);
    }
    tagger.set(VertexHandle(2), TestEnum::Two);
    tagger.set(VertexHandle(1), TestEnum::Zero);
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], vh.idx() == 2 ? TestEnum::Two : TestEnum::Zero);
    }
    tagger.reset();
    for (auto vh: mesh_.vertices()) {
      EXPECT_EQ(tagger[vh], TestEnum::Zero);
    }
    for (int i = 0; i < 100; ++i) {
      tagger.reset();
      for (auto vh: mesh_.vertices()) {
        EXPECT_EQ(tagger[vh], TestEnum::Zero);
      }
    }
}


// TODO: test SmartTagger overflow
