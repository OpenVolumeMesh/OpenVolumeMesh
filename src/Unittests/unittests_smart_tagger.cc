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

TEST_F(TetrahedralMeshBase, SmartTaggerCurrentBase) {

    generateTetrahedralMesh(mesh_);

    auto tagger = SmartTagger<Entity::Vertex, uint8_t, bool>(mesh_, 2);
    auto vh_true = VH(0);
    auto vh_false = VH(1);

    tagger.set(vh_true, true);
    tagger.set(vh_false, false);
    EXPECT_EQ(tagger.current_base(), 0);
    EXPECT_EQ(tagger.get(vh_true), true);
    EXPECT_EQ(tagger.get(vh_false), false);

    for (uint8_t i = 0; i <= std::numeric_limits<uint8_t>::max() - 2*2; i += 2) {
        tagger.reset();
        EXPECT_EQ(tagger.current_base(), i + 2);

        EXPECT_EQ(tagger.get(vh_true), false);
        EXPECT_EQ(tagger.get(vh_false), false);
        tagger.set(vh_true, true);
        tagger.set(vh_false, false);
        EXPECT_EQ(tagger.get(vh_true), true);
        EXPECT_EQ(tagger.get(vh_false), false);
    }

    tagger.reset();
    EXPECT_EQ(tagger.current_base(), 0);

    EXPECT_EQ(tagger.get(vh_true), false);
    EXPECT_EQ(tagger.get(vh_false), false);
    tagger.set(vh_true, true);
    tagger.set(vh_false, false);
    EXPECT_EQ(tagger.get(vh_true), true);
    EXPECT_EQ(tagger.get(vh_false), false);
}


// TODO: test SmartTagger overflow
