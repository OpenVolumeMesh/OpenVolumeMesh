#include "unittests_common.hh"
#include <gtest/gtest-typed-test.h>

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>

using namespace OpenVolumeMesh;

using namespace Geometry;

template<bool Binary>
class LoadSave {
public:
    template<typename MeshT>
    bool readFile(const char *filename,
              MeshT &mesh,
              bool _topo_check = true,
              bool _bottom_up = true) const
    {
        if constexpr(Binary) {
            std::string ovmb_filename = filename;
            ovmb_filename += 'b';
            OpenVolumeMesh::IO::ReadOptions options;
            options.topology_check = _topo_check;
            options.bottom_up_incidences = _bottom_up;
            auto result = OpenVolumeMesh::IO::ovmb_read(ovmb_filename.c_str(), mesh, options);
            return result == OpenVolumeMesh::IO::ReadResult::Ok;
        } else {
            OpenVolumeMesh::IO::FileManager fileManager;
            return fileManager.readFile(filename, mesh, _topo_check, _bottom_up);
        }
    }
    template<typename MeshT>
    bool writeFile(const char *filename,
              MeshT const &mesh) const
    {
        if constexpr(Binary) {
            std::string ovmb_filename = filename;
            ovmb_filename += 'b';
            auto result = OpenVolumeMesh::IO::ovmb_write(ovmb_filename.c_str(), mesh);
            return result == OpenVolumeMesh::IO::WriteResult::Ok;
        } else {
            OpenVolumeMesh::IO::FileManager fileManager;
            return fileManager.writeFile(filename, mesh);
        }
    }
    bool isHexahedralMesh(const char *filename) const
    {
        if constexpr(Binary) {
            std::string ovmb_filename = filename;
            ovmb_filename += 'b';
            std::ifstream stream(ovmb_filename, std::ios::binary);
            auto result = OpenVolumeMesh::IO::make_ovmb_reader(
                        stream,
                        OpenVolumeMesh::IO::ReadOptions(),
                        OpenVolumeMesh::IO::g_default_property_codecs);
            auto topo_type = result->topo_type();
            return topo_type.has_value() &&
                    topo_type.value() == OpenVolumeMesh::IO::TopoType::Hexahedral;
        } else {
            OpenVolumeMesh::IO::FileManager fileManager;
            return fileManager.isHexahedralMesh(filename);
        }
    }

};

template<typename use_binary>
class PolyhedralFileTest : public PolyhedralMeshBase
                         , public LoadSave<use_binary::value>
{
public:
    using LoadSave<use_binary::value>::readFile;
    using LoadSave<use_binary::value>::writeFile;
    using LoadSave<use_binary::value>::isHexahedralMesh;
};

template<typename use_binary>
class HexahedralFileTest : public HexahedralMeshBase
                         , public LoadSave<use_binary::value>
{
public:
    using LoadSave<use_binary::value>::readFile;
    using LoadSave<use_binary::value>::writeFile;
    using LoadSave<use_binary::value>::isHexahedralMesh;
};

using TrueFalse = testing::Types<std::true_type, std::false_type>;
TYPED_TEST_SUITE(PolyhedralFileTest, TrueFalse);
TYPED_TEST_SUITE(HexahedralFileTest, TrueFalse);

TYPED_TEST(PolyhedralFileTest, LoadFile)
{
  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));


  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());
  auto v0 = mesh.vertex(VH(0));
  EXPECT_DOUBLE_EQ(v0[0], 4.97203);
  EXPECT_DOUBLE_EQ(v0[1], 0.504645);
  EXPECT_DOUBLE_EQ(v0[2], 0.913112);
}

TYPED_TEST(PolyhedralFileTest, LoadNonManifoldMesh) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("NonManifold.ovm", mesh));

  EXPECT_EQ(12u, mesh.n_vertices());
  EXPECT_EQ(20u, mesh.n_edges());
  EXPECT_EQ(11u, mesh.n_faces());
  EXPECT_EQ(2u, mesh.n_cells());
}

TYPED_TEST(HexahedralFileTest, LoadFile) {

  HexahedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->isHexahedralMesh("Cylinder.ovm"));
  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());
}

TYPED_TEST(PolyhedralFileTest, SaveFile) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  // Write file
  ASSERT_TRUE(this->writeFile("Cylinder.copy.ovm", mesh));

  mesh.clear();

  ASSERT_TRUE(this->readFile("Cylinder.copy.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());
  auto v0 = mesh.vertex(VH(0));
  EXPECT_DOUBLE_EQ(v0[0], 4.97203);
  EXPECT_DOUBLE_EQ(v0[1], 0.504645);
  EXPECT_DOUBLE_EQ(v0[2], 0.913112);
}

TYPED_TEST(PolyhedralFileTest, SaveFileWithProps) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  // Attach non-persistent properties
  HalfFacePropertyT<float> hfprop = mesh.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop = mesh.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      hfprop[HalfFaceHandle(i)] = (float)i/2.0f;
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      vprop[VertexHandle(i)] = i;
  }

  mesh.set_persistent(hfprop);
  mesh.set_persistent(vprop);

  // Write file
  ASSERT_TRUE(this->writeFile("Cylinder.copy.ovm", mesh));

  mesh.clear();

  EXPECT_EQ(0u, mesh.n_persistent_props<Entity::HalfFace>());
  EXPECT_EQ(0u, mesh.n_persistent_props<Entity::Vertex>());

  ASSERT_TRUE(this->readFile("Cylinder.copy.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  EXPECT_EQ(1u, mesh.n_persistent_props<Entity::HalfFace>());
  EXPECT_EQ(1u, mesh.n_persistent_props<Entity::Vertex>());

  HalfFacePropertyT<float> hfprop2 = mesh.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop2 = mesh.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      EXPECT_FLOAT_EQ((float)i/2.0f, hfprop2[HalfFaceHandle(i)]);
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      EXPECT_EQ(i, vprop2[VertexHandle(i)]);
  }
}

TYPED_TEST(PolyhedralFileTest, SaveFileWithVectorProps) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  // Attach non-persistent properties
  HalfFacePropertyT<Vec3d> hfprop = mesh.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop = mesh.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      hfprop[HalfFaceHandle(i)] = Vec3d((double)i/2.0, (double)i/2.0, (double)i/2.0);
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      vprop[VertexHandle(i)] = Vec2i(i, i);
  }

  mesh.set_persistent(hfprop);
  mesh.set_persistent(vprop);

  // Write file
  ASSERT_TRUE(this->writeFile("Cylinder.copy.ovm", mesh));

  mesh.clear();
  EXPECT_EQ(0u, mesh.n_persistent_props<Entity::HalfFace>());
  EXPECT_EQ(0u, mesh.n_persistent_props<Entity::Vertex>());

  ASSERT_TRUE(this->readFile("Cylinder.copy.ovm", mesh, true, false));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  EXPECT_EQ(1u, mesh.n_persistent_props<Entity::HalfFace>());
  EXPECT_EQ(1u, mesh.n_persistent_props<Entity::Vertex>());

  HalfFacePropertyT<Vec3d> hfprop2 = mesh.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop2 = mesh.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][0]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][1]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][2]);
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      VertexHandle vh(i);
      EXPECT_EQ((int)i, vprop2[vh][0]);
      EXPECT_EQ((int)i, vprop2[vh][1]);
  }
}

TYPED_TEST(PolyhedralFileTest, SerializeVectorValuedProperties) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cylinder.ovm", mesh));

  EXPECT_EQ(399u, mesh.n_vertices());
  EXPECT_EQ(1070u, mesh.n_edges());
  EXPECT_EQ(960u, mesh.n_faces());
  EXPECT_EQ(288u, mesh.n_cells());

  // Attach persistent properties
  HalfFacePropertyT<Vec3d> hfprop = mesh.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop = mesh.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      hfprop[hfh] = Vec3d((double)i/2.0, (double)i/2.0, (double)i/2.0);
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      VertexHandle vh(i);
      vprop[vh] = Vec2i(i, i);
  }

  mesh.set_persistent(hfprop);
  mesh.set_persistent(vprop);

  std::ofstream ofs1("hfVecPropTest");
  std::ofstream ofs2("vVecPropTest");

  hfprop.serialize(ofs1);
  vprop.serialize(ofs2);

  ofs1.close();
  ofs2.close();

  /*
   * Change property values
   */
  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
	  hfprop[HalfFaceHandle(i)] = Vec3d((double)i/3.0, (double)i/3.0, (double)i/3.0);
  }

  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
	  vprop[VertexHandle(i)] = Vec2i(2*i, 2*i);
  }

  std::ifstream ifs1("hfVecPropTest");
  std::ifstream ifs2("vVecPropTest");

  hfprop.deserialize(ifs1);
  vprop.deserialize(ifs2);

  ifs1.close();
  ifs2.close();

  for(unsigned int i = 0; i < mesh.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][0]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][1]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][2]);
  }
  for(unsigned int i = 0; i < mesh.n_vertices(); ++i) {
      VertexHandle vh(i);
      EXPECT_EQ((int)i, vprop[vh][0]);
      EXPECT_EQ((int)i, vprop[vh][1]);
  }
}

TYPED_TEST(PolyhedralFileTest, LoadFileWithProps) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cube_with_props.ovm", mesh, true, false));

  EXPECT_EQ(8u, mesh.n_vertices());
  EXPECT_EQ(12u, mesh.n_edges());
  EXPECT_EQ(6u, mesh.n_faces());
  EXPECT_EQ(1u, mesh.n_cells());

  EXPECT_EQ(2u, mesh.n_props<Entity::Vertex>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Edge>());
  EXPECT_EQ(0u, mesh.n_props<Entity::HalfEdge>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Face>());
  EXPECT_EQ(1u, mesh.n_props<Entity::HalfFace>());
  EXPECT_EQ(0u, mesh.n_props<Entity::Cell>());
}

TYPED_TEST(PolyhedralFileTest, SaveFileWithProps2) {

  PolyhedralMesh &mesh = this->mesh_;

  ASSERT_TRUE(this->readFile("Cube_with_props.ovm", mesh, true, false));

  EXPECT_EQ(8u, mesh.n_vertices());
  EXPECT_EQ(12u, mesh.n_edges());
  EXPECT_EQ(6u, mesh.n_faces());
  EXPECT_EQ(1u, mesh.n_cells());

  EXPECT_EQ(2u, mesh.n_props<Entity::Vertex>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Edge>());
  EXPECT_EQ(0u, mesh.n_props<Entity::HalfEdge>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Face>());
  EXPECT_EQ(1u, mesh.n_props<Entity::HalfFace>());
  EXPECT_EQ(0u, mesh.n_props<Entity::Cell>());

  ASSERT_TRUE(this->writeFile("Cube_with_props.copy.ovm", mesh));

  mesh.clear();

  ASSERT_TRUE(this->readFile("Cube_with_props.copy.ovm", mesh, true, false));

  EXPECT_EQ(8u, mesh.n_vertices());
  EXPECT_EQ(12u, mesh.n_edges());
  EXPECT_EQ(6u, mesh.n_faces());
  EXPECT_EQ(1u, mesh.n_cells());

  EXPECT_EQ(2u, mesh.n_props<Entity::Vertex>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Edge>());
  EXPECT_EQ(0u, mesh.n_props<Entity::HalfEdge>());
  EXPECT_EQ(1u, mesh.n_props<Entity::Face>());
  EXPECT_EQ(1u, mesh.n_props<Entity::HalfFace>());
  EXPECT_EQ(0u, mesh.n_props<Entity::Cell>());
}

