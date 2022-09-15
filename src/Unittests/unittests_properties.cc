
#include <iostream>

#include "unittests_common.hh"

#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>

using namespace OpenVolumeMesh;
using namespace Geometry;

TEST_F(PolyhedralMeshBase, PropertySmartPointerTestSimple) {

    generatePolyhedralMesh(mesh_);

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("SimpleTestProp");
    }
}

TEST_F(PolyhedralMeshBase, PropertySmartPointerTest1) {

    generatePolyhedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());
    const size_t nvp = mesh_.n_vertex_props();

    VertexPropertyT<double> v_prop_d = mesh_.request_vertex_property<double>("MyVPropDbl");
    EXPECT_EQ(nvp + 1, mesh_.n_vertex_props());

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("MyVProp");

        v_prop[VertexHandle(0)] = 1.4f;

        VertexPropertyT<float> v_prop2(v_prop);

        VertexPropertyT<float> v_prop3 = v_prop;

        EXPECT_EQ(12u, v_prop3.size());

        EXPECT_EQ(nvp + 2, mesh_.n_vertex_props());

        VertexPropertyT<float> v_prop_duplicate = mesh_.request_vertex_property<float>("MyVProp");

        EXPECT_EQ(nvp + 2, mesh_.n_vertex_props());

        EXPECT_FLOAT_EQ(1.4f, v_prop3[VertexHandle(0)]);

        VertexPropertyT<std::string> v_prop_duplicate_2 = mesh_.request_vertex_property<std::string>("MyVProp");

        EXPECT_EQ(nvp + 3, mesh_.n_vertex_props());
    }

    mesh_.set_persistent(v_prop_d);

    EXPECT_EQ(nvp + 1, mesh_.n_vertex_props());

    mesh_.add_vertex(Vec3d(0.0, 0.0, 0.0));

    EXPECT_EQ(13u, mesh_.n_vertices());

    VertexPropertyT<double> v_prop_d2 = mesh_.request_vertex_property<double>("MyVPropDbl");

    EXPECT_EQ(nvp + 1, mesh_.n_vertex_props());

    HalfEdgePropertyT<int> he_prop = mesh_.request_halfedge_property<int>("MyHEProp");

    EXPECT_EQ(40u, he_prop.size());

    mesh_.add_edge(VertexHandle(0), VertexHandle(2));

    EXPECT_EQ(42u, he_prop.size());
}

TEST_F(HexahedralMeshBase, PropertySmartPointerPersistencyTest1) {

    generateHexahedralMesh(mesh_);

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("FloatVProp");

        v_prop[VertexHandle(0)] = 24.5f;
        v_prop[VertexHandle(11)] = 2.34f;

        mesh_.set_persistent(v_prop);
    }

    VertexPropertyT<float> v_prop2 = mesh_.request_vertex_property<float>("FloatVProp");

    EXPECT_FLOAT_EQ(24.5f, v_prop2[VertexHandle(0)]);
    EXPECT_FLOAT_EQ(2.34f, v_prop2[VertexHandle(11)]);
}

TEST_F(HexahedralMeshBase, PropertySmartPointerPersistencyTest2) {

    generateHexahedralMesh(mesh_);
    const size_t nvp = mesh_.n_vertex_props();

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("FloatVProp");

        mesh_.set_persistent(v_prop);
    }

    EXPECT_EQ(nvp + 1, mesh_.n_vertex_props());

    for(int i = 0; i < 1; ++i) {

        VertexPropertyT<float> v_prop = mesh_.request_vertex_property<float>("FloatVProp");

        EXPECT_EQ(nvp + 1, mesh_.n_vertex_props());

        mesh_.set_persistent(v_prop, false);
    }

    EXPECT_EQ(nvp, mesh_.n_vertex_props());
}

TEST_F(HexahedralMeshBase, AnonymousPropertiesTest1) {

    generateHexahedralMesh(mesh_);

    CellPropertyT<float> c_prop = mesh_.request_cell_property<float>();
    CellPropertyT<float> c_prop2 = c_prop;

    EXPECT_EQ(1u, mesh_.n_cell_props());

    // only shared props can be persistent
    EXPECT_THROW(mesh_.set_persistent(c_prop), std::runtime_error);

    // anonymous properties cannot be be shared
    EXPECT_THROW(mesh_.set_shared(c_prop), std::runtime_error);

    c_prop.set_name("test");

    EXPECT_EQ(c_prop2.name(), "test");


    // only shared props can be persistent
    EXPECT_THROW(mesh_.set_persistent(c_prop), std::runtime_error);
    mesh_.set_shared(c_prop);
    mesh_.set_persistent(c_prop);

    EXPECT_EQ(1u, mesh_.n_cell_props());
    EXPECT_EQ(1u, mesh_.n_persistent_props<Entity::Cell>());

    mesh_.set_persistent(c_prop2, false);

    EXPECT_EQ(1u, mesh_.n_cell_props());
    EXPECT_EQ(0u, mesh_.n_persistent_props<Entity::Cell>());
}

TEST_F(HexahedralMeshBase, AnonymousPropertiesTest2) {

    generateHexahedralMesh(mesh_);

    CellPropertyT<float> c_prop1 = mesh_.request_cell_property<float>();

    for(int i = 0; i < 1; ++i) {

        CellPropertyT<float> c_prop2 = mesh_.request_cell_property<float>();
        EXPECT_EQ(2u, mesh_.n_cell_props());
    }

    EXPECT_EQ(1u, mesh_.n_cell_props());
}

TEST_F(PolyhedralMeshBase, StatusTest) {

    generatePolyhedralMesh(mesh_);

    EXPECT_EQ(12u, mesh_.n_vertices());

    StatusAttrib status(mesh_);
}

TEST_F(PolyhedralMeshBase, PropertyOnConstMesh)
{
    auto test = [](auto const &mesh) {
        auto prop = mesh.template create_private_property<int, Entity::Vertex>();
    };


    generatePolyhedralMesh(mesh_);
    test(mesh_);
}


TEST_F(PolyhedralMeshBase, PropertyAtAccess)
{
    generatePolyhedralMesh(mesh_);
    auto prop = mesh_.create_private_property<int, Entity::Vertex>();
    prop.at(VertexHandle(1)) = 12345;
    EXPECT_EQ(prop.at(VertexHandle(1)), 12345);
    try {
        prop.at(VertexHandle(654321)) = 1;
        FAIL() << "Expected exception";
    } catch (std::out_of_range &) {
    }

}

TEST_F(PolyhedralMeshBase, EntitySwapInt)
{
    generatePolyhedralMesh(mesh_);
    auto prop = mesh_.create_private_property<int, Entity::Vertex>();
    auto v0 = VertexHandle(0);
    auto v1 = VertexHandle(1);
    prop[v0] = 23;
    prop[v1] = 42;
    mesh_.swap_vertex_indices(v0, v1);
    EXPECT_EQ(prop[v0], 42);
    EXPECT_EQ(prop[v1], 23);
}

TEST_F(PolyhedralMeshBase, EntitySwapBool)
{
    generatePolyhedralMesh(mesh_);
    auto prop = mesh_.create_private_property<bool, Entity::Vertex>();
    auto v0 = VertexHandle(0);
    auto v1 = VertexHandle(1);
    prop[v0] = false;
    prop[v1] = true;
    mesh_.swap_vertex_indices(v0, v1);
    EXPECT_EQ(prop[v0], true);
    EXPECT_EQ(prop[v1], false);
}



