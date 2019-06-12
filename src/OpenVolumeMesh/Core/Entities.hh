#pragma once

#include <type_traits>

namespace OpenVolumeMesh {

namespace Entity {
    struct Vertex   {};
    struct Edge     {};
    struct HalfEdge {};
    struct Face     {};
    struct HalfFace {};
    struct Cell     {};
    struct Mesh     {};
}

template<typename T>
struct is_entity : std::false_type {};

template<> struct is_entity<Entity::Vertex>   : std::true_type {};
template<> struct is_entity<Entity::Edge>     : std::true_type {};
template<> struct is_entity<Entity::HalfEdge> : std::true_type {};
template<> struct is_entity<Entity::Face>     : std::true_type {};
template<> struct is_entity<Entity::HalfFace> : std::true_type {};
template<> struct is_entity<Entity::Cell>     : std::true_type {};
template<> struct is_entity<Entity::Mesh>     : std::true_type {};

} // namespace OpenVolumeMesh
