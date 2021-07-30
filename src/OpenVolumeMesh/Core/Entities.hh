#pragma once

#include "OpenVolumeMesh/Config/Export.hh"
#include <type_traits>

namespace OpenVolumeMesh {

/// Entity tags
namespace Entity {
    struct OVM_EXPORT Vertex   { Vertex()   = delete;};
    struct OVM_EXPORT Edge     { Edge()     = delete;};
    struct OVM_EXPORT HalfEdge { HalfEdge() = delete;};
    struct OVM_EXPORT Face     { Face()     = delete;};
    struct OVM_EXPORT HalfFace { HalfFace() = delete;};
    struct OVM_EXPORT Cell     { Cell()     = delete;};
    struct OVM_EXPORT Mesh     { Mesh()     = delete;};
}

template<typename T>
struct is_entity : std::false_type {};

template<> struct OVM_EXPORT is_entity<Entity::Vertex>   : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Edge>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfEdge> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Face>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfFace> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Cell>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Mesh>     : std::true_type {};

} // namespace OpenVolumeMesh
