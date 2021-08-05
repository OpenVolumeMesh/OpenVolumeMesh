#pragma once

#include "OpenVolumeMesh/Config/Export.hh"
#include <type_traits>

namespace OpenVolumeMesh {

enum class EntityType {
    Vertex, Edge, HalfEdge, Face, HalfFace, Cell, Mesh
};

/// EntityTag:
namespace Entity {
struct OVM_EXPORT Vertex   { Vertex()   = delete; static EntityType type() {return EntityType::Vertex;}};
struct OVM_EXPORT Edge     { Edge()     = delete; static EntityType type() {return EntityType::Edge;}};
struct OVM_EXPORT HalfEdge { HalfEdge() = delete; static EntityType type() {return EntityType::HalfEdge;}};
struct OVM_EXPORT Face     { Face()     = delete; static EntityType type() {return EntityType::Face;}};
struct OVM_EXPORT HalfFace { HalfFace() = delete; static EntityType type() {return EntityType::HalfFace;}};
struct OVM_EXPORT Cell     { Cell()     = delete; static EntityType type() {return EntityType::Cell;}};
struct OVM_EXPORT Mesh     { Mesh()     = delete; static EntityType type() {return EntityType::Mesh;}};
} // namespace Entity


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
