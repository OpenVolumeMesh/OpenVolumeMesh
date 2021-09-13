#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <type_traits>
#include <cstddef>
#include <cstdint>

namespace OpenVolumeMesh {

enum class EntityType : uint8_t {
    Vertex, Edge, HalfEdge, Face, HalfFace, Cell, Mesh
};
inline static const size_t n_entity_types = 7;

/// EntityTag:
namespace Entity {
struct OVM_EXPORT Vertex   { static EntityType type() {return EntityType::Vertex;}};
struct OVM_EXPORT Edge     { static EntityType type() {return EntityType::Edge;}};
struct OVM_EXPORT HalfEdge { static EntityType type() {return EntityType::HalfEdge;}};
struct OVM_EXPORT Face     { static EntityType type() {return EntityType::Face;}};
struct OVM_EXPORT HalfFace { static EntityType type() {return EntityType::HalfFace;}};
struct OVM_EXPORT Cell     { static EntityType type() {return EntityType::Cell;}};
struct OVM_EXPORT Mesh     { static EntityType type() {return EntityType::Mesh;}};
} // namespace Entity


template<typename T>
struct is_entity : std::false_type {};

template<typename T>
inline constexpr bool is_entity_v = is_entity<T>::value;


template<> struct OVM_EXPORT is_entity<Entity::Vertex>   : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Edge>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfEdge> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Face>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfFace> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Cell>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Mesh>     : std::true_type {};

} // namespace OpenVolumeMesh
