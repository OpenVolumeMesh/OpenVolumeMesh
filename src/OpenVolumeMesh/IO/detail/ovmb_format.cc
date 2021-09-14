#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <string>

namespace OpenVolumeMesh::IO::detail {

const std::array<uint8_t, 8> ovmb_magic {'O', 'V', 'M', 'B', '\n', '\r', '\n', 0xff};

EntityType as_entity_type(PropertyEntity pe)
{
    switch (pe) {
    case PropertyEntity::Vertex:   return EntityType::Vertex;
    case PropertyEntity::Edge:     return EntityType::Edge;
    case PropertyEntity::Face:     return EntityType::Face;
    case PropertyEntity::Cell:     return EntityType::Cell;
    case PropertyEntity::HalfEdge: return EntityType::HalfEdge;
    case PropertyEntity::HalfFace: return EntityType::HalfFace;
    case PropertyEntity::Mesh:     return EntityType::Mesh;
    }
    throw std::runtime_error(std::string("unknown property entity: ") +
                             std::to_string(static_cast<size_t>(pe)));
}

PropertyEntity as_prop_entity(EntityType et)
{
    switch (et) {
    case EntityType::Vertex:   return PropertyEntity::Vertex;
    case EntityType::Edge:     return PropertyEntity::Edge;
    case EntityType::Face:     return PropertyEntity::Face;
    case EntityType::Cell:     return PropertyEntity::Cell;
    case EntityType::HalfEdge: return PropertyEntity::HalfEdge;
    case EntityType::HalfFace: return PropertyEntity::HalfFace;
    case EntityType::Mesh:     return PropertyEntity::Mesh;
    }
    throw std::runtime_error(std::string("unknown property entity: ") +
                             std::to_string(static_cast<size_t>(et)));
}

IntEncoding suitable_int_encoding(uint32_t max_value)
{
   if (max_value <= std::numeric_limits<uint8_t>::max())  return IntEncoding::U8;
   if (max_value <= std::numeric_limits<uint16_t>::max())  return IntEncoding::U16;
   return IntEncoding::U32;
}

} //namespace OpenVolumeMesh::IO::detail
