#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/Properties/AddEntityHandleAccess.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>

namespace OpenVolumeMesh {
template<typename T, typename EntityTag>
// TODO: need better name, the special thing about this is
//       the lack of shared_ptr indirection overhead
class PrivateProperty : public AddEntityHandleAccess<EntityTag, PropertyStorageT<T>>
{
    static_assert(is_entity<EntityTag>::value);
public:
    PrivateProperty(const ResourceManager *resman, T const& _def = T())
        : AddEntityHandleAccess<EntityTag, PropertyStorageT<T>>(
        &resman->storage_tracker<EntityTag>(),
              "", EntityTag::type(), _def)
    {}

};
} // namespace OpenVolumeMesh
