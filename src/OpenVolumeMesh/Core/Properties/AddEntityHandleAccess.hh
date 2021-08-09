#pragma once

#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

namespace OpenVolumeMesh {
template<typename EntityTag, typename Parent>
class AddEntityHandleAccess : public Parent
{
    static_assert(is_entity<EntityTag>::value);

    using EntityHandleT = HandleT<EntityTag>;

public:
    using Parent::Parent;

    using reference = typename Parent::reference;
    using const_reference = typename Parent::const_reference;

    /// No range check performed!
    reference operator[](const EntityHandleT& _h) { return Parent::operator[](_h.uidx()); }

    /// No range check performed!
    const_reference operator[](const EntityHandleT& _h) const { return Parent::operator[](_h.uidx()); }

    reference       at(const EntityHandleT& _h) { return Parent::at(_h.uidx()); }
    const_reference at(const EntityHandleT& _h) const { return Parent::at(_h.uidx()); }


};
} // namespace OpenVolumeMesh
