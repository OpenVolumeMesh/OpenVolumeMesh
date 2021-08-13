#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

namespace OpenVolumeMesh {

/// Given a "Parent" class, likely an std::vector,
/// replace operator[](size_t) and at(size_t) with handle-type-safe versions

template<typename EntityTag, typename Parent>
class HandleIndexing : public Parent
{
    static_assert(is_entity<EntityTag>::value);

    using EntityHandleT = HandleT<EntityTag>;

public:
    using Parent::Parent;

    using reference = typename Parent::reference;
    using const_reference = typename Parent::const_reference;

    reference       operator[](size_t) = delete;
    reference       operator[](const EntityHandleT _h)       { return Parent::operator[](_h.uidx()); }
    const_reference operator[](size_t) const = delete;
    const_reference operator[](const EntityHandleT _h) const { return Parent::operator[](_h.uidx()); }

    reference       at(const EntityHandleT _h)       { return Parent::at(_h.uidx()); }
    const_reference at(const EntityHandleT _h) const { return Parent::at(_h.uidx()); }
};

} // namespace OpenVolumeMesh
