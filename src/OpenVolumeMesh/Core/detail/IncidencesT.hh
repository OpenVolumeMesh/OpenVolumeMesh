#pragma once

#include <vector>
#include <optional>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/Properties/PrivateProperty.hh>

namespace OpenVolumeMesh {

class ResourceManager;

template<typename Derived, typename Entity, typename _Incidences>
class IncidencesT
{
    static_assert(is_entity<Entity>::value);
public:
    using Handle = HandleT<Entity>;
    using Incidences = _Incidences;

    IncidencesT() = default;
    // TODO IMPORTANT: implement copy/move constructor/= for this and for ve-hf::ordered_
    IncidencesT(IncidencesT<Derived, Entity, _Incidences> const &other);
    IncidencesT(IncidencesT<Derived, Entity, _Incidences> &&other);
    IncidencesT &operator=(IncidencesT<Derived, Entity, _Incidences> const &other);
    IncidencesT &operator=(IncidencesT<Derived, Entity, _Incidences> &&other);

    bool enabled() const {return incident_.has_value();}
    void set_enabled(bool enable);

    Incidences const& incident(Handle _h) const;

protected:
    const Derived *topo() const {return static_cast<const Derived*>(this);}
    Incidences & incident(Handle _h);
    /// HACK: this is only for lazy ordering of HE-HF incidences!
    Incidences & incident_mutable(Handle _h) const;

    bool valid(Handle vh) const;
    virtual void recompute() = 0;

private:
    // HACK: mutable;
    // we need this for the lazy ordering of HE-HF incidences :(
    // I wish we could have mutable inheritance from IncidencesT instead.
    mutable std::optional<PrivateProperty<Incidences, Entity>> incident_;
};

} // namespace OpenVolumeMesh
