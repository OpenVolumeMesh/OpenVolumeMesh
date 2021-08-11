#pragma once

#include <vector>
#include <optional>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/Properties/PrivateProperty.hh>

namespace OpenVolumeMesh {

class ResourceManager;

template<typename Derived, typename _Entity, typename _Incidences>
class IncidencesT
{
public:
    using Entity = _Entity;
    using Incidences = _Incidences;
    static_assert(is_entity<_Entity>::value);
    using Handle = HandleT<Entity>;

    IncidencesT() = default;
    bool enabled() const {return enabled_;}
    void set_enabled(bool enable);
    void deleted(Handle);

    Incidences const& incident(Handle _h) const;

protected:
    const Derived *topo() const {return static_cast<const Derived*>(this);}
    Incidences & incident(Handle _h);
    /// HACK: this is only for lazy ordering of HE-HF incidences!
    Incidences & incident_mutable(Handle _h) const;

    bool valid(Handle vh) const;
    virtual void recompute() = 0;
    void resize();

    void swap(Handle _h1, Handle _h2);

private:
    bool enabled_ = false;
    // HACK: mutable;
    // we need this for the lazy ordering of HE-HF incidences :(
    // I wish we could have mutable inheritance from IncidencesT instead.
    mutable std::vector<_Incidences> incident_;
};

} // namespace OpenVolumeMesh
