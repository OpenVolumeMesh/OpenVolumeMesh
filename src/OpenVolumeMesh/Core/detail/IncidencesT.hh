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

    bool enabled() const {return incident_.has_value();}
    void setEnabled(bool enable);

    Incidences const& incident(Handle _h) const;

protected:
    const Derived *topo() const {return static_cast<const Derived*>(this);}
    Incidences & incident(Handle _h);

    bool valid(Handle vh) const;
    virtual void recompute() = 0;

private:
    std::optional<PrivateProperty<Incidences, Entity>> incident_;
};


} // namespace OpenVolumeMesh
