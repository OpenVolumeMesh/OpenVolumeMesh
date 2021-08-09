#pragma once
#include <OpenVolumeMesh/Core/TopologyKernel.hh>


namespace OpenVolumeMesh {

template<typename Derived, typename Entity, typename _Incidences>
_Incidences &
IncidencesT<Derived, Entity, _Incidences>::
incident(Handle _h)
{
    assert(enabled());
    assert(_h.uidx() < incident_->size());
    return (*incident_)[_h];
}

template<typename Derived, typename Entity, typename _Incidences>
_Incidences const &
IncidencesT<Derived, Entity, _Incidences>::
incident(Handle _h) const
{
    assert(enabled());
    assert(_h.uidx() < incident_->size());
    return (*incident_)[_h];
}

template<typename Derived, typename Entity, typename _Incidences>
void
IncidencesT<Derived, Entity, _Incidences>::
setEnabled(bool enable)
{
    if (enabled() == enable)
        return;
    if (enable) {
        incident_ = PrivateProperty<Incidences, Entity>(topo());
        recompute();
    } else {
        incident_ = {};
    }
}

} // namespace OpenVolumeMesh
