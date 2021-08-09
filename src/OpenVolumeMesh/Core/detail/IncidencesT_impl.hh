#pragma once
#include <OpenVolumeMesh/Core/TopologyKernel.hh>


namespace OpenVolumeMesh {

template<typename Entity, typename _Incidences>
_Incidences &
IncidencesT<Entity, _Incidences>::
incident(Handle _h)
{
    assert(enabled_);
    assert(_h.uidx() < incident_.size());
    return incident_[_h.uidx()];
}

template<typename Entity, typename _Incidences>
_Incidences const &
IncidencesT<Entity, _Incidences>::
incident(Handle _h) const
{
    assert(enabled_);
    assert(_h.uidx() < incident_.size());
    return incident_[_h.uidx()];
}

template<typename Entity, typename _Incidences>
void
IncidencesT<Entity, _Incidences>::
setEnabled(bool enable)
{
    if (enabled_ == enable)
        return;
    enabled_ = enable;
    if (enable) {
        assert(incident_.empty());
        recompute();
    } else {
        clear();
    }
}

template<typename Entity, typename _Incidences>
void
IncidencesT<Entity, _Incidences>::
clear() {
    // TODO: once we store incident_ as prop, we do not need clear()
    // anymore
    incident_.clear();
}

} // namespace OpenVolumeMesh
