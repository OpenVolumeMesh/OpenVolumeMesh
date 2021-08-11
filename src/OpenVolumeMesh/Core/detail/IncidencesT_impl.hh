#pragma once
#include <OpenVolumeMesh/Core/TopologyKernel.hh>


namespace OpenVolumeMesh {

template<typename Derived, typename Entity, typename _Incidences>
_Incidences &
IncidencesT<Derived, Entity, _Incidences>::
incident(Handle _h)
{
    assert(enabled());
    assert(_h.uidx() < incident_.size());
    return incident_[_h.idx()];
}

template<typename Derived, typename Entity, typename _Incidences>
_Incidences const &
IncidencesT<Derived, Entity, _Incidences>::
incident(Handle _h) const
{
    assert(enabled());
    assert(_h.uidx() < incident_.size());
    return incident_[_h.idx()];
}

template<typename Derived, typename Entity, typename _Incidences>
_Incidences &
IncidencesT<Derived, Entity, _Incidences>::
incident_mutable(Handle _h) const
{
    assert(enabled());
    assert(_h.uidx() < incident_.size());
    return incident_[_h.idx()];
}

template<typename Derived, typename Entity, typename _Incidences>
void
IncidencesT<Derived, Entity, _Incidences>::
set_enabled(bool enable)
{
    if (enabled() == enable)
        return;
    enabled_ = enable;
    if (enabled_) {
        resize();
        recompute();
    } else {
        incident_.clear();
        incident_.shrink_to_fit();
    }
}

template<typename Derived, typename Entity, typename _Incidences>
void IncidencesT<Derived, Entity, _Incidences>::deleted(Handle _h)
{
    if (!enabled()) return;
    incident(_h) = {};
}


template<typename Derived, typename Entity, typename _Incidences>
void IncidencesT<Derived, Entity, _Incidences>::resize()
{
    if (!enabled()) return;
    incident_.resize(topo()->template n<Entity>());
}

template<typename Derived, typename _Entity, typename _Incidences>
void IncidencesT<Derived, _Entity, _Incidences>::swap(Handle _h1, Handle _h2)
{
    if (!enabled()) return;
    std::swap(incident_[_h1.idx()], incident_[_h2.idx()]);
}


} // namespace OpenVolumeMesh
