#pragma once
#include <OpenVolumeMesh/Core/TopologyKernel.hh>


namespace OpenVolumeMesh {

template<typename Derived, typename Entity, typename _Incidences>
IncidencesT<Derived, Entity, _Incidences>&
IncidencesT<Derived, Entity, _Incidences>::
operator=(IncidencesT<Derived, Entity, _Incidences> const &other)
{
    incident_ = other.incident_;
    if (enabled()) {
        incident_->attach_to(topo());
    }
    return *this;
}

template<typename Derived, typename Entity, typename _Incidences>
IncidencesT<Derived, Entity, _Incidences>&
IncidencesT<Derived, Entity, _Incidences>::
operator=(IncidencesT<Derived, Entity, _Incidences> &&other)
{
    incident_ = std::move(other.incident_);
    if (enabled()) {
        incident_->attach_to(topo());
    }
    return *this;

}
template<typename Derived, typename Entity, typename _Incidences>
IncidencesT<Derived, Entity, _Incidences>::
IncidencesT(IncidencesT<Derived, Entity, _Incidences> &&other)
{
    *this = std::move(other);
}
template<typename Derived, typename Entity, typename _Incidences>
IncidencesT<Derived, Entity, _Incidences>::
IncidencesT(IncidencesT<Derived, Entity, _Incidences> const &other)
{
    *this = other;
}



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
_Incidences &
IncidencesT<Derived, Entity, _Incidences>::
incident_mutable(Handle _h) const
{
    assert(enabled());
    assert(_h.uidx() < incident_->size());
    return (*incident_)[_h];
}

template<typename Derived, typename Entity, typename _Incidences>
void
IncidencesT<Derived, Entity, _Incidences>::
set_enabled(bool enable)
{
    if (enabled() == enable)
        return;
    if (enable) {
        incident_.emplace(topo(), "bottom-up incidences");
        recompute();
    } else {
        incident_.reset();
    }
}

template<typename Derived, typename Entity, typename _Incidences>
void IncidencesT<Derived, Entity, _Incidences>::deleted(Handle _h)
{
    if (!enabled()) return;
    incident(_h) = {};
}


} // namespace OpenVolumeMesh
