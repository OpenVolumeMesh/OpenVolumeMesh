#pragma once
#include <OpenVolumeMesh/Core/IncidencesT.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>


namespace OpenVolumeMesh {

template<typename SubEntity, typename SuperEntity>
void
IncidencesT<SubEntity, SuperEntity>::
setEnabled(bool enable)
{
    if (enabled_ == enable)
        return;
    enabled_ = enable;
    if (enable) {
        incident_.resize(topo()->template n<SubEntity>());
        recompute();
    } else {
        clear();
    }
}

template<typename SubEntity, typename SuperEntity>
void
IncidencesT<SubEntity, SuperEntity>::
clear() {
    incident_.clear();
}

template<typename SubEntity, typename SuperEntity>
const TopologyKernel *
IncidencesT<SubEntity, SuperEntity>::
topo() const {
    return static_cast<const TopologyKernel*>(this);
}


} // namespace OpenVolumeMesh
