#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Derived> // CRDT
class OVM_EXPORT HalfFaceCellIncidence
        : public IncidencesT<Derived, Entity::HalfFace, CellHandle>
{
protected:
    using Parent = IncidencesT<Derived, Entity::HalfFace, CellHandle>;
    using Incidences = typename Parent::Incidences;
    using Parent::incident;
    using Parent::enabled;
    void add_cell(CellHandle _ch, OpenVolumeMeshCell const &_cell);
    void delete_cell(CellHandle _ch, OpenVolumeMeshCell const &_cell);

    void swap(CellHandle _h1, CellHandle _h2);
    using Parent::swap;
    using Parent::resize;
private:
    using Parent::topo;
    void recompute() override;
};

} // namespace OpenVolumeMesh

