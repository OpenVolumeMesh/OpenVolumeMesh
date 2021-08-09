#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Derived> // CRDT
class OVM_EXPORT VertexHalfEdgeIncidence
        : public IncidencesT<Derived, Entity::Vertex, std::vector<HalfEdgeHandle>>
{
protected:
    using Parent = IncidencesT<Derived, Entity::Vertex, std::vector<HalfEdgeHandle>>;
    using Incidences = typename Parent::Incidences;
    using Parent::incident;
    using Parent::enabled;
    void add_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
    void delete_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
    void swap(EdgeHandle _h1, EdgeHandle _h2);
private:
    using Parent::topo;
    void recompute() override;
};

} // namespace OpenVolumeMesh

