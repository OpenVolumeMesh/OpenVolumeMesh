#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Derived> // CRDT
class OVM_EXPORT VertexHalfEdgeIncidence
        : public IncidencesT<Entity::Vertex, std::vector<HalfEdgeHandle>>
{
public:
    using IncidencesT::IncidencesT;
protected:
    void add_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
    void delete_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
private:
    const Derived *topo() const {return static_cast<const Derived*>(this);}
    void recompute() override;
};

} // namespace OpenVolumeMesh

