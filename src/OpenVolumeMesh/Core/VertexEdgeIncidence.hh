#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

class OVM_EXPORT VertexEdgeIncidence
        : public IncidencesT<Entity::Vertex, Entity::HalfEdge>
{
public:
protected:
    void add_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
    void delete_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge);
private:
    void recompute() override;

private:
    bool enabled_;
};

} // namespace OpenVolumeMesh

