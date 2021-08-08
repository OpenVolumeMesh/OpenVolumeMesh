#include <OpenVolumeMesh/Core/VertexEdgeIncidence.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/IncidencesT_impl.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

#include <algorithm>


namespace OpenVolumeMesh {

template class IncidencesT<Entity::Vertex, Entity::HalfEdge>;

void VertexEdgeIncidence::recompute()
{
    // TODO
}

void VertexEdgeIncidence::add_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge)
{
    if (enabled_) return;
    incident(_edge.from_vertex()).push_back(topo()->halfedge_handle(_eh, 0));
    incident(_edge.to_vertex()).push_back(topo()->halfedge_handle(_eh, 1));
}

void VertexEdgeIncidence::delete_edge(EdgeHandle _eh, const OpenVolumeMeshEdge &_edge)
{
    if (enabled_) return;
    auto rm_v_heh = [this](VertexHandle vh, HalfEdgeHandle heh)
    {
        Incidences &vec = incident(vh);
        vec.erase(std::remove( vec.begin(), vec.end(),
                      heh),
                  vec.end());
    };

    rm_v_heh(_edge.from_vertex(), topo()->halfedge_handle(_eh, 0));
    rm_v_heh(_edge.to_vertex(),   topo()->halfedge_handle(_eh, 1));
}

} // namespace OpenVolumeMesh
