#include <OpenVolumeMesh/Core/detail/VertexHalfEdgeIncidence.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT_impl.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

#include <algorithm>


namespace OpenVolumeMesh {

template class IncidencesT<Entity::Vertex, Entity::HalfEdge>;

template<typename Derived>
void VertexHalfEdgeIncidence<Derived>::recompute()
{
    resize(topo()->n_vertices());
    // count
    std::vector<int> n_edges_per_vertex(topo()->n_vertices(), 0);
    for (const auto &eh: topo()->edges()) {
        auto &e = topo()->edge(eh);
        ++n_edges_per_vertex[e.from_vertex().idx()];
        ++n_edges_per_vertex[e.to_vertex().idx()];
    }
    // reserve
    for (const auto &vh: topo()->vertices()) {
        incident(vh).reserve(n_edges_per_vertex[vh.idx()]);
    }

    // store
    for (const auto &eh: topo()->edges())
    {
        auto &e = topo()->edge(eh);
        incident(e.from_vertex()).push_back(topo()->halfedge_handle(eh, 0));
        incident(e.to_vertex()).push_back(topo()->halfedge_handle(eh, 1));
    }

}

template<typename Derived>
void VertexHalfEdgeIncidence<Derived>::add_edge(EdgeHandle _eh, OpenVolumeMeshEdge const &_edge)
{
    if (!enabled()) return;
    incident(_edge.from_vertex()).push_back(topo()->halfedge_handle(_eh, 0));
    incident(_edge.to_vertex()).push_back(topo()->halfedge_handle(_eh, 1));
}

template<typename Derived>
void VertexHalfEdgeIncidence<Derived>::delete_edge(EdgeHandle _eh, const OpenVolumeMeshEdge &_edge)
{
    if (!enabled()) return;
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


// instantiate for the only class that derives from this:
template class VertexHalfEdgeIncidence<TopologyKernel>;

} // namespace OpenVolumeMesh
