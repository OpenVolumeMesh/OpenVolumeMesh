#include <OpenVolumeMesh/Core/detail/VertexHalfEdgeIncidence.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT_impl.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

#include <algorithm>
#include <stdexcept>
#include <unordered_set>


namespace OpenVolumeMesh {

//template class IncidencesT<Entity::Vertex, Entity::HalfEdge>;

template<typename Derived>
void VertexHalfEdgeIncidence<Derived>::recompute()
{
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
        if (!vh.is_valid()) return;
        Incidences &vec = incident(vh);
        // PERF: we can swap in the last element and resize instead of using remove
        vec.erase(std::remove( vec.begin(), vec.end(),
                      heh),
                  vec.end());
    };

    rm_v_heh(_edge.from_vertex(), topo()->halfedge_handle(_eh, 0));
    rm_v_heh(_edge.to_vertex(),   topo()->halfedge_handle(_eh, 1));
}

template<typename Derived>
void VertexHalfEdgeIncidence<Derived>::swap(EdgeHandle _h1, EdgeHandle _h2)
{
    if (!enabled()) return;
    if (_h1 == _h2) return;

    auto swap_indices = [this, _h1, _h2](VertexHandle vh)
    {
        for (auto &heh: incident(vh)) {
            if      (heh.full() == _h1) heh = _h2.half(heh.subidx());
            else if (heh.full() == _h2) heh = _h1.half(heh.subidx());
        }
    };

    std::unordered_set<VertexHandle> processed;
    processed.reserve(4);


    for (const auto eh: {_h1, _h2}) {
        const auto &edge = topo()->edge(eh);
        for (const auto vh: {edge.from_vertex(), edge.to_vertex()}) {
            if (!vh.is_valid() || processed.find(vh) != processed.end())
                continue;
            processed.insert(vh);
            swap_indices(vh);
        }
    }


#if 0
    const auto &e1 = topo()->edge(_h1);
    const auto &e2 = topo()->edge(_h2);

    auto heh1 = topo()->halfedge_handle(_h1, 0);
    auto heh2 = topo()->halfedge_handle(_h2, 0);

    if (e1.from_vertex().is_valid()) {
        for (auto &heh: incident(e1.from_vertex())) {
            if (heh == heh1) {
                heh = heh2;
            } else if (heh == heh2) {
                heh = heh1;
            }
        }
    }
    if (e1.from_vertex() != e2.from_vertex() && e2.from_vertex().is_valid()) {
        for (auto &heh: incident(e2.from_vertex())) {
            if (heh == heh2) {
                heh = heh1;
            }
        }
    }
    auto opp1 = topo()->halfedge_handle(_h1, 1);
    auto opp2 = topo()->halfedge_handle(_h2, 1);

    if (e1.to_vertex().is_valid()) {
        for (auto &heh: incident(e1.to_vertex())) {
            if (heh == opp1) {
                heh = opp2;
            } else if (heh == opp2) {
                heh = opp1;
            }
        }
    }
    if (e1.to_vertex() != e2.to_vertex() && e2.to_vertex().is_valid()) {
        for (auto &heh: incident(e2.to_vertex())) {
            if (heh == opp2) {
                heh = opp1;
            }
        }
    }
#endif
}


// instantiate for the only class that derives from this:
template class IncidencesT<TopologyKernel, Entity::Vertex, std::vector<HalfEdgeHandle>>;
template class VertexHalfEdgeIncidence<TopologyKernel>;

} // namespace OpenVolumeMesh
