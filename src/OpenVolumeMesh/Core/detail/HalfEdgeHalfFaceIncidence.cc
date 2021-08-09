#include <OpenVolumeMesh/Core/detail/HalfEdgeHalfFaceIncidence.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT_impl.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

#include <algorithm>


namespace OpenVolumeMesh {

template class IncidencesT<Entity::Vertex, Entity::HalfEdge>;

template<typename Derived>
void HalfEdgeHalfFaceIncidence<Derived>::add_face(FaceHandle _fh, const OpenVolumeMeshFace &_face)
{
    auto hfh0 = topo()->halfface_handle(_fh, 0);
    auto hfh1 = topo()->halfface_handle(_fh, 1);
    for (const auto &heh: _face.halfedges()) {
        incident(heh).push_back(hfh0);
        auto opp = topo()->opposite_halfedge_handle(heh);
        incident(opp).push_back(hfh1);
    }
    // we added a single face, so the configuration cannot be manifold.
    // no need to call reorder_incident_halffaces()
}

template<typename Derived>
void HalfEdgeHalfFaceIncidence<Derived>::delete_face(FaceHandle _fh, const OpenVolumeMeshFace &_face)
{
    auto rm = [this](HalfEdgeHandle heh, HalfFaceHandle hfh)
    {
        Incidences &vec = incident(heh);
        vec.erase(std::remove( vec.begin(), vec.end(),
                      hfh),
                  vec.end());
    };

    auto hfh0 = topo()->halfface_handle(_fh, 0);
    auto hfh1 = topo()->halfface_handle(_fh, 1);
    for (const auto &heh: _face.halfedges()) {
        rm(heh, hfh0);
        auto opp = topo()->opposite_halfedge_handle(heh);
        rm(opp, hfh1);
    }

}

template<typename Derived>
void HalfEdgeHalfFaceIncidence<Derived>::reorder_halffaces(const EdgeHandle &_eh)
{
    /* Put halffaces in clockwise order via the
     * same cell property which now exists.
     * Note, this only works for manifold configurations though.
     * Proceed as follows: Pick one starting halfface. Assuming
     * that all halfface normals point into the incident cell,
     * we find the adjacent halfface within the incident cell
     * along the considered halfedge. We set the found halfface
     * to be the one to be processed next. If we reach an outside
     * region, we try to go back from the starting halfface in reverse
     * order. If the complex is properly connected (the pairwise
     * intersection of two adjacent 3-dimensional cells is always
     * a 2-dimensional entity, namely a facet), such an ordering
     * always exists and will be found. If not, a correct order
     * can not be given and, as a result, the related iterators
     * will address the related entities in an arbitrary fashion.
     */

    assert(topo()->is_valid(_eh));
    HalfEdgeHandle heh = topo()->halfedge_handle(_eh, 0);
    auto &incident_hfs = incident(heh);

    const size_t n_hfs = incident_hfs.size();

    if(n_hfs < 2)
        return;

    std::vector<HalfFaceHandle> new_halffaces;
    new_halffaces.reserve(n_hfs);

    // Start with one incident halfface and go into the first direction
    auto start_hf = incident_hfs.front();
    auto cur_hf = start_hf;

    do {
        new_halffaces.push_back(cur_hf);
        if (new_halffaces.size() > incident_hfs.size()) {
            //std::cerr << "reorder_incident_halffaces(" << _eh.idx() << "): weird topology, aborting." << std::endl;
            return;
        };

        auto ch = topo()->incident_cell(cur_hf);
        if (!ch.is_valid() || topo()->is_deleted(ch))
            break;

        cur_hf = topo()->adjacent_halfface_in_cell(cur_hf, heh);
        if(!cur_hf.is_valid()) {
            return;
        }
        cur_hf = topo()->opposite_halfface_handle(cur_hf);

    } while (cur_hf != start_hf);

    // First direction has terminated
    // If new_halffaces has the same size as old (unordered)
    // vector of incident halffaces, we are done here
    // If not, try the other way round
    // (this must be a boundary edge)
    if(new_halffaces.size() != incident_hfs.size()) {

        cur_hf = start_hf;

        while(true) {
            cur_hf = topo()->opposite_halfface_handle(cur_hf);

            auto ch = topo()->incident_cell(cur_hf);
            if (!ch.is_valid() || topo()->is_deleted(ch))
                break;

            cur_hf = topo()->adjacent_halfface_in_cell(cur_hf, heh);
            if(!cur_hf.is_valid()) {
                return;
            }

            // TODO PERF: just move everything we already have to the end *once* and fill backwards
            new_halffaces.insert(new_halffaces.begin(), cur_hf);
            if(new_halffaces.size() > incident_hfs.size()) {
                //std::cerr << "reorder_incident_halffaces(" << _eh.idx() << ") #2: weird topology, aborting" << std::endl;
                return;
            }
        }
    }

    // Everything worked just fine, set the new ordered vector
    if(new_halffaces.size() == incident_hfs.size()) {
        incident_hfs = std::move(new_halffaces);
        auto &opp_incident = incident(topo()->opposite_halfedge_handle(heh));
        // update incident halffaces of the opposite halfedge:
        std::transform(incident_hfs.rbegin(), incident_hfs.rend(),
                       opp_incident.begin(),
                       Derived::opposite_halfface_handle);
    }
#if 0
    else {
        std::cerr << "reorder_incident_halffaces: found " << new_halffaces.size() << " of " << incident_hfs.size()
            << " incident halffaces, likely the edge has more than one boundary! Currently not supported, not reordering." << std::endl;
        // TODO FIXME: we should support this case.
    }
#endif


}

template<typename Derived>
void HalfEdgeHalfFaceIncidence<Derived>::recompute()
{
    resize(topo()->n_halfedges());
    // count
    std::vector<int> n_faces_per_edge(topo()->n_edges(), 0);
    for (const auto &fh: topo()->faces()) {
        for (const auto &heh: topo()->face(fh).halfedges()) {
            ++n_faces_per_edge[topo()->edge_handle(heh).idx()];
        }
    }

    // reserve
    for (const auto &eh: topo()->edges()) {
        auto n = n_faces_per_edge[eh.idx()];
        auto heh = topo()->halfedge_handle(eh, 0);
        auto opp = topo()->halfedge_handle(eh, 1);
        incident(heh).clear();
        incident(heh).reserve(n);
        incident(opp).clear();
        incident(opp).reserve(n);
    }

    // store
    for (const auto &fh: topo()->faces()) {
        auto hfh0 = topo()->halfface_handle(fh, 0);
        auto hfh1 = topo()->halfface_handle(fh, 1);
        for(const auto &heh: topo()->face_halfedges(fh)) {
            incident(heh).push_back(hfh0);
            auto opp = topo()->opposite_halfedge_handle(heh);
            incident(opp).push_back(hfh1);
        }
    }
    // TODO: if we have hf->ch incidences, reorder halffaces
}

// instantiate for the only class that derives from this:
template class HalfEdgeHalfFaceIncidence<TopologyKernel>;

} // namespace OpenVolumeMesh
