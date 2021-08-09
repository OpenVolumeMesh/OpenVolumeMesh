#include <OpenVolumeMesh/Core/detail/HalfFaceCellIncidence.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT_impl.hh>

#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

#include <algorithm>
#include <stdexcept>


namespace OpenVolumeMesh {

//template class IncidencesT<Entity::Vertex, Entity::HalfEdge>;

template<typename Derived>
void HalfFaceCellIncidence<Derived>::add_cell(CellHandle _ch, const OpenVolumeMeshCell &_cell)
{
    if (!enabled()) return;
    for (const auto &hfh: _cell.halffaces()) {
        incident(hfh) = _ch;
    }
    // TODO: if we have eh->hfh incidences, reorder halffaces
}

template<typename Derived>
void HalfFaceCellIncidence<Derived>::delete_cell(CellHandle _ch, const OpenVolumeMeshCell &_cell)
{
    if (!enabled()) return;
    for (const auto &hfh: _cell.halffaces()) {
        if (incident(hfh) == _ch) {
            incident(hfh) = CellHandle{};
        }
    }

}

template<typename Derived>
void HalfFaceCellIncidence<Derived>::swap(CellHandle _h1, CellHandle _h2) {
    if(!enabled()) return;
    throw std::runtime_error("unimplemented");


#if 0
    // correct pointers to those cells
    std::vector<HalfFaceHandle> hfhs1 = c1.halffaces();
    for (unsigned int i = 0; i < hfhs1.size(); ++i)
    {
        HalfFaceHandle hfh = hfhs1[i];
        if (incident_cell_per_hf_[hfh.idx()] == _h1)
            incident_cell_per_hf_[hfh.idx()] = _h2;
    }

    std::vector<HalfFaceHandle> hfhs2 = c2.halffaces();
    for (unsigned int i = 0; i < hfhs2.size(); ++i)
    {
        HalfFaceHandle hfh = hfhs2[i];
        if (incident_cell_per_hf_[hfh.idx()] == _h2)
            incident_cell_per_hf_[hfh.idx()] = _h1;
    }
#endif
}

template<typename Derived>
void HalfFaceCellIncidence<Derived>::recompute()
{
    for (const auto ch: topo()->cells())
    {
        for (const auto hfh: topo()->cell_halffaces(ch)) {
#ifndef NDEBUG
            if(incident(hfh).is_valid()) {
                std::cerr << "compute_face_bottom_up_incidences(): Detected non-three-manifold configuration!" << std::endl;
                std::cerr << "Connectivity probably won't work." << std::endl;
                continue;
            }
#endif
            incident(hfh) = ch;
        }
    }
    // TODO: if we have eh->hfh incidences, reorder halffaces
}

// instantiate for the only class that derives from this:
template class IncidencesT<TopologyKernel, Entity::HalfFace, CellHandle>;
template class HalfFaceCellIncidence<TopologyKernel>;

} // namespace OpenVolumeMesh
