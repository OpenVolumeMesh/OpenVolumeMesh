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
    // We assume there is no common halfface between the cells (implied by manifoldness)

    for (const auto &ch: {_h1, _h2}) {
        for (const auto hfh: topo()->cell_halffaces(ch)) {
            auto &inc = incident(hfh);
            if (inc == _h1) {
                inc = _h2;
            } else if (inc == _h2) {
                inc = _h1;
            }
        }
    }
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
