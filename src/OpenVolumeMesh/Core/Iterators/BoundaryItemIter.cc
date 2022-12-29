#include <OpenVolumeMesh/Core/Iterators/BoundaryItemIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

template <class Iter, class Handle>
size_t BoundaryItemIter<Iter, Handle>::n_items() const {
    return BaseIter::mesh()->template n<typename Handle::EntityTag>();
}

template <>
bool BoundaryItemIter<VertexIter, VertexHandle>::has_incidences() const {
    return BaseIter::mesh()->has_full_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<EdgeIter, EdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<FaceIter, FaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
bool BoundaryItemIter<CellIter, CellHandle>::has_incidences() const {
    return true;
}

template class OVM_EXPORT BoundaryItemIter<VertexIter, VertexHandle>;
template class OVM_EXPORT BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>;
template class OVM_EXPORT BoundaryItemIter<EdgeIter, EdgeHandle>;
template class OVM_EXPORT BoundaryItemIter<HalfFaceIter, HalfFaceHandle>;
template class OVM_EXPORT BoundaryItemIter<FaceIter, FaceHandle>;
template class OVM_EXPORT BoundaryItemIter<CellIter, CellHandle>;


} // namespace OpenVolumeMesh
