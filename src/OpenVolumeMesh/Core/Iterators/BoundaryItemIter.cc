#include <OpenVolumeMesh/Core/Iterators/BoundaryItemIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

template <class Iter, class Handle>
OVM_EXPORT
size_t BoundaryItemIter<Iter, Handle>::n_items() const {
    return BaseIter::mesh()->template n<typename Handle::EntityTag>();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<VertexIter, VertexHandle>::has_incidences() const {
    return BaseIter::mesh()->has_full_bottom_up_incidences();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<EdgeIter, EdgeHandle>::has_incidences() const {
    const TopologyKernel *mesh = BaseIter::mesh();
    return mesh->has_edge_bottom_up_incidences() && mesh->has_face_bottom_up_incidences();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<FaceIter, FaceHandle>::has_incidences() const {
    return BaseIter::mesh()->has_face_bottom_up_incidences();
}

template <>
OVM_EXPORT
bool BoundaryItemIter<CellIter, CellHandle>::has_incidences() const {
    return true;
}

template class BoundaryItemIter<VertexIter, VertexHandle>;
template class BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>;
template class BoundaryItemIter<EdgeIter, EdgeHandle>;
template class BoundaryItemIter<HalfFaceIter, HalfFaceHandle>;
template class BoundaryItemIter<FaceIter, FaceHandle>;
template class BoundaryItemIter<CellIter, CellHandle>;

} // namespace OpenVolumeMesh
