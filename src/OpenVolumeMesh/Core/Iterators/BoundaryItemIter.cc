#include <OpenVolumeMesh/Core/Iterators/BoundaryItemIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

template <>
size_t BoundaryItemIter<VertexIter, VertexHandle>::n_items() const {
    return BaseIter::mesh()->n_vertices();
}

template <>
size_t BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::n_items() const {
    return BaseIter::mesh()->n_halfedges();
}

template <>
size_t BoundaryItemIter<EdgeIter, EdgeHandle>::n_items() const {
    return BaseIter::mesh()->n_edges();
}

template <>
size_t BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::n_items() const {
    return BaseIter::mesh()->n_halffaces();
}

template <>
size_t BoundaryItemIter<FaceIter, FaceHandle>::n_items() const {
    return BaseIter::mesh()->n_faces();
}

template <>
size_t BoundaryItemIter<CellIter, CellHandle>::n_items() const {
    return BaseIter::mesh()->n_cells();
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


} // namespace OpenVolumeMesh
