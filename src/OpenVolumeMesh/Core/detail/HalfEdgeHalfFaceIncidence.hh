#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Derived> // CRDT
class OVM_EXPORT HalfEdgeHalfFaceIncidence
        : public IncidencesT<Derived, Entity::HalfEdge, std::vector<HalfFaceHandle>>
{
protected:
    using Parent = IncidencesT<Derived, Entity::HalfEdge, std::vector<HalfFaceHandle>>;
    using Incidences = typename Parent::Incidences;
    using Parent::incident;
    using Parent::enabled;

    void add_face(FaceHandle _fh, OpenVolumeMeshFace const &_face);
    void delete_face(FaceHandle _fh, OpenVolumeMeshFace const &_face);

    /// Try to reorder the vector of incident half-faces to cyclic order.
    /// This requires HalfFace-Cell incidences and a manifold configuration around the edge.
    void reorder_halffaces(const EdgeHandle &eh);
    void invalidate_order(const EdgeHandle &) { /* TODO */ };
    void invalidate_order(const FaceHandle &) { /* TODO */ };
    void invalidate_order(const CellHandle &) { /* TODO */ };

    void swap(FaceHandle _h1, FaceHandle _h2);
private:
    using Parent::topo;
    //std::vector<bool> order_valid_;
    void recompute() override;
    // TODO: maybe keep a prop<bool, EH> "ordered" to invalidate?
};

} // namespace OpenVolumeMesh

