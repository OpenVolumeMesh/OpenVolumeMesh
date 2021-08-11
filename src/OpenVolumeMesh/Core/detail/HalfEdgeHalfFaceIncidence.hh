#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

// TODO: it suffices to store the incident halffaces for he0, and compute them on-the-fly for he1
template<typename Derived> // CRDT
class OVM_EXPORT HalfEdgeHalfFaceIncidence
        : public IncidencesT<Derived, Entity::HalfEdge, std::vector<HalfFaceHandle>>
{
protected:
    using Parent = IncidencesT<Derived, Entity::HalfEdge, std::vector<HalfFaceHandle>>;
    using Incidences = typename Parent::Incidences;
    using Parent::incident;
    using Parent::incident_mutable;
    using Parent::enabled;

    void add_face(FaceHandle _fh, OpenVolumeMeshFace const &_face);
    void delete_face(FaceHandle _fh, OpenVolumeMeshFace const &_face);

    /// Try to reorder the vector of incident half-faces to cyclic order.
    /// This requires HalfFace-Cell incidences and a manifold configuration around the edge.
    void reorder_halffaces(EdgeHandle) const;

    void invalidate_order(EdgeHandle);
    void invalidate_order(FaceHandle);
    void invalidate_order(CellHandle);
    void ensure_ordered(EdgeHandle) const;

    void swap(FaceHandle _h1, FaceHandle _h2);
    void swap(EdgeHandle _h1, EdgeHandle _h2);

    void set_enabled(bool enable);

    void debug_check(EdgeHandle) const;
    void debug_check(HalfEdgeHandle) const;
    size_t count(HalfEdgeHandle _h) const {return incident(_h).size();};
    void resize();
private:
    using Parent::topo;
    mutable std::vector<bool> edge_ordered_;
    void recompute() override;
};

} // namespace OpenVolumeMesh

