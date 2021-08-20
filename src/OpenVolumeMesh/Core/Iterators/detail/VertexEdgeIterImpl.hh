#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexOHalfEdgeIter.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT VertexEdgeIterImpl : public BaseCirculator<VertexHandle, EdgeHandle> {
public:

    typedef BaseCirculator<VertexHandle, EdgeHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexEdgeIterImpl& operator++();
    VertexEdgeIterImpl& operator--();

private:
    VertexOHalfEdgeIter voh_iter_;
};



} // namespace OpenVolumeMesh::detail
