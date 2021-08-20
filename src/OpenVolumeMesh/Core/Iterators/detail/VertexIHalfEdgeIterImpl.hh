#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexOHalfEdgeIter.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT VertexIHalfEdgeIterImpl : public BaseCirculator<VertexHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<VertexHandle, HalfEdgeHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexIHalfEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexIHalfEdgeIterImpl& operator++();
    VertexIHalfEdgeIterImpl& operator--();

private:
    VertexOHalfEdgeIter voh_iter_;
};


} // namespace OpenVolumeMesh::detail
