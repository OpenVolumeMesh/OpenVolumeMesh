#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT CellHalfEdgeIterImpl : public BaseCirculator<CellHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<CellHandle, HalfEdgeHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellHalfEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellHalfEdgeIterImpl& operator++();
    CellHalfEdgeIterImpl& operator--();

private:
    std::vector<HalfEdgeHandle> halfedges_;
    size_t cur_index_;
};



} // namespace OpenVolumeMesh::detail
