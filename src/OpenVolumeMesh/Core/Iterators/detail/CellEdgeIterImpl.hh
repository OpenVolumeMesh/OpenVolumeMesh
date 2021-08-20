#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT CellEdgeIterImpl : public BaseCirculator<CellHandle, EdgeHandle> {
public:

    typedef BaseCirculator<CellHandle, EdgeHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellEdgeIterImpl& operator++();
    CellEdgeIterImpl& operator--();

private:
    std::vector<EdgeHandle> edges_;
    size_t cur_index_;
};



} // namespace OpenVolumeMesh::detail
