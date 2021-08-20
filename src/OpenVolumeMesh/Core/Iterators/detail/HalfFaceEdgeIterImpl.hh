#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT HalfFaceEdgeIterImpl : public BaseCirculator<HalfFaceHandle, EdgeHandle> {
public:

    typedef BaseCirculator<HalfFaceHandle, EdgeHandle> BaseIter;
    typedef HalfFaceHandle CenterEntityHandle;

    HalfFaceEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfFaceEdgeIterImpl& operator++();
    HalfFaceEdgeIterImpl& operator--();

private:
    size_t cur_index_;
};


} // namespace OpenVolumeMesh::detail
