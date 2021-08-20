#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT VertexHalfFaceIterImpl : public BaseCirculator<VertexHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<VertexHandle, HalfFaceHandle> BaseIter;
    typedef VertexHandle CenterEntityHandle;

    VertexHalfFaceIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    VertexHalfFaceIterImpl& operator++();
    VertexHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle> halffaces_;
    size_t cur_index_;
};


} // namespace OpenVolumeMesh::detail
