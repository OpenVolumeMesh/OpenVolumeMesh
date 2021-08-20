#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT HalfEdgeFaceIterImpl : public BaseCirculator<HalfEdgeHandle, FaceHandle> {
public:

    typedef BaseCirculator<HalfEdgeHandle, FaceHandle> BaseIter;
    typedef HalfEdgeHandle CenterEntityHandle;

    HalfEdgeFaceIterImpl(const HalfEdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfEdgeFaceIterImpl& operator++();
    HalfEdgeFaceIterImpl& operator--();

private:
    std::vector<FaceHandle> faces_;
    size_t cur_index_;
};


} // namespace OpenVolumeMesh::detail
