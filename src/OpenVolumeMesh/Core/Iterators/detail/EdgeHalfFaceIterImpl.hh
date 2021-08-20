#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT EdgeHalfFaceIterImpl : public BaseCirculator<EdgeHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<EdgeHandle, HalfFaceHandle> BaseIter;
    typedef EdgeHandle CenterEntityHandle;

    EdgeHalfFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    EdgeHalfFaceIterImpl& operator++();
    EdgeHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle> halffaces_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh::detail
