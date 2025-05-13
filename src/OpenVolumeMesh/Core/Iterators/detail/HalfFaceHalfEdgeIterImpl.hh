#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT HalfFaceHalfEdgeIterImpl : public BaseCirculator<HalfFaceHandle, HalfEdgeHandle> {
public:

    typedef BaseCirculator<HalfFaceHandle, HalfEdgeHandle> BaseIter;
    typedef HalfFaceHandle CenterEntityHandle;

    HalfFaceHalfEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    HalfFaceHalfEdgeIterImpl& operator++();
    HalfFaceHalfEdgeIterImpl& operator--();

private:
    size_t cur_index_;

    const std::vector<HalfEdgeHandle>& f_hehs() const;

    HalfEdgeHandle cur_heh() const;
};


} // namespace OpenVolumeMesh::detail
