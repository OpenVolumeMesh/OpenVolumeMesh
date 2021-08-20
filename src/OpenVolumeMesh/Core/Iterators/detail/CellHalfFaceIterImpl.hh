#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT CellHalfFaceIterImpl : public BaseCirculator<CellHandle, HalfFaceHandle> {
public:

    typedef BaseCirculator<CellHandle, HalfFaceHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellHalfFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellHalfFaceIterImpl& operator++();
    CellHalfFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle>::const_iterator hf_begin_, hf_iter_, hf_end_;
};



} // namespace OpenVolumeMesh::detail
