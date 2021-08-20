#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh::detail {

class OVM_EXPORT CellFaceIterImpl : public BaseCirculator<CellHandle, FaceHandle> {
public:

    typedef BaseCirculator<CellHandle, FaceHandle> BaseIter;
    typedef CellHandle CenterEntityHandle;

    CellFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps = 1);

    CellFaceIterImpl& operator++();
    CellFaceIterImpl& operator--();

private:
    std::vector<HalfFaceHandle>::const_iterator hf_iter_;
};


} // namespace OpenVolumeMesh::detail
