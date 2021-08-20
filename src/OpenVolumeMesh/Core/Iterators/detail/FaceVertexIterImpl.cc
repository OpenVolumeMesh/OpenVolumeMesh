#include <OpenVolumeMesh/Core/Iterators/detail/FaceVertexIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

FaceVertexIterImpl::FaceVertexIterImpl(const FaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfFaceVertexIter(_mesh->halfface_handle(_ref_h, 0), _mesh, _max_laps) {}


} // namespace OpenVolumeMesh::detail
