#include <OpenVolumeMesh/Core/Iterators/detail/EdgeCellIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

EdgeCellIterImpl::EdgeCellIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    HalfEdgeCellIter(_mesh->halfedge_handle(_ref_h, 0), _mesh, _max_laps) {}



} // namespace OpenVolumeMesh::detail
