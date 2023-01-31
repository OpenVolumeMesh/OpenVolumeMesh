#include <OpenVolumeMesh/Unstable/Topology/TriangleTopology.hh>

namespace OpenVolumeMesh {

TriangleTopology::TriangleTopology(
        const TetrahedralMeshTopologyKernel &_mesh,
        HFH _hfh,
        VH _a)
    : hfh_(_hfh)
{
    int idx = 0;
    for (const auto heh: _mesh.halfface_halfedges(_hfh, 2))
    {
        if (idx == 0 && _mesh.from_vertex_handle(heh) != _a) {
            continue;
        }
        heh_[idx] = heh;
        vh_[idx] = _mesh.from_vertex_handle(heh);
        ++idx;
        if (idx == 3) {
            break;
        }
    }
    assert(idx == 3);
}
TriangleTopology::TriangleTopology(
        const TetrahedralMeshTopologyKernel &_mesh,
        HFH _hfh)
    : hfh_(_hfh)
{
    int idx = 0;
    for (const auto heh: _mesh.halfface_halfedges(_hfh))
    {
        heh_[idx] = heh;
        vh_[idx] = _mesh.from_vertex_handle(heh);
        ++idx;
    }
    assert(idx == 3);
}

bool TriangleTopology::operator==(TriangleTopology const& other) const
{
    // std::array::operator== is not constexpr before C++20, compare elements individually:
    //return vh_[A] == other.vh_[A] && .....
    return vh_ == other.vh_
       && heh_ == other.heh_;
}

} // namespace OpenVolumeMesh
