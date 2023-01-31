#pragma once
#include <OpenVolumeMesh/Unstable/Topology/TetTopology.hh>
#include <OpenVolumeMesh/Unstable/Topology/TriangleTopology.hh>

namespace OpenVolumeMesh {

template<TetTopology::HalfFaceLabel HFL>
constexpr TriangleTopology TetTopology::triangle_topology() const
{
    static_assert(has_start(HFL));

    return TriangleTopology(hfh<inner<HFL>()>(),
            {
                vh<hfl_vl<HFL,0>()>(),
                vh<hfl_vl<HFL,1>()>(),
                vh<hfl_vl<HFL,2>()>()
            },
            {
                heh<hfl_hel<HFL,0>()>(),
                heh<hfl_hel<HFL,1>()>(),
                heh<hfl_hel<HFL,2>()>()
            });
}


} // namespace OpenVolumeMesh
