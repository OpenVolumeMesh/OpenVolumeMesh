#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>

namespace OpenVolumeMesh {

class TetTopology;

class OVM_EXPORT TriangleTopology
{
public:
    enum VertexLabel : uint8_t {
        A, B, C
    };
    enum HalfEdgeLabel : uint8_t {
        AB, BC, CA
    };
    // TODO: template on mesh type - this should work on mixed meshes too.
    // TODO: constructors that takes HFH, HFH+VH, ...
    TriangleTopology(TetrahedralMeshTopologyKernel const &mesh, HFH hfh);
    TriangleTopology(TetrahedralMeshTopologyKernel const &mesh, HFH hfh, VH _a);
    bool operator==(TriangleTopology const& other) const;

    template<VertexLabel I>
    constexpr VH vh() const {return vh_[I];}

    constexpr VH a() const {return vh<A>();}
    constexpr VH b() const {return vh<B>();}
    constexpr VH c() const {return vh<C>();}

    template<HalfEdgeLabel I>
    constexpr HEH heh() const {return heh_[I];}

    constexpr HEH ab() const {return heh<AB>();}
    constexpr HEH bc() const {return heh<BC>();}
    constexpr HEH ca() const {return heh<CA>();}
private:
    constexpr TriangleTopology(HFH _hfh, std::array<VH, 3> _vh, std::array<HEH, 3> _heh)
        : hfh_(_hfh)
          , vh_(std::move(_vh))
          , heh_(std::move(_heh))
    {
    }
    HFH hfh_;
    std::array<VH, 3> vh_;
    std::array<HEH, 3> heh_;
    friend class TetTopology;
};

} // namespace OpenVolumeMesh

