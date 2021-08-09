#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/detail/IncidencesT.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Derived> // CRDT
class OVM_EXPORT HalfFaceCellIncidence
        : public IncidencesT<Entity::HalfFace, CellHandle>
{
public:
    using IncidencesT::IncidencesT;
protected:
    void add_cell(CellHandle _ch, OpenVolumeMeshCell const &_cell);
    void delete_cell(CellHandle _ch, OpenVolumeMeshCell const &_cell);

    using IncidencesT::swap;
    void swap(CellHandle _h1, CellHandle _h2);
private:
    const Derived *topo() const {return static_cast<const Derived*>(this);}
    void recompute() override;
};

} // namespace OpenVolumeMesh

