#pragma once

#include <vector>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename Entity, typename _Incidences>
class IncidencesT
{
    static_assert(is_entity<Entity>::value);
public:
    using Handle = HandleT<Entity>;
    using Incidences = _Incidences;

    IncidencesT() = default;

    bool enabled() const {return enabled_;}
    void setEnabled(bool enable);

    Incidences const& incident(Handle _h) const;

protected:
    void reserve(size_t n) {
        if(!enabled_) return;
        incident_.reserve(n);
    }
    void resize(size_t n) {
        if(!enabled_) return;
        incident_.resize(n, Incidences{});
    }
#if 0
    void remove(Handle _h) {
        if(!enabled_) return;
        incident(_h).clear();
    }
#endif
    // TODO: no need to do this if we keep our vector as prop :)
    void swap(Handle _h1, Handle _h2) {
        if(!enabled_) return;
        std::swap(incident(_h1), incident(_h2));
    }
    void clear();

    Incidences & incident(Handle _h);

    bool valid(Handle vh) const;
    virtual void recompute() = 0;

private:
    std::vector<Incidences> incident_;
    bool enabled_ = false;
};


} // namespace OpenVolumeMesh
