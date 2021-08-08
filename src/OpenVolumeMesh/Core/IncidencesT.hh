#pragma once

#include <vector>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename SubEntity, typename SuperEntity>
class IncidencesT
{
    static_assert(is_entity<SubEntity>::value);
    static_assert(is_entity<SuperEntity>::value);
    using SubHandle = HandleT<SubEntity>;
    using SuperHandle = HandleT<SuperEntity>;
public:
    using Incidences = std::vector<SuperHandle>;
    bool enabled() const {return enabled_;}
    void setEnabled(bool enable);
    void clear();
    size_t valence(SubHandle _h) const {
        assert(enabled_);
        return incident(_h).size();
    }
    Incidences const& incident(SubHandle _h) const {
        assert(enabled_);
        assert(topo()->is_valid(_h));
        return incident_[_h.uidx()];
    }
    void reserve(size_t n) {
        if(!enabled_) return;
        incident_.reserve(n);
    }
    void resize(size_t n) {
        if(!enabled_) return;
        incident_.resize(n);
    }
    void remove(SubHandle _h) {
        if(!enabled_) return;
        incident(_h).clear();
    }
    // TODO: no need to do this if we keep our vector as prop :)
    void swap(SubHandle _h1, SubHandle _h2) {
        if(!enabled_) return;
        std::swap(incident(_h1), incident(_h2));
    }
protected:
    Incidences & incident(SubHandle _h) {
        assert(enabled_);
        assert(topo()->is_valid(_h));
        return incident_[_h.uidx()];
    }

    TopologyKernel const* topo() const;
    bool valid(SubHandle vh) const;
    virtual void recompute() = 0;

private:
    bool enabled_;
    std::vector<Incidences> incident_;
};

} // namespace OpenVolumeMesh
