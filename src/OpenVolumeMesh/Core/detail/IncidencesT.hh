#pragma once

#include <vector>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

namespace OpenVolumeMesh {

class TopologyKernel;

template<typename SubEntity, typename _Incidences>
class IncidencesT
{
    static_assert(is_entity<SubEntity>::value);
public:
    using SubHandle = HandleT<SubEntity>;
    using Incidences = _Incidences;

    IncidencesT() = default;

#if 0
    IncidencesT(IncidencesT const&) = default;
    IncidencesT(IncidencesT &&) = default;
    IncidencesT& operator=(IncidencesT &&) = default;
    IncidencesT& operator=(IncidencesT const &) = default;
#endif

    bool enabled() const {return enabled_;}
    void setEnabled(bool enable);

    Incidences const& incident(SubHandle _h) const;

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
    void remove(SubHandle _h) {
        if(!enabled_) return;
        incident(_h).clear();
    }
#endif
    // TODO: no need to do this if we keep our vector as prop :)
    void swap(SubHandle _h1, SubHandle _h2) {
        if(!enabled_) return;
        std::swap(incident(_h1), incident(_h2));
    }
    void clear();

    Incidences & incident(SubHandle _h);

    bool valid(SubHandle vh) const;
    virtual void recompute() = 0;

private:
    std::vector<Incidences> incident_;
    bool enabled_ = false;
};


} // namespace OpenVolumeMesh
