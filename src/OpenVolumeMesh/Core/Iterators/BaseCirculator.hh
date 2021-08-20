#pragma once

#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

namespace OpenVolumeMesh {

template <
class IH /*  Input handle type */,
class OH /* Output handle type */>
class BaseCirculator : public BaseIterator<OH> {
public:

    typedef BaseIterator<OH> BaseIter;

    BaseCirculator(const TopologyKernel* _mesh, const IH& _ih, const OH& _oh, int _max_laps = 1) :
        BaseIter(_mesh, _oh),
        lap_(0),
        max_laps_(_max_laps),
        ref_handle_(_ih)
    {}

    BaseCirculator(const TopologyKernel* _mesh, const IH& _ih, int _max_laps = 1) :
        BaseIter(_mesh, OH()),
        lap_(0),
        max_laps_(_max_laps),
        ref_handle_(_ih)
    {}

    // STL compliance (needs to have default constructor)
    BaseCirculator() :
        BaseIter(),
        lap_(0),
        max_laps_(1)
    {}
    BaseCirculator(const BaseCirculator& _c) = default;

    virtual ~BaseCirculator() = default;

    bool operator== (const BaseCirculator& _c) const {
        return (BaseIter::operator==(_c) &&
                this->lap() == _c.lap() &&
                this->ref_handle() == _c.ref_handle());
    }
    bool operator!= (const BaseCirculator& _c) const {
        return !this->operator==(_c);
    }

    bool operator< (const BaseCirculator& _c) const {
        if (lap_ == _c.lap_)
            return BaseIter::operator<(_c);
        else
            return lap_ < _c.lap_;
    }

    BaseCirculator& operator=(const BaseCirculator& _c) = default;

    const IH& ref_handle() const {
        return ref_handle_;
    }

    void lap(int _lap) {
        lap_ = _lap;
    }
    int lap() const {
        return lap_;
    }

    void max_laps(int _max_laps) {
        max_laps_ = _max_laps;
    }
    int max_laps() const {
        return max_laps_;
    }

protected:
    int lap_;
    int max_laps_;
    IH ref_handle_;

};


} // namespace OpenVolumeMesh
