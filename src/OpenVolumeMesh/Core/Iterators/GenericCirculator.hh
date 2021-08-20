#pragma once
namespace OpenVolumeMesh {

class TopologyKernel;

template <class CirculatorImpl>
class GenericCirculator : public CirculatorImpl {
public:

    GenericCirculator(const typename CirculatorImpl::CenterEntityHandle& _ref_h,
                      const TopologyKernel* _mesh, int _max_laps = 1) :
        CirculatorImpl(_ref_h, _mesh, _max_laps) {}

    GenericCirculator& operator++() {
        CirculatorImpl::operator++();
        return *this;
    }

    GenericCirculator& operator--() {
        CirculatorImpl::operator--();
        return *this;
    }

    // Post increment/decrement operator
    GenericCirculator operator++(int) {
        GenericCirculator cpy = *this;
        ++(*this);
        return cpy;
    }
    GenericCirculator operator--(int) {
        GenericCirculator cpy = *this;
        --(*this);
        return cpy;
    }
    GenericCirculator operator+(int _n) {
        GenericCirculator cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    GenericCirculator operator-(int _n) {
        GenericCirculator cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    GenericCirculator& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    GenericCirculator& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

};

} // namespace OpenVolumeMesh
