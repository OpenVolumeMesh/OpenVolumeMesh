#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif


namespace OpenVolumeMesh {

class OVM_EXPORT EdgeIter : public BaseIterator<EdgeHandle> {
public:
    typedef BaseIterator<EdgeHandle> BaseIter;


    EdgeIter(const TopologyKernel* _mesh, const EdgeHandle& _eh = EdgeHandle(0));

    // Post increment/decrement operator
    EdgeIter operator++(int) {
        EdgeIter cpy = *this;
        ++(*this);
        return cpy;
    }
    EdgeIter operator--(int) {
        EdgeIter cpy = *this;
        --(*this);
        return cpy;
    }
    EdgeIter operator+(int _n) {
        EdgeIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    EdgeIter operator-(int _n) {
        EdgeIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    EdgeIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    EdgeIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    EdgeIter& operator++();
    EdgeIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
