#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif


namespace OpenVolumeMesh {

class OVM_EXPORT CellIter : public BaseIterator<CellHandle> {
public:
    typedef BaseIterator<CellHandle> BaseIter;


    CellIter(const TopologyKernel* _mesh, const CellHandle& _ch = CellHandle(0));

    // Post increment/decrement operator
    CellIter operator++(int) {
        CellIter cpy = *this;
        ++(*this);
        return cpy;
    }
    CellIter operator--(int) {
        CellIter cpy = *this;
        --(*this);
        return cpy;
    }
    CellIter operator+(int _n) {
        CellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    CellIter operator-(int _n) {
        CellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    CellIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    CellIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    CellIter& operator++();
    CellIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
