#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT VertexIter : public BaseIterator<VertexHandle> {
public:
    typedef BaseIterator<VertexHandle> BaseIter;


    VertexIter(const TopologyKernel* _mesh, const VertexHandle& _vh = VertexHandle(0));

    // Post increment/decrement operator
    VertexIter operator++(int) {
        VertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    VertexIter operator--(int) {
        VertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    VertexIter operator+(int _n) {
        VertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    VertexIter operator-(int _n) {
        VertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    VertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    VertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    VertexIter& operator++();
    VertexIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
