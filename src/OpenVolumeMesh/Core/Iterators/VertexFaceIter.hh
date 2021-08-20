#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT VertexFaceIter : public BaseCirculator<
  VertexHandle,
  FaceHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      FaceHandle> BaseIter;

    VertexFaceIter(const VertexHandle& _vIdx, const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexFaceIter operator++(int) {
    VertexFaceIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexFaceIter operator--(int) {
    VertexFaceIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexFaceIter operator+(int _n) {
    VertexFaceIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexFaceIter operator-(int _n) {
    VertexFaceIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexFaceIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexFaceIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexFaceIter& operator++();
  VertexFaceIter& operator--();

private:
    std::vector<FaceHandle> faces_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
