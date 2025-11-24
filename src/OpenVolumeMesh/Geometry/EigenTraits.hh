#pragma once
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <Eigen/Core>

namespace OpenVolumeMesh::Geometry
{

template<typename Scalar, int DIM>
struct vector_dim<Eigen::Vector<Scalar, DIM>>
{ constexpr static int  dim = DIM;};

}
