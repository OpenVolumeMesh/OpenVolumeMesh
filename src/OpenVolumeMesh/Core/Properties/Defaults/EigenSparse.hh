#pragma once
#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <Eigen/SparseCore>

namespace OpenVolumeMesh {

template<typename Scalar, int Options, typename StorageIndex>
struct default_prop<Eigen::SparseMatrix<Scalar, Options, StorageIndex>>
{
    static inline const Eigen::SparseMatrix<Scalar, Options, StorageIndex> value = {};
};

}
