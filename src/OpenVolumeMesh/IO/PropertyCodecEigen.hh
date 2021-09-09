#pragma once

#include <OpenVolumeMesh/IO/PropertyCodec.hh>
#include <OpenVolumeMesh/IO/PropertyCodecT_impl.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>

#include <Eigen/Core>


namespace OpenVolumeMesh::IO::Codecs {

template<typename _Scalar, int _Rows, int _Cols>
struct EigenDenseFixedMatrix
{
    static_assert(_Rows > 0);
    static_assert(_Cols > 0);
    using T = Eigen::Matrix<_Scalar, _Rows, _Cols>;

    static void encode(detail::Encoder &enc, const T &val) {
        for (const auto x: val.reshaped()) {
            enc.write(x);
        }
    }
    static void decode(detail::Decoder &reader, T &val) {
        for (const auto x: val.reshaped()) {
            reader.read(x);
        }
    }

};

} // namespace OpenVolumeMesh::IO::Codecs

namespace OpenVolumeMesh::IO {

void register_eigen_codecs(PropertyCodecs &_codecs)
{
    using namespace Codecs;
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 2, 1>>>("2d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 3, 1>>>("3d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 4, 1>>>("4d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 9, 1>>>("9d");

    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 2, 1>>>("2f");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 3, 1>>>("3f");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 4, 1>>>("4f");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 9, 1>>>("9f");

    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 2, 2>>>("2x2d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 3, 3>>>("3x3d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 4, 4>>>("4x4d");

    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 2, 2>>>("2x2d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 3, 3>>>("3x3d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 4, 4>>>("4x4d");

    // TODO: dynamic (dense/sparse) types
}

} // namespace OpenVolumeMesh::IO
