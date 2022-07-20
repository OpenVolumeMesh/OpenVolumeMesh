#pragma once

#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/PropertyCodecsT_impl.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>

#include <Eigen/Core>


namespace OpenVolumeMesh::IO::Codecs {

template<typename _Scalar, int _Rows, int _Cols>
struct EigenDenseFixedMatrix
{
    // TODO: replace with register_matrixlike
    static_assert(_Rows > 0);
    static_assert(_Cols > 0);
    using T = Eigen::Matrix<_Scalar, _Rows, _Cols>;

    static void encode(detail::Encoder &enc, const T &val) {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                enc.write(val(r, c));
            }
        }
    }

    static void decode(detail::Decoder &reader, T &val) {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                reader.read(val(r, c));
            }
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

    // matrices stored as column-major:

    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 2, 2>>>("2x2d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 3, 3>>>("3x3d");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<double, 4, 4>>>("4x4d");

    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 2, 2>>>("2x2f");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 3, 3>>>("3x3f");
    _codecs.register_codec<SimplePropCodec<EigenDenseFixedMatrix<float, 4, 4>>>("4x4f");

    // TODO: dynamic (dense/sparse) types
}

} // namespace OpenVolumeMesh::IO
