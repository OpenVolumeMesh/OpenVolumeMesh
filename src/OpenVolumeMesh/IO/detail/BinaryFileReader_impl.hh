#pragma once

#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/Core/Handles.hh>

#include <istream>
#include <cassert>
#include <numeric>
#include <iostream>
#include <type_traits>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

template<typename MeshT>
ReadResult BinaryFileReader::read_file(MeshT &out)
{
    try {
        static_assert(std::is_base_of_v<TopologyKernel, MeshT>);
        if (compatibility<MeshT>() != ReadCompatibility::Ok) {
            state_ = ReadState::ErrorIncompatible;
            return ReadResult::IncompatibleMesh;
        }
        if (state_ != ReadState::HeaderRead) {
            return ReadResult::InvalidFile;
        }
        using GeomReader = detail::GeometryReaderT<typename MeshT::Point>;
        geometry_reader_ = std::make_unique<GeomReader>(out.vertex_positions());
        return internal_read_file(out);
    } catch (parse_error &e) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = std::string("parse_error: ") + e.what();
        return ReadResult::InvalidFile;
    } catch (std::exception &e) {
        state_ = ReadState::Error;
        error_msg_ = std::string("exception: ") + e.what();
        return ReadResult::OtherError;
    }
}

constexpr const inline auto max_handle_idx = static_cast<size_t>(std::numeric_limits<int>::max());

template<typename MeshT>
ReadCompatibility BinaryFileReader::compatibility()
{
    static_assert(std::is_base_of_v<TopologyKernel, MeshT>);
    read_header();
    if (state_ == ReadState::Error) { // TODO: check for more generic error state!
        error_msg_ = "Error reading header.";
        return ReadCompatibility::InvalidFile;
    }
    if (file_header_.header_version != 1) {
        error_msg_ = "Header version unsupported";
        return ReadCompatibility::FileVersionUnsupported;
    }
    if (file_header_.file_version != 1) {
#ifndef NDEBUG
        std::cerr << "OVM::IO: file version is too new, still trying to read." << std::endl;
#endif
    }
    // TODO: allow converting vertex dimension?
    if (file_header_.vertex_dim != MeshT::Point::dim()) {
        return ReadCompatibility::MeshVertexDimensionIncompatible;
    }
    if (std::is_base_of<TetrahedralMeshTopologyKernel, MeshT>::value) {
        if (file_header_.topo_type != TopoType::Tetrahedral) {
            return ReadCompatibility::MeshTopologyIncompatible;
        }
    }
    if (std::is_base_of<HexahedralMeshTopologyKernel, MeshT>::value) {
        if (file_header_.topo_type != TopoType::Hexahedral) {
            return ReadCompatibility::MeshTopologyIncompatible;
        }
    }
    uint64_t max_index = std::max(
                std::max(file_header_.n_verts, file_header_.n_edges),
                std::max(file_header_.n_faces, file_header_.n_cells));

    if (max_index > max_handle_idx)
    {
        return ReadCompatibility::MeshHandleIncompatible;
    }
    return ReadCompatibility::Ok;
}




} // namespace OpenVolumeMesh::IO::detail
