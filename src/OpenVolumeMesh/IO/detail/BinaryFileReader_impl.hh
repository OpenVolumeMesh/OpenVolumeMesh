#pragma once

#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>
#include <OpenVolumeMesh/Core/Handles.hh>

#include <istream>
#include <cassert>
#include <numeric>
#include <iostream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

template<typename MeshT>
ReadResult BinaryFileReader::read_file(MeshT &out)
{
    if (compatibility<MeshT>() != ReadCompatibility::Ok) {
        state_ = ReadState::ErrorIncompatible;
        return ReadResult::IncompatibleMesh;
    }
    if (state_ != ReadState::HeaderRead) {
        return ReadResult::InvalidFile;
    }
    return internal_read_file(out);
}

constexpr const inline auto max_handle_idx = static_cast<size_t>(std::numeric_limits<int>::max());

template<typename MeshT>
ReadCompatibility BinaryFileReader::compatibility()
{
    read_header();
    if (state_ == ReadState::Error) { // TODO: check for more generic error state!
        return ReadCompatibility::InvalidFile;
    }
    if (file_header_.header_version > 0) {
        return ReadCompatibility::FileVersionUnsupported;
    }
    if (file_header_.file_version > 0) {
        std::cerr << "OVM::IO: file version is too new, still trying to read." << std::endl;
    }
    // TODO: allow converting vertex dimension?
    if (file_header_.vertex_dim != geometry_reader_->dim()) {
        return ReadCompatibility::MeshVertexDimensionIncompatible;
    }
    // TODO: check geometry type (float, double)

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
