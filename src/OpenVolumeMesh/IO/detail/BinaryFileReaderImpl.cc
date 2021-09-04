#include <OpenVolumeMesh/IO/detail/BinaryFileReaderImpl.hh>
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

void BinaryFileReaderImpl::read_header()
{
    if (state_ == ReadState::Init) {
        //stream_.seek(0);
        auto decoder = stream_.make_decoder(ovmb_size<FileHeader>);
        auto valid_magic = read(decoder, file_header_);
        if (valid_magic) {
            state_ = ReadState::HeaderRead;
        } else {
            state_ = ReadState::ErrorInvalidMagic;
        }
    }
}


bool BinaryFileReaderImpl::validate_span(uint64_t total, uint64_t read, ArraySpan const&span)
{
    if (span.base != read) {
        state_ = ReadState::Error;
        return false;
    }
    if (total - read < span.count) {
        state_ = ReadState::Error;
        return false;
    }
    return true;
}

std::vector<uint32_t> BinaryFileReaderImpl::read_valences(Decoder &reader, IntEncoding enc, size_t count)
{
    if (!is_valid(enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return {};
    }

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto read_all = [&](auto read_one)
    {
        if (reader.remaining_bytes() < count * elem_size(enc)) {
            state_ = ReadState::Error;
            return;
        }
        for (size_t i = 0; i < count; ++i) {
            valences.push_back(read_one(reader));
        }
    };

    call_with_decoder(enc, read_all);

    return valences;
}

void BinaryFileReaderImpl::readPropChunk(Decoder &reader)
{
    PropChunkHeader header;
    read(reader, header);
    if (header.idx >= props_.size()) {
        state_ = ReadState::Error; // TODO more specific error
        return;
    }
    auto &prop = props_[header.idx];
    if (!prop.decoder) {
        reader.skip();
        return;
    }
    uint64_t n = 0;
    switch (prop.entity) {
        case PropertyEntity::Vertex:   n = n_verts_read_; break;
        case PropertyEntity::Edge:     n = n_edges_read_; break;
        case PropertyEntity::Face:     n = n_faces_read_; break;
        case PropertyEntity::Cell:     n = n_cells_read_; break;
        case PropertyEntity::HalfEdge: n = 2 * n_edges_read_; break;
        case PropertyEntity::HalfFace: n = 2 * n_faces_read_; break;
        case PropertyEntity::Mesh:     n = 1; break;
    }
    if (header.span.base >= n || n - header.span.base < header.span.count) {
        state_ = ReadState::ErrorHandleRange;
        return;
    }
    prop.decoder->deserialize(prop.prop.get(),
            reader,
            static_cast<size_t>(header.span.base),
            static_cast<size_t>(header.span.base + header.span.count));
}




} // namespace OpenVolumeMesh::IO::detail
