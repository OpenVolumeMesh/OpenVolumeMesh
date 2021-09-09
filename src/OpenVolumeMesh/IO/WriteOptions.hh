#pragma once

namespace OpenVolumeMesh::IO {

struct WriteOptions {
    /// try to determine if a polyhedral mesh can in fact be stored as tet or hex mesh
    //bool detect_specialized_mesh = true;

    enum class TopologyType {
        AutoDetect,
        Polyhedral,
        Tetrahedral,
        Hexahedral,
    } topology_type = TopologyType::AutoDetect;
};

} // namespace OpenVolumeMesh::IO
