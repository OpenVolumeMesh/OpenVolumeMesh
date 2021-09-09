#pragma once

namespace OpenVolumeMesh::IO {


struct ReadOptions {
    bool topology_check = true;
    bool bottom_up_incidences = true;
};

} // namespace OpenVolumeMesh::IO
