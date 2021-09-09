#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Geometry/Vector11T.hh>

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>

namespace OVM = OpenVolumeMesh;

using MeshT = OVM::GeometryKernel<OVM::Geometry::Vec3d, OVM::TopologyKernel>;

static bool ends_with(std::string &s, std::string needle)
{
    if (s.size() < needle.size())
        return false;
    return s.substr(s.size()-needle.size()) == needle;
}

int main(int argc, char**argv) {
    if (argc != 3 && argc != 2) {
        std::cout << "OpenVolumeMesh .ovm <-> .ovmb converter\n"
                  << "Usage: " << argv[0] << " <infile> [outfile]\n"
                     "If only an infile is given, nothing happens after reading. Useful for read speed benchmarking." << std::endl;
        return 1;
    }
    std::string fname_in = argv[1];

    bool in_binary = ends_with(fname_in, ".ovmb");
    if (!in_binary && !ends_with(fname_in, ".ovm"))
    {
        std::cerr << "Error: Input filename needs to be .ovm or .ovmb." << std::endl;
        return 1;
    }
    MeshT mesh;

    bool topo_check = false;
    bool bottom_up = false;

    std::ifstream stream_in(fname_in, std::ios::binary);

    if (!stream_in.good()) {
        std::cerr << "Error: Could not open input file." << std::endl;
        return 7;
    }

    stream_in.seekg(0, stream_in.end);
    uint64_t input_size = stream_in.tellg();
    stream_in.seekg(0, stream_in.beg);
    auto time_start = std::chrono::steady_clock::now();
    if (in_binary) {
        auto result = OVM::IO::ovmb_read(stream_in, mesh);
        if (result != OVM::IO::ReadResult::Ok) {
            std::cout << "Error: Reading binary file failed, state: "
            << to_string(result)
            << std::endl;
            return 3;
        }
    } else {
        OVM::IO::FileManager read_man;
        bool success = read_man.readStream(stream_in, mesh, topo_check, bottom_up);
        if (!success) {
            std::cout << "Error: Reading ascii file failed." << std::endl;
            return 4;
        }
    }
    std::cout << "Read file with " << mesh.n_vertices() << " vertices." << std::endl;

    auto time_end = std::chrono::steady_clock::now();
    std::cout << "Reading took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count()
              << "ms." << std::endl;

    if (argc == 2)
        return 0;
    std::string fname_out = argv[2];

    bool out_binary = ends_with(fname_out, ".ovmb");
    if (!out_binary && !ends_with(fname_out, ".ovm"))
    {
        std::cerr << "Error: Output filename needs to be .ovm or .ovmb." << std::endl;
        return 1;
    }

    std::ofstream stream_out(fname_out, std::ios::binary);

    if (!stream_out.good()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 7;
    }

    time_start = std::chrono::steady_clock::now();

    if (out_binary) {
        auto result = OVM::IO::ovmb_write(stream_out, mesh);
        if (result != OVM::IO::WriteResult::Ok) {
            std::cerr << "Error: Writing binary file failed: " << to_string(result) << std::endl;
            return 5;
        }
    } else {
        OVM::IO::FileManager write_man;
        write_man.writeStream(stream_out, mesh);
    }
    if (!stream_out.good()) {
        std::cerr << "Error: Writing file failed, stream not good." << std::endl;
        return 6;
    }
    uint64_t output_size = stream_out.tellp();
    time_end = std::chrono::steady_clock::now();
    std::cout << "Writing took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count()
              << "ms." << std::endl;
    std::cout << "Input size:  " << input_size << " bytes,\n"
              << "output size: " << output_size << " bytes"
              << std::setprecision(3)
              << " (" << double(100*output_size)/input_size << "% of input size)"
              << std::endl;


    return 0;
}
