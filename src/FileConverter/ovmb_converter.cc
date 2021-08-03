#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/IO/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/BinaryFileWriter.hh>

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Geometry/Vector11T.hh>

#include <iostream>
#include <string>

namespace OVM = OpenVolumeMesh;

using Mesh = OVM::GeometryKernel<OVM::Geometry::Vec3d, OVM::TopologyKernel>;

int main(int argc, char**argv) {
    if (argc != 3 && argc != 2) {
        std::cout << "OpenVolumeMesh .ovm <-> .ovmb converter\n"
                  << "Usage: " << argv[0] << " <infile> [outfile]\n"
                     "If only an infile is given, nothing happens after reading. Useful for read speed benchmarking." << std::endl;
        return 1;
    }
    std::string fname_in = argv[1];

    bool in_binary = fname_in.ends_with(".ovmb");
    if (!in_binary && !fname_in.ends_with(".ovm"))
    {
        std::cerr << "Error: Input filename needs to be .ovm or .ovmb." << std::endl;
        return 1;
    }
    Mesh mesh;

    bool topo_check = false;
    bool bottom_up = false;

    std::ifstream stream_in(fname_in);

    if (in_binary) {
        OVM::IO::BinaryFileReader reader(stream_in, mesh);
        reader.enable_topology_check(topo_check);
        reader.enable_bottom_up_incidences(bottom_up);
        auto compat = reader.compatibility();
        if (compat != OVM::IO::ReadCompatibility::Ok) {
            std::cout << "Error: Input file not compatible: "
                      << to_string(compat)
                      << std::endl;
            return 2;
        }
        bool success = reader.read_file();
        if (!success) {
            auto state = reader.state();
            std::cout << "Error: Reading binary file failed, state: "
            << to_string(state)
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

    if (argc == 2)
        return 0;
    std::string fname_out = argv[2];

    bool out_binary = fname_out.ends_with(".ovmb");
    if (!out_binary && !fname_out.ends_with(".ovm"))
    {
        std::cerr << "Error: Output filename needs to be .ovm or .ovmb." << std::endl;
        return 1;
    }

    std::ofstream stream_out(fname_out);
    if (out_binary) {
        OVM::IO::BinaryFileWriter writer(stream_out, mesh); // TODO: make API a free function?
        bool success = writer.write();
        if (!success) {
            std::cout << "Error: Writing binary file failed." << std::endl;
            return 5;
        }
    } else {
        OVM::IO::FileManager write_man;
        write_man.writeStream(stream_out, mesh);
    }

    return 0;
}
