// C++ includes
#include <iostream>
#include <vector>

// Include VtkColorReader class
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>

// Include vector classes
#include <OpenVolumeMesh/Geometry/VectorT.hh>

// Include polyhedral mesh kernel
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

// Make some typedefs to facilitate your life
typedef OpenVolumeMesh::Geometry::Vec3f         Vec3f;
typedef OpenVolumeMesh::GeometryKernel<Vec3f>   PolyhedralMeshV3f;


template<class IteratorT>
void stateProperties(const IteratorT& _begin, const IteratorT& _end)  {
  for (IteratorT p_it = _begin; p_it != _end; ++p_it) {
    auto _it = *p_it;
    if (!_it->persistent()) {
      std::cout << "N O T \t P E R S I S T E N T \t ? !" << std::endl;
      continue;
    }
    if (_it->anonymous()) std::cout << "A N O N Y M O U S \t ? !" << std::endl;
    std::cout << _it->entityType() << " " << _it->typeNameWrapper() << " "  << _it->name() << std::endl;
    _it->serialize(std::cout);
  }

}

int main (int argc, char* argv[]) {

  OpenVolumeMesh::Reader::VtkColorReader vtk;
  
  OpenVolumeMesh::GeometricPolyhedralMeshV3d mesh;

  vtk.readFile("../examples/vtk_datafile_2/vtk_files/s01c_cube.vtk", mesh, true, true);

  std::cout << mesh.n_vertices() << " vertices" << std::endl;
  std::cout << mesh.n_cells() << " tetrahedra" << std::endl;

  //stateProperties(mesh.vertex_props_begin(), mesh.vertex_props_end());
  //stateProperties(mesh.edge_props_begin(), mesh.edge_props_end());
  //stateProperties(mesh.face_props_begin(), mesh.face_props_end());
  //stateProperties(mesh.cell_props_begin(), mesh.cell_props_end());
  
  return 0;
}
