/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/


/*              A U T H O R               */
/*       Pierre-Alexandre BEAUFORT        */
/* pierre-alexandre.beaufort@inf.unibe.ch */


/* Class to read vtk DataFile Version 2.0. */
/* Exclusively read UNSTRUCTURED_GRID with */
/* CELL_DATA defined for **every** cells   */


#ifndef VTKCOLORREADER_HH_
#define VTKCOLORREADER_HH_

#include <string>
#include <fstream>
#include <istream>
#include <ostream>

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>
#include <stdint.h>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>


namespace OpenVolumeMesh {

  namespace Reader {

    using namespace OpenVolumeMesh::Geometry;
    
    //=======================================
    
    /**
     * \class VtkColorReader
     * \brief Read mesh data from vtk DataFile Version 2.0
     *
     * \todo Heriting from FileManager class would possibly be better
     * but FileManager class hasn't been design to handle this kind of format
     * since it appears it has been specialized for OVM format. 
     * A recommendation would be to rewrite FileManager class in order to be generic.
     */

    class VtkColorReader {

    public:
      
      class CellWithColor {

      public:

	CellWithColor() :
	  type(0),
	  index(-1),
	  color(-1)
	{
	}

	int type;

	int index;

	int color;

	std::vector<uint32_t> nodes;
      }; // class CellWithColor
      
      
      VtkColorReader() {
	
      }
      
      ~VtkColorReader() {
	
      }

      inline void nextLine(std::istream &_ifs, std::string &_string) const {
	while (true) {
	  std::getline(_ifs, _string);
	  if (_string.size() != 0)
	    break;
	  if(_ifs.eof())
	    std::cerr << "End of file reached while searching for input!" << std::endl;
	}
	// Trim Both leading and trailing spaces
	size_t start = _string.find_first_not_of(" \t\r\n");
	size_t end = _string.find_last_not_of(" \t\r\n");

	if((std::string::npos == start) || (std::string::npos == end)) {
	  _string = "";
	} else {
	  _string = _string.substr(start, end - start + 1);
	}
      }// inline void nextLine

      template < typename MeshT>
      inline int addVertex(MeshT &_mesh, CellWithColor & _cell) {
	
	if (_cell.nodes.size() != 1)
	  std::cerr << "Expected at least 1 index to add a vertex but got " << _cell.nodes.size() << std::endl;	
	
	return _cell.nodes[0];
      } // int addVertex

      template < typename MeshT>
      inline int addEdge(MeshT &_mesh, CellWithColor & _cell) {
	
	if (_cell.nodes.size() != 2)
	  std::cerr << "Expected 2 indices to add an edge but got " << _cell.nodes.size() << std::endl;

	std::vector<VertexHandle> edge;
	for (const auto & node:  _cell.nodes)
	  edge.push_back(OpenVolumeMesh::VertexHandle(node));


	HalfEdgeHandle heh;
	heh = _mesh.find_halfedge(edge[0], edge[1]);
	if (!heh.is_valid())
	  {
	    EdgeHandle eh = _mesh.add_edge(edge[0], edge[1]);
	    heh = _mesh.halfedge_handle(eh, 0);
	  }


	
	return heh.idx();
      } // int addEdge

      template < typename MeshT>
      inline int addTriangle(MeshT &_mesh, CellWithColor & _cell) {
	
	if (_cell.nodes.size() != 3)
	  std::cerr << "Expected 3 indices to add a triangle but got " << _cell.nodes.size() << std::endl;

	std::vector<VertexHandle> triangle;
	for (const auto & node:  _cell.nodes)
	  triangle.push_back(OpenVolumeMesh::VertexHandle(node));


	HalfFaceHandle hfh;
	hfh = _mesh.find_halfface(triangle);
	if (!hfh.is_valid())
	  {
	    FaceHandle fh = _mesh.add_face(triangle);
	    hfh = _mesh.halfface_handle(fh, 0);
	  }

	FacePropertyT<int> triangle_colors = _mesh.template request_face_property<int>("face_colors");
	FaceHandle fh = _mesh.face_handle(hfh);
	triangle_colors[fh] = _cell.color;
	_mesh.set_persistent(triangle_colors, true);
	
	return hfh.idx();
      } // int addTriangle

      template < typename MeshT>
      inline int addTetra(MeshT &_mesh, CellWithColor & _cell) {

	if (_cell.nodes.size() != 4)
	  std::cerr << "Expected 4 indices to add a tetra but got " << _cell.nodes.size() << std::endl;

	std::vector< std::vector<VertexHandle> > faces;
	faces.push_back(std::vector<VertexHandle>());
	faces[0].push_back(VertexHandle(_cell.nodes[0]));
	faces[0].push_back(VertexHandle(_cell.nodes[1]));
	faces[0].push_back(VertexHandle(_cell.nodes[2]));
	faces.push_back(std::vector<VertexHandle>());
	faces[1].push_back(VertexHandle(_cell.nodes[0]));
	faces[1].push_back(VertexHandle(_cell.nodes[3]));
	faces[1].push_back(VertexHandle(_cell.nodes[1]));
	faces.push_back(std::vector<VertexHandle>());
	faces[2].push_back(VertexHandle(_cell.nodes[1]));
	faces[2].push_back(VertexHandle(_cell.nodes[3]));
	faces[2].push_back(VertexHandle(_cell.nodes[2]));
	faces.push_back(std::vector<VertexHandle>());
	faces[3].push_back(VertexHandle(_cell.nodes[0]));
	faces[3].push_back(VertexHandle(_cell.nodes[2]));
	faces[3].push_back(VertexHandle(_cell.nodes[3]));

	std::vector<HalfFaceHandle> halffacehandles;
	for (unsigned int i = 0; i < faces.size(); i++)
	  halffacehandles.push_back(_mesh.find_halfface(faces[i]));

	for (unsigned int i = 0; i < halffacehandles.size(); i++)
	  if (!halffacehandles[i].is_valid())
	    {
	      FaceHandle fh = _mesh.add_face(faces[i]);
	      halffacehandles[i] = _mesh.halfface_handle(fh, 0);
	    }

	
	return _mesh.add_cell(halffacehandles).idx();
	
      } // int addTetra 


      template <class MeshT>
      inline bool paintMesh(std::vector<CellWithColor> &_cells, MeshT &_mesh) {

	VertexPropertyT<int> vertex_colors = _mesh.template request_vertex_property<int>("vertex_colors");
	EdgePropertyT<int> edge_colors = _mesh.template request_edge_property<int>("edge_colors");
	FacePropertyT<int> face_colors = _mesh.template request_face_property<int>("face_colors");
	CellPropertyT<int> cell_colors = _mesh.template request_cell_property<int>("cell_colors"); 	
	for (auto & current_cell: _cells) {
	  switch (current_cell.type) {
	  case 1: // VERTEX	    
	    vertex_colors[OpenVolumeMesh::VertexHandle(current_cell.index)] = current_cell.color;
	    break;
	  case 3:  // EDGE
	    edge_colors[_mesh.edge_handle(OpenVolumeMesh::HalfEdgeHandle(current_cell.index))] = current_cell.color;
	    break;
	  case 5: // TRIANGLE
	    face_colors[_mesh.face_handle(OpenVolumeMesh::HalfFaceHandle(current_cell.index))] = current_cell.color;
	    break;
	  case 10: // TETRAHEDRON
	    cell_colors[OpenVolumeMesh::CellHandle(current_cell.index)] = current_cell.color;
	    break;
	  default:
	    std::cerr << "Trying to add an unsuppoerted type of cell: " << current_cell.type << std::endl;
	    return false;
	  }
	}
	_mesh.set_persistent(vertex_colors, true);
	_mesh.set_persistent(edge_colors, true);
	_mesh.set_persistent(face_colors, true);
	_mesh.set_persistent(cell_colors, true);
	
	return true;
	
      }// bool paintMesh

      
      template <class MeshT>
      inline bool cells2mesh(std::vector<CellWithColor> &_cells, MeshT &_mesh) {
	_mesh.enable_bottom_up_incidences(true);
	for (auto & current_cell: _cells) {
	  switch (current_cell.type) {
	  case 1: // VERTEX
	    current_cell.index = addVertex(_mesh, current_cell);
	    break;
	  case 3:  // EDGE
	    current_cell.index = addEdge(_mesh, current_cell);
	    break;
	  case 5: // TRIANGLE
	    current_cell.index = addTriangle(_mesh, current_cell);
	    break;
	  case 10: // TETRAHEDRON
	    current_cell.index = addTetra(_mesh, current_cell);
	    break;
	  default:
	    return false;
	  }
	}
	return true;
	
      }// bool cells2mesh

      // QUESTION: is it really usefull to initiliaze the properties?
      template < typename MeshT>
      void initializeColor(MeshT &_mesh) {

	VertexPropertyT<int> vertex_colors = _mesh.template request_vertex_property<int>("vertex_colors", 0);
	//for (auto & vertex: _mesh.vertices())
	//  vertex_colors[vertex] = 0;

	EdgePropertyT<int> edge_colors = _mesh.template request_edge_property<int>("edge_colors", 0);
	//for (auto & edge: _mesh.edges())
	//  edge_colors[edge] = 0;

	FacePropertyT<int> face_colors = _mesh.template request_face_property<int>("face_colors", 0);
	//for (auto & face: _mesh.faces())
	//  face_colors[face] = 0;

	CellPropertyT<int> cell_colors = _mesh.template request_cell_property<int>("cell_colors", 0);
	//for (auto & cell: _mesh.cells())
	//  cell_colors[cell] = -1;

	_mesh.set_persistent(vertex_colors, true);
	_mesh.set_persistent(edge_colors, true);
	_mesh.set_persistent(face_colors, true);
	_mesh.set_persistent(cell_colors, true);

      } // void initializecolor
      
      /**
       * \brief Read vtk DataFile Version 2.0 ASCCI from an std::istream
       *
       * Returns true if the file was successfully read. The mesh is stored in parameter _mesh. If something goes wrong, false is returned.
       *
       * @param _istream        The stream to read the mesh from
       * @param _mesh           A reference to an OpenVolumeMesh instance
       * @param _topologyCheck  Pass true if you want to perform a topology check
       *                        each time an entity is added (slower performance)
       * @param _computeBottomUpIncidences Pass true if you want the file manager
       *                                    to directly compute the bottom-up incidences
       *                                    for the mesh. (Note: These are needed for
       *                                    some iterators to work, see documentation)
       */	      
      template <class MeshT>
      bool readStream(std::istream &_istream, MeshT &_mesh, bool _topologyCheck = true, bool _computeBottomUpIncidences = true)
      {

	std::stringstream sstr;
	std::string current_line;
	std::string s_tmp;
	uint64_t c = 0u;	
	typedef typename MeshT::PointT Point;
	Point pt = Point(0.0, 0.0, 0.0);	
	std::vector<CellWithColor> cells;	
	  
	_mesh.clear(false);
	// Temporarily disable bottom-up incidences
	// since it's way faster to first add all the
	// geometry and compute them in one pass afterwards
	_mesh.enable_bottom_up_incidences(false);
    
	/*    H E A D E R   */
	nextLine(_istream, current_line);
	if (current_line != "# vtk DataFile Version 2.0") {
	  std::cerr << "The mesh file does not correspond to _vtk DataFile Version 2.0_!" << std::endl;
	  return false;
	}

	/*    A S C I I   */
	std::getline(_istream, current_line); // info line that is ignored	
	nextLine(_istream, current_line);
	if (current_line != "ASCII") {
	  std::cerr << "Only ASCII files are supported!" << std::endl;
	  return false;
	}

	/*    U N S T R U C T U R E D _ G R I D   */
	nextLine(_istream, current_line);
	sstr.str(current_line);
	sstr >> s_tmp;
	if (s_tmp != "DATASET") {
	  std::cerr << "Expected _DATASET_ section at this point!" << std::endl;
	  return false;
	}
	sstr >> s_tmp;
	if (s_tmp != "UNSTRUCTURED_GRID") {
	  std::cerr << "Only unstructured grid are supported!" << std::endl;
	  return false;
	}
	sstr.clear();
	
	/*    P O I N T S    */
	nextLine(_istream, current_line);
	sstr.str(current_line);
	sstr >> s_tmp;
	if (s_tmp != "POINTS") {
	  std::cerr << "Expected _POINTS_ section at this point!" << std::endl;
	  return false;
	}
	sstr >> c; // number of points
	for (uint64_t i = 0u; i < c; i++) {
	  sstr.clear();
	  nextLine(_istream, current_line);
	  sstr.str(current_line);
	  for (auto j = 0; j<3; j++)
	    sstr >> pt[j];
	  _mesh.add_vertex(pt);
	}
	sstr.clear();
	
	/*    C E L L S    */
	nextLine(_istream, current_line);
	sstr.str(current_line);
	sstr >> s_tmp;
	if (s_tmp != "CELLS") {
	  std::cerr << "Expected _CELLS_ section at this point!" << std::endl;
	  return false;
	}
	sstr >> c; // number of cells
	for (uint64_t i=0u; i < c; i++) {
	  sstr.clear();
	  CellWithColor current_cell;
	  nextLine( _istream, current_line);
	  sstr.str(current_line);
	  uint16_t npt = 0u;
	  sstr >> npt;
	  for (uint16_t j = 0u; j < npt; j++) {
	    uint32_t tmp = 0u;
	    sstr >> tmp;
	    current_cell.nodes.push_back(tmp);
	  }
	  cells.push_back(current_cell);
	}
	sstr.clear();
	
	/*    C E L L _ T Y P E S    */
	nextLine(_istream, current_line);
	sstr.str(current_line);
	sstr >> s_tmp;
	if (s_tmp != "CELL_TYPES") {
	  std::cerr << "Expected _CELL_TYPES_ section at this point!" << std::endl;
	  return false;
	}
	sstr >> c; // number of cell_types: it should be equal to the number of cells
	for (auto & current_cell: cells) {
	  sstr.clear();
	  nextLine( _istream, current_line);
	  sstr.str(current_line);
	  sstr >> current_cell.type;
	}
	sstr.clear();

	/*    C E L L _ D A T A    */
	nextLine(_istream, current_line);
	sstr.str(current_line);
	sstr >> s_tmp;
	if (s_tmp != "CELL_DATA") {
	  std::cerr << "Expected _CELL_DATA_ section at this point!" << std::endl;
	  return false;
	}
	sstr >> c; // number of cell_data: it should be equal to the number of cells
	nextLine(_istream, current_line); // SCALARS subsection
	nextLine(_istream, current_line); // LOOKUP_TABLE subsection
	for (auto & current_cell: cells) {
	  sstr.clear();
	  nextLine( _istream, current_line);
	  sstr.str(current_line);
	  sstr >> current_cell.color;
	}


	std::cout << "Data are about to be exported" << std::endl;
	
	bool ok = cells2mesh(cells, _mesh);
	if (!ok) {
	  std::cerr << "Error while creating the mesh structure!" << std::endl;
	  return false;
	}
	std::cout << "Colors are about to be initialized" << std::endl;
	
	initializeColor(_mesh);

	std::cout << "Colors are about to be exported" << std::endl;
	
	return paintMesh(cells, _mesh);

      } // bool readDataFile

      
      template <class MeshT>
      bool readFile(const std::string& _filename, MeshT& _mesh,
				 bool _topologyCheck=true, bool _computeBottomUpIncidences=true) {

	std::ifstream iff(_filename.c_str(), std::ios::in);
		
	if(!iff.good()) {
	  std::cerr << "Error: Could not open file " << _filename << " for reading!" << std::endl;        
	  iff.close();
	  return false;
	}
	return readStream(iff, _mesh);
      } // bool readFile
      
    }; // class VtkColorReader
    
  } // namespace Reader

} // namespace OpenVolumemesh


#endif // #ifndef VTKCOLORREADER_HH_
