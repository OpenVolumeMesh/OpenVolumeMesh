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


#define FILEMANAGERT_CC

#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <typeinfo>
#include <stdint.h>
#include <memory>
#include <array>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

namespace OpenVolumeMesh {

namespace IO {

using namespace OpenVolumeMesh::Geometry;

//==================================================

template<class MeshT>
bool FileManager::readStream(std::istream &_istream, MeshT &_mesh,
    bool _topologyCheck, bool _computeBottomUpIncidences) const
{
    std::stringstream sstr;
    std::string line;
    std::string s_tmp;
    typedef typename MeshT::PointT Point;
    Point v = Point(0.0, 0.0, 0.0);

    _istream.imbue(std::locale::classic());

    _mesh.clear(false);
    // Temporarily disable bottom-up incidences
    // since it's way faster to first add all the
    // geometry and compute them in one pass afterwards
    _mesh.enable_bottom_up_incidences(false);

    /*
     * Header
     */

    bool header_found = true;

    // Get first line
    getCleanLine(_istream, line);
    sstr.str(line);

    // Check header
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "OVM") {
        //_istream.close();
        header_found = false;
        if (verbosity_level_ >= 1) {
            std::cerr << "The specified file might not be in OpenVolumeMesh format!" << std::endl;
        }
        //return false;
    }

    // Get ASCII/BINARY string
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp == "BINARY") {
        if (verbosity_level_ >= 1) {
            std::cerr << "Binary files are not supported at the moment!" << std::endl;
        }
        return false;
    }

    /*
     * Vertices
     */
    if(!header_found) {
        sstr.clear();
        sstr.str(line);
    } else {
        getCleanLine(_istream, line);
        sstr.clear();
        sstr.str(line);
    }

    size_t n_vertices = 0;
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "VERTICES") {
        if (verbosity_level_ >= 1) {
            std::cerr << "No vertex section defined!" << std::endl;
        }
        return false;
    } else {

        // Read in number of vertices
        getCleanLine(_istream, line);
        sstr.clear();
        sstr.str(line);
        sstr >> n_vertices;

        _mesh.reserve_vertices(static_cast<size_t>(n_vertices));

        // Read in vertices
        for(uint64_t i = 0u; i < n_vertices; ++i) {

            getCleanLine(_istream, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v[0];
            sstr >> v[1];
            sstr >> v[2];
            _mesh.add_vertex(v);
        }
    }

    /*
     * Edges
     */
    size_t n_edges = 0;
    getCleanLine(_istream, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "EDGES") {
        if (verbosity_level_ >= 2) {
            std::cerr << "No edge section defined!" << std::endl;
        }
        return false;
    } else {

        // Read in number of edges
        getCleanLine(_istream, line);
        sstr.clear();
        sstr.str(line);
        sstr >> n_edges;

        _mesh.reserve_edges(static_cast<size_t>(n_edges));

        // Read in edges
        for(uint64_t i = 0u; i < n_edges; ++i) {

            unsigned int v1 = 0;
            unsigned int v2 = 0;
            getCleanLine(_istream, line);
            sstr.clear();
            sstr.str(line);
            sstr >> v1;
            sstr >> v2;
            if (v1 >= n_vertices || v2 >= n_vertices) {
                std::cerr << "OVM file loading error: invalid vertex for edge #" << i
                          << " = (" << v1 << ", " << v2 << ") - there are only" << n_vertices << " vertices." << std::endl;
                return false;
            }
            // TODO: check index validity
            _mesh.add_edge(VertexHandle(v1), VertexHandle(v2), true);
        }
    }

    size_t n_halfedges = 2 * n_edges;
    /*
     * Faces
     */
    size_t n_faces = 0;
    getCleanLine(_istream, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "FACES") {
        if (verbosity_level_ >= 2) {
            std::cerr << "No face section defined!" << std::endl;
        }
        return false;
    } else {

        // Read in number of faces
        getCleanLine(_istream, line);
        sstr.clear();
        sstr.str(line);
        sstr >> n_faces;

        _mesh.reserve_faces(static_cast<size_t>(n_faces));

        std::vector<HalfEdgeHandle> hes;

        // Read in faces
        for(uint64_t i = 0u; i < n_faces; ++i) {

            getCleanLine(_istream, line);
            sstr.clear();
            sstr.str(line);


            // Get face valence
            uint64_t val = 0u;
            sstr >> val;

            hes.clear();
            hes.reserve(static_cast<size_t>(val));
            // Read half-edge indices
            for(unsigned int e = 0; e < val; ++e) {

                unsigned int v1 = 0;
                sstr >> v1;
                if (v1 >= n_halfedges) {
                    std::cerr << "OVM File loading error: invalid halfedge #" << v1 << " in face #" << i
                              << " - there are only" << n_halfedges << " halfedges." << std::endl;
                    return false;

                }
                // TODO: check index validity
                hes.emplace_back(v1);
            }

            _mesh.add_face(hes, _topologyCheck);
        }
    }
    size_t n_halffaces = 2 * n_faces;

    /*
     * Cells
     */
    size_t n_cells;
    getCleanLine(_istream, line);
    sstr.clear();
    sstr.str(line);
    sstr >> s_tmp;
    std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
    if(s_tmp != "POLYHEDRA") {
        if (verbosity_level_ >= 2) {
            std::cerr << "No polyhedra section defined!" << std::endl;
        }
        return false;
    } else {

        // Read in number of cells
        getCleanLine(_istream, line);
        sstr.clear();
        sstr.str(line);
        sstr >> n_cells;

        _mesh.reserve_cells(static_cast<size_t>(n_cells));

        std::vector<HalfFaceHandle> hfs;
        // Read in cells
        for(uint64_t i = 0u; i < n_cells; ++i) {

            getCleanLine(_istream, line);
            sstr.clear();
            sstr.str(line);


            // Get cell valence
            uint64_t val = 0u;
            sstr >> val;

            hfs.clear();
            hfs.reserve(static_cast<size_t>(val));

            // Read half-face indices
            for(unsigned int f = 0; f < val; ++f) {

                unsigned int v1 = 0;
                sstr >> v1;
                if (v1 >= n_halffaces) {
                    std::cerr << "OVM File loading error: invalid halfface #" << v1 << " in cell #" << i
                              << " - there are only" << n_halffaces << " halffaces." << std::endl;
                    return false;

                }
                // TODO: check index validity
                hfs.emplace_back(v1);
            }

            _mesh.add_cell(hfs, _topologyCheck);
        }
    }

    while(!_istream.eof()) {
        // "End of file reached while searching for input!"
        // is thrown here. \TODO Fix it!

        // Read property
        readProperty(_istream, _mesh);
    }

    if(_computeBottomUpIncidences) {
        // Compute bottom-up incidences
        _mesh.enable_bottom_up_incidences(true);
    }

    if (verbosity_level_ >= 2) {
        std::cerr << "######## openvolumemesh info #########" << std::endl;
        std::cerr << "#vertices: " << _mesh.n_vertices() << std::endl;
        std::cerr << "#edges:    " << _mesh.n_edges() << std::endl;
        std::cerr << "#faces:    " << _mesh.n_faces() << std::endl;
        std::cerr << "#cells:    " << _mesh.n_cells() << std::endl;
        std::cerr << "######################################" << std::endl;
    }

    return true;
}

template <class MeshT>
bool FileManager::readFile(const std::string& _filename, MeshT& _mesh,
    bool _topologyCheck, bool _computeBottomUpIncidences) const {

    std::ifstream iff(_filename.c_str(), std::ios::in);

    auto read_buf = std::make_unique<std::array<char, 0x10000>>();
    iff.rdbuf()->pubsetbuf(read_buf->data(), read_buf->size());

    if(!iff.good()) {
        if (verbosity_level_ >= 1) {
            std::cerr << "Error: Could not open file " << _filename << " for reading!" << std::endl;
        }
        iff.close();
        return false;
    }
    return readStream(iff, _mesh, _topologyCheck,_computeBottomUpIncidences);
}

//==================================================

template <class MeshT>
void FileManager::readProperty(std::istream& _iff, MeshT& _mesh) const {

    std::string line, entity_t, prop_t, name;
    std::stringstream sstr;

    getCleanLine(_iff, line);

    if(line.empty()) return;

    sstr.clear();
    sstr.str(line);
    sstr >> entity_t;
    std::transform(entity_t.begin(), entity_t.end(), entity_t.begin(), ::tolower);
    sstr >> prop_t;
    std::transform(prop_t.begin(), prop_t.end(), prop_t.begin(), ::tolower);
    name = line;
    extractQuotedText(name);

    if (verbosity_level_ >= 2) {
        std::cerr << "OVM read property " << name << " of type " << prop_t << std::endl;
    }

    if(prop_t == typeName<int>()) generateGenericProperty<int, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned int>()) generateGenericProperty<unsigned int, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<short>()) generateGenericProperty<short, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<long>()) generateGenericProperty<long, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned long>()) generateGenericProperty<unsigned long, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<char>()) generateGenericProperty<char, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<unsigned char>()) generateGenericProperty<unsigned char, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<bool>()) generateGenericProperty<bool, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<float>()) generateGenericProperty<float, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<double>()) generateGenericProperty<double, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::string>()) generateGenericProperty<std::string, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::map<HalfEdgeHandle, int> >()) generateGenericProperty<std::map<HalfEdgeHandle, int>, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::vector<double> >()) generateGenericProperty<std::vector<double> , MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::vector<VertexHandle> >()) generateGenericProperty<std::vector<VertexHandle> , MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::vector<HalfFaceHandle> >()) generateGenericProperty<std::vector<HalfFaceHandle> , MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<std::vector<std::vector<HalfFaceHandle> > >()) generateGenericProperty<std::vector<std::vector<HalfFaceHandle> > , MeshT>(entity_t, name, _iff, _mesh);

    else if(prop_t == typeName<Vec2f>()) generateGenericProperty<Vec2f, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec2d>()) generateGenericProperty<Vec2d, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec2i>()) generateGenericProperty<Vec2i, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec2ui>()) generateGenericProperty<Vec2ui, MeshT>(entity_t, name, _iff, _mesh);

    else if(prop_t == typeName<Vec3f>()) generateGenericProperty<Vec3f, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec3d>()) generateGenericProperty<Vec3d, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec3i>()) generateGenericProperty<Vec3i, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec3ui>()) generateGenericProperty<Vec3ui, MeshT>(entity_t, name, _iff, _mesh);

    else if(prop_t == typeName<Vec4f>()) generateGenericProperty<Vec4f, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec4d>()) generateGenericProperty<Vec4d, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec4i>()) generateGenericProperty<Vec4i, MeshT>(entity_t, name, _iff, _mesh);
    else if(prop_t == typeName<Vec4ui>()) generateGenericProperty<Vec4ui, MeshT>(entity_t, name, _iff, _mesh);



}

//==================================================

template <class PropT, class MeshT>
void FileManager::generateGenericProperty(const std::string& _entity_t, const std::string& _name,
                                          std::istream& _iff, MeshT& _mesh) const {

    if(_entity_t == "vprop") {
        VertexPropertyT<PropT> prop = _mesh.template request_vertex_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "eprop") {
        EdgePropertyT<PropT> prop = _mesh.template request_edge_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "heprop") {
        HalfEdgePropertyT<PropT> prop = _mesh.template request_halfedge_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "fprop") {
        FacePropertyT<PropT> prop = _mesh.template request_face_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "hfprop") {
        HalfFacePropertyT<PropT> prop = _mesh.template request_halfface_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "cprop") {
        CellPropertyT<PropT> prop = _mesh.template request_cell_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    } else if(_entity_t == "mprop") {
        MeshPropertyT<PropT> prop = _mesh.template request_mesh_property<PropT>(_name);
        prop.deserialize(_iff);
        _mesh.set_persistent(prop);
    }
}

//==================================================


template<class MeshT>
void FileManager::writeStream(std::ostream &_ostream, const MeshT &_mesh) const
{
    _ostream.imbue(std::locale::classic());
    // Write header
    _ostream << "OVM ASCII" << std::endl;

    uint64_t n_vertices(_mesh.n_vertices());
    _ostream << "Vertices" << std::endl;
    _ostream << n_vertices << std::endl;

    typedef typename MeshT::PointT Point;

    // write vertices
    for(VertexIter v_it = _mesh.v_iter(); v_it; ++v_it) {

        Point v = _mesh.vertex(*v_it);
        _ostream << v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    uint64_t n_edges(_mesh.n_edges());
    _ostream << "Edges" << std::endl;
    _ostream << n_edges << std::endl;

    // write edges
    for(EdgeIter e_it = _mesh.e_iter(); e_it; ++e_it) {

        VertexHandle from_vertex = _mesh.edge(*e_it).from_vertex();
        VertexHandle to_vertex = _mesh.edge(*e_it).to_vertex();
        _ostream << from_vertex << " " << to_vertex << std::endl;
    }

    uint64_t n_faces(_mesh.n_faces());
    _ostream << "Faces" << std::endl;
    _ostream << n_faces << std::endl;

    // write faces
    for(FaceIter f_it = _mesh.f_iter(); f_it; ++f_it) {

        _ostream << static_cast<uint64_t>(_mesh.face(*f_it).halfedges().size()) << " ";

        std::vector<HalfEdgeHandle> halfedges = _mesh.face(*f_it).halfedges();

        for(typename std::vector<HalfEdgeHandle>::const_iterator it = halfedges.begin(); it
                                                                                         != halfedges.end(); ++it) {

            _ostream << it->idx();

            if((it + 1) != halfedges.end())
                _ostream << " ";
        }

        _ostream << std::endl;
    }

    uint64_t n_cells(_mesh.n_cells());
    _ostream << "Polyhedra" << std::endl;
    _ostream << n_cells << std::endl;

    for(CellIter c_it = _mesh.c_iter(); c_it; ++c_it) {

        _ostream << static_cast<uint64_t>(_mesh.cell(*c_it).halffaces().size()) << " ";

        std::vector<HalfFaceHandle> halffaces = _mesh.cell(*c_it).halffaces();

        for(typename std::vector<HalfFaceHandle>::const_iterator it = halffaces.begin(); it
                                                                                         != halffaces.end(); ++it) {

            _ostream << it->idx();

            if((it + 1) != halffaces.end())
                _ostream << " ";
        }

        _ostream << std::endl;
    }

    writeProps<Entity::Vertex>  (_ostream, _mesh);
    writeProps<Entity::Edge>    (_ostream, _mesh);
    writeProps<Entity::HalfEdge>(_ostream, _mesh);
    writeProps<Entity::Face>    (_ostream, _mesh);
    writeProps<Entity::HalfFace>(_ostream, _mesh);
    writeProps<Entity::Cell>    (_ostream, _mesh);
    writeProps<Entity::Mesh>    (_ostream, _mesh);
}

template<class MeshT>
bool FileManager::writeFile(const std::string& _filename, const MeshT& _mesh) const {

    std::ofstream off(_filename.c_str(), std::ios::out);

    if(!off.good()) {
        if (verbosity_level_ >= 1) {
            std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        }
        return false;
    }
    writeStream(off, _mesh);
    return off.good();
}

//==================================================
template<typename EntityTag, class MeshT>
void FileManager::writeProps(std::ostream &_ostream, const MeshT& _mesh) const {
    writeProps(_ostream,
               _mesh.template persistent_props_begin<EntityTag>(),
               _mesh.template persistent_props_end<EntityTag>(),
               entityTypeName<EntityTag>()
               );
}


template<class IteratorT>
void FileManager::writeProps(
        std::ostream& _ostr,
        const IteratorT& _begin, const IteratorT& _end,
        std::string const &_entityType
        ) const {

    // write props
    for(IteratorT p_it = _begin; p_it != _end; ++p_it) {
        auto prop = *p_it;
        if(!prop->persistent()) continue;
        if(prop->anonymous()) {
            if (verbosity_level_ >= 2) {
                std::cerr << "Serialization of anonymous properties is not supported!" << std::endl;
            }
            continue;
        }

        std::string type_name;
        try {
            type_name = prop->typeNameWrapper();
        } catch (std::runtime_error &e) { // type not serializable
            if (verbosity_level_ >= 1) {
                std::cerr << "Failed to save property " << prop->name() << " , skipping: " << e.what() << std::endl;
            }
            continue;
        }
        _ostr << _entityType << " ";
        _ostr << type_name << " ";
        _ostr << "\"" << prop->name() << "\"" << std::endl;

        prop->serialize(_ostr);
    }
}

//==================================================

} // Namespace IO

} // Namespace FileManager
