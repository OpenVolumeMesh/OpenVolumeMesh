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

#define RESOURCEMANAGERT_CC

#include "ResourceManager.hh"
#include "PropertyDefines.hh"
#include "TypeName.hh"

namespace OpenVolumeMesh {


template<class T>
VertexPropertyT<T> ResourceManager::request_vertex_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,VertexPropertyT<T>,VertexPropHandle,T>(vertex_props_, _name, n_vertices(), _def);
}

template<class T>
EdgePropertyT<T> ResourceManager::request_edge_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,EdgePropertyT<T>,EdgePropHandle,T>(edge_props_, _name, n_edges(), _def);
}

template<class T>
HalfEdgePropertyT<T> ResourceManager::request_halfedge_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,HalfEdgePropertyT<T>,HalfEdgePropHandle,T>(halfedge_props_, _name, n_edges()*2u, _def);
}

template<class T>
FacePropertyT<T> ResourceManager::request_face_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,FacePropertyT<T>,FacePropHandle,T>(face_props_, _name, n_faces(), _def);
}

template<class T>
HalfFacePropertyT<T> ResourceManager::request_halfface_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,HalfFacePropertyT<T>,HalfFacePropHandle,T>(halfface_props_, _name, n_faces()*2u, _def);
}

template<class T>
CellPropertyT<T> ResourceManager::request_cell_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,CellPropertyT<T>,CellPropHandle,T>(cell_props_, _name, n_cells(), _def);
}

template<class T>
MeshPropertyT<T> ResourceManager::request_mesh_property(const std::string& _name, const T _def) {

    return internal_request_property<std::vector<BaseProperty*>,MeshPropertyT<T>,MeshPropHandle,T>(mesh_props_, _name, 1, _def);
}

template<class StdVecT, class PropT, class HandleT, class T>
PropT ResourceManager::internal_request_property(StdVecT& _vec, const std::string& _name, size_t _size, const T _def)
{
    auto type_name = get_type_name<T>();

    if(!_name.empty()) {
        for(typename StdVecT::iterator it = _vec.begin();
                it != _vec.end(); ++it) {
            if((*it)->name() == _name
                && (*it)->internal_type_name() == type_name)
            {
                return *static_cast<PropT*>(*it);
            }
        }
    }

    auto handle = HandleT::from_unsigned(_vec.size());

    PropT* prop = new PropT(_name, type_name, *this, handle, _def);
    prop->resize(_size);

    // Store property pointer
    _vec.push_back(prop);

    return *prop;
}

// request_property: work around C++ currently not allowing partial specialisation on functions by using structs:

template<class T, typename EntityTag>
struct request_property_impl {
    static PropertyTT<T, EntityTag> _(ResourceManager* /*resman*/, const std::string& /*_name*/, const T /*_def*/);
};

template<class T>
struct request_property_impl<T, Entity::Vertex>{
    static PropertyTT<T, Entity::Vertex> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_vertex_property<T>(_name, _def);
    }
};
template<class T>
struct request_property_impl<T, Entity::Edge>{
    static PropertyTT<T, Entity::Edge> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_edge_property<T>(_name, _def);
    }
};
template<class T>
struct request_property_impl<T, Entity::HalfEdge>{
    static PropertyTT<T, Entity::HalfEdge> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_halfedge_property<T>(_name, _def);
    }
};
template<class T>
struct request_property_impl<T, Entity::Face>{
    static PropertyTT<T, Entity::Face> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_face_property<T>(_name, _def);
    }
};
template<class T>
struct request_property_impl<T, Entity::HalfFace>{
    static PropertyTT<T, Entity::HalfFace> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_halfface_property<T>(_name, _def);
    }
};
template<class T>
struct request_property_impl<T, Entity::Cell>{
    static PropertyTT<T, Entity::Cell> _(ResourceManager *resman, const std::string &_name, const T _def) {
        return resman->request_cell_property<T>(_name, _def);
    }
};

template<typename T, typename EntityTag>
PropertyTT<T, EntityTag> ResourceManager::request_property(const std::string& _name, const T _def)
{
    return request_property_impl<T, EntityTag>::_(this, _name, _def);
}


template<typename T, class EntityTag>
void ResourceManager::set_persistent(PropertyTT<T, EntityTag>& _prop, bool _flag)
{
    if(_flag == _prop->persistent()) return;
    _prop->set_persistent(_flag);
}

template<class StdVecT>
void ResourceManager::remove_property(StdVecT& _vec, size_t _idx) {

    auto prop_ptr = _vec[_idx];
    prop_ptr->setResMan(nullptr);
    delete prop_ptr;
    _vec.erase(_vec.begin() + _idx);
    updatePropHandles(_vec);
}

template<class StdVecT>
void ResourceManager::resize_props(StdVecT& _vec, size_t _n) {

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        (*it)->resize(_n);
    }
}

template<class StdVecT>
void ResourceManager::entity_deleted(StdVecT& _vec, const OpenVolumeMeshHandle& _h) {

    for(typename StdVecT::iterator it = _vec.begin();
            it != _vec.end(); ++it) {
        (*it)->delete_element(_h.uidx());
    }
}

template<class StdVecT>
void ResourceManager::clearVec(StdVecT& _vec) {

    for (auto prop: _vec) {
        prop->setResMan(nullptr);
        delete prop;
    }
    _vec.clear();
}

template<class StdVecT>
void ResourceManager::updatePropHandles(StdVecT &_vec)
{
    size_t n = _vec.size();
    for(size_t i = 0; i < n; ++i) {
        _vec[i]->set_handle(OpenVolumeMeshHandle(static_cast<int>(i)));
    }
}

} // Namespace OpenVolumeMesh
