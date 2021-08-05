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

#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Core/TypeName.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>

namespace OpenVolumeMesh {


template<typename EntityTag>
void ResourceManager::clear_props()
{
    for (auto &prop: all_props_.get<EntityTag>()) {
        prop->setResMan(nullptr);
    }
    all_props_.get<EntityTag>().clear();
    persistent_props_.get<EntityTag>().clear();
}

template <typename EntityTag>
void ResourceManager::swap_property_elements(HandleT<EntityTag> const &_idx_a, HandleT<EntityTag> &_idx_b)
{
    for (auto &prop: all_props_.get<EntityTag>()) {
        prop->swap(_idx_a.uidx(), _idx_b.uidx());
    }
}

template<class T>
VertexPropertyT<T> ResourceManager::request_vertex_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::Vertex>(_name, _def);
}

template<class T>
EdgePropertyT<T> ResourceManager::request_edge_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::Edge>(_name, _def);
}

template<class T>
HalfEdgePropertyT<T> ResourceManager::request_halfedge_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::HalfEdge>(_name, _def);
}

template<class T>
FacePropertyT<T> ResourceManager::request_face_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::Face>(_name, _def);
}

template<class T>
HalfFacePropertyT<T> ResourceManager::request_halfface_property(const std::string& _name, const T _def) {
    return request_property<T, Entity::HalfFace>(_name, _def);
}

template<class T>
CellPropertyT<T> ResourceManager::request_cell_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::Cell>(_name, _def);
}

template<class T>
MeshPropertyT<T> ResourceManager::request_mesh_property(const std::string& _name, const T _def) {

    return request_property<T, Entity::Mesh>(_name, _def);
}

template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::internal_find_property(const std::string& _name) const
{
    if(_name.empty()) {
        return {};
    }

    auto type_name = get_type_name(typeid(T));

    for(auto &prop: all_props_.get<EntityTag>())
    {
        if(prop->name() == _name
            && prop->internal_type_name() == type_name)
        {
            auto ps = std::static_pointer_cast<PropertyStorageT<T>>(
                        prop->shared_from_this());
            return PropertyPtr<T, EntityTag>(std::move(ps));
        }
    }
    return {};
}

template<class T, class EntityTag>
PropertyPtr<T, EntityTag> ResourceManager::internal_create_property(const std::string& _name, const T _def)
{
    auto type_name = get_type_name(typeid(T));
    auto &propVec = all_props_.get<EntityTag>();
    auto storage = std::make_shared<PropertyStorageT<T>>(_name, type_name, EntityTag::type(), _def);
    storage->resize(n<EntityTag>());
    storage->setResMan(this);
    propVec.insert(storage.get());
    return PropertyPtr<T, EntityTag>(std::move(storage));
}

template<typename T, typename EntityTag>
PropertyPtr<T, EntityTag> ResourceManager::request_property(const std::string& _name, const T _def)
{
    auto prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return *prop;
    return internal_create_property<T, EntityTag>(_name, _def);
}


template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::create_property(const std::string& _name, const T _def)
{
    auto *prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return {};
    return {*internal_create_property<T, EntityTag>(_name, _def)};
}

template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::get_property(const std::string& _name)
{
    auto *prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return {*prop};
    return {};
}


template<typename T, class EntityTag>
void ResourceManager::set_persistent(PropertyPtr<T, EntityTag>& _prop, bool _flag)
{
    if(_flag == _prop.persistent()) return;

    auto sptr = std::static_pointer_cast<PropertyStorageBase>(_prop.ptr());
    if (_flag) {
        persistent_props_.get<EntityTag>().insert(sptr);
    } else {
        persistent_props_.get<EntityTag>().erase(sptr);
    }
    sptr->set_persistent(_flag);
}


template<class EntityTag>
void ResourceManager::delete_multiple_entities(const std::vector<bool>& _tags)
{
    for (auto &prop: all_props_.get<EntityTag>()) {
        prop->delete_multiple_entries(_tags);
    }
}

template<class Container>
void ResourceManager::resize_props(Container& _vec, size_t _n)
{
    for (auto &prop: _vec) {
        prop->resize(_n);
    }
}

template<class Container>
void ResourceManager::reserve_props(Container& _vec, size_t _n)
{
    for (auto &prop: _vec) {
        prop->reserve(_n);
    }
}


template<class Container>
void ResourceManager::entity_deleted(Container& _vec, const OpenVolumeMeshHandle& _h) {

    for (auto &prop: _vec) {
        prop->delete_element(_h.uidx());
    }
}


template<typename EntityTag>
size_t ResourceManager::n_props() const {
    auto const &props = all_props_.get<EntityTag>();
    return props.size();
}



template <class T, typename Entity>
PropertyPtr<T, Entity>::PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def)
{
    *this = mesh->request_property<T, Entity>(_name, _def);
}
} // Namespace OpenVolumeMesh
