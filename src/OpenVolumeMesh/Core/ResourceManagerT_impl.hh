#pragma once
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

#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Core/TypeName.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>

namespace OpenVolumeMesh {


template<typename EntityTag>
detail::Tracker<PropertyStorageBase> &
ResourceManager::storage_tracker() const
{
    return storage_trackers_.get<EntityTag>();
}


template<typename EntityTag>
void ResourceManager::clear_props()
{
    for (auto prop: persistent_props_.get<EntityTag>()) {
        prop->set_persistent(false);
    }
    persistent_props_.get<EntityTag>().clear();
}

template <typename Handle>
void ResourceManager::swap_property_elements(Handle _idx_a, Handle _idx_b)
{
    static_assert(is_handle_v<Handle>);
    for (auto &prop: storage_tracker<typename Handle::EntityTag>()) {
        prop->swap(_idx_a.uidx(), _idx_b.uidx());
    }
}

template <typename Handle>
void ResourceManager::copy_property_elements(Handle _idx_a, Handle _idx_b)
{
    static_assert(is_handle_v<Handle>);
    for (auto &prop: storage_tracker<typename Handle::EntityTag>()) {
        prop->copy(_idx_a.uidx(), _idx_b.uidx());
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

    // TODO: maybe we should only look for persistent props!
    //       also make sure persistent props always have a name!
    for(auto &prop: storage_tracker<EntityTag>())
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
PropertyPtr<T, EntityTag> ResourceManager::internal_create_property(const std::string _name, const T _def) const
{
    auto storage = std::make_shared<PropertyStorageT<T>> (&storage_tracker<EntityTag>(), std::move(_name), EntityTag::type(), _def);
    storage->resize(n<EntityTag>());
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
    return internal_create_property<T, EntityTag>(_name, _def);
}
template<typename T, typename EntityTag>
PropertyPtr<T, EntityTag>
ResourceManager::create_private_property(std::string _name, const T _def) const
{
    return internal_create_property<T, EntityTag>(std::move(_name), _def);
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

    auto sptr = std::static_pointer_cast<PropertyStorageBase>(_prop.storage());
    if (_flag) {
        persistent_props_.get<EntityTag>().insert(sptr);
    } else {
        persistent_props_.get<EntityTag>().erase(sptr);
    }
    sptr->set_persistent(_flag);
}


template<class EntityTag>
void ResourceManager::resize_props(size_t _n)
{
    static_assert(is_entity_v<EntityTag>);
    for (auto &prop: storage_tracker<EntityTag>()) {
        prop->resize(_n);
    }
}

template<class EntityTag>
void ResourceManager::reserve_props(size_t _n)
{
    static_assert(is_entity_v<EntityTag>);
    for (auto &prop: storage_tracker<EntityTag>()) {
        prop->reserve(_n);
    }
}


template<typename  Handle>
void ResourceManager::entity_deleted(Handle _h)
{
    static_assert(is_handle_v<Handle>);
    for (auto &prop: storage_tracker<typename Handle::EntityTag>()) {
        prop->delete_element(_h.uidx());
    }
}


template<typename EntityTag>
size_t ResourceManager::n_props() const {
    return storage_tracker<EntityTag>().size();
}

template<typename EntityTag>
size_t ResourceManager::n_persistent_props() const {
    return persistent_props_.get<EntityTag>().size();
}


template<bool Move, typename EntityTag>
void ResourceManager::assignProperties(typename std::conditional<Move, ResourceManager&, const ResourceManager&>::type src)
{
    auto &src_all = src.persistent_props_.template get<EntityTag>();
    auto &dst_all = persistent_props_.get<EntityTag>();

    // If possible, re-use existing properties instead of copying
    // everything blindly.

    PersistentProperties out;

    // TODO OPT: this will be slow for many props (quadratic) - we could do this in nlogn!
    //           sort both sets by key (name, type), then traverse in parallel
    for (const auto &srcprop: src_all) {
        bool found = false;
        for (auto it = dst_all.begin(); it != dst_all.end(); ++it)
        {
            auto &dstprop = *it;
            if (dstprop->name() == srcprop->name()
                    && dstprop->internal_type_name() == srcprop->internal_type_name())
            {
                // found a correspondence!
                out.insert(dstprop);
                if (Move) {
                    dstprop->move_values_from(srcprop.get());
                } else {
                    dstprop->assign_values_from(srcprop.get());
                }
                dst_all.erase(it);
                found = true;
                break;
            }
        }

        if (!found)
        {
            std::shared_ptr<PropertyStorageBase> prop;
            if (Move) {
                prop = srcprop;
            } else {
                prop = srcprop->clone();
            }
            prop->set_tracker(&storage_tracker<EntityTag>());
            out.insert(std::move(prop));
        }
    }
    dst_all = std::move(out);
}

template<bool Move>
void ResourceManager::assignAllPropertiesFrom(typename std::conditional<Move, ResourceManager&, const ResourceManager&>::type src)
{
    // TODO: we need some helper function to clean up this dispatch:
    for_each_entity([&](auto entity_tag) {
        assignProperties<Move, decltype(entity_tag)>(src);

    });
    assignProperties<Move, Entity::Vertex>(src);
    assignProperties<Move, Entity::Vertex>(src);
    assignProperties<Move, Entity::Edge>(src);
    assignProperties<Move, Entity::HalfEdge>(src);
    assignProperties<Move, Entity::Face>(src);
    assignProperties<Move, Entity::HalfFace>(src);
    assignProperties<Move, Entity::Cell>(src);
    assignProperties<Move, Entity::Mesh>(src);
}



template <class T, typename Entity>
PropertyPtr<T, Entity>::PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def)
{
    *this = mesh->request_property<T, Entity>(_name, _def);
}


} // Namespace OpenVolumeMesh
