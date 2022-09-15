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
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>

namespace OpenVolumeMesh {


template<typename EntityTag>
detail::Tracker<PropertyStorageBase> &
ResourceManager::storage_tracker() const
{
    return storage_trackers_.get<EntityTag>();
}


// Change all properties to be private, such that they are
// not visible anymore from the mesh API.
// Crucially, this does not delete props that are in use,
// i.e., they are kept in storage_tracker for index updates.
template<typename EntityTag>
void ResourceManager::clear_props()
{
    for (auto prop: persistent_props_.get<EntityTag>()) {
        prop->set_persistent(false);
    }
    persistent_props_.get<EntityTag>().clear();

    for (auto prop: storage_tracker<EntityTag>()) {
        prop->set_shared(false);
    }
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

    auto type_name = detail::internal_type_name<T>();

    for(auto &prop: storage_tracker<EntityTag>())
    {
        if(prop->shared()
                && prop->name() == _name
                && prop->internal_type_name() == type_name)
        {
            return prop_ptr_from_storage<T, EntityTag>(prop);
        }
    }
    return {};
}

template<class T, class EntityTag>
PropertyPtr<T, EntityTag> ResourceManager::internal_create_property(
        std::string _name, const T _def, bool _shared) const
{
    auto storage = std::make_shared<PropertyStorageT<T>>(
                                     &storage_tracker<EntityTag>(),
                                     std::move(_name),
                                     EntityTag::type(),
                                     std::move(_def),
                                     _shared);
    storage->resize(n<EntityTag>());
    return PropertyPtr<T, EntityTag>(std::move(storage));
}

template<typename T, typename EntityTag>
PropertyPtr<T, EntityTag> ResourceManager::request_property(const std::string& _name, const T _def)
{
    auto prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return *prop;
    bool shared = !_name.empty();
    return internal_create_property<T, EntityTag>(_name, _def, shared);
}

template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::create_persistent_property(std::string _name, const T _def)
{
    auto prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return {};
    auto ptr =  internal_create_property<T, EntityTag>(std::move(_name), _def, true);
    set_persistent(ptr);
    return ptr;
}

template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::create_shared_property(std::string _name, const T _def)
{
    auto prop = internal_find_property<T, EntityTag>(_name);
    if (prop)
        return {};
    return internal_create_property<T, EntityTag>(std::move(_name), _def, true);
}
template<typename T, typename EntityTag>
PropertyPtr<T, EntityTag>
ResourceManager::create_private_property(std::string _name, const T _def) const
{
    return internal_create_property<T, EntityTag>(std::move(_name), _def, false);
}

template<typename T, typename EntityTag>
std::optional<PropertyPtr<T, EntityTag>>
ResourceManager::get_property(const std::string& _name)
{
    return internal_find_property<T, EntityTag>(_name);
}

template<typename T, typename EntityTag>
std::optional<const PropertyPtr<T, EntityTag>>
ResourceManager::get_property(const std::string& _name) const
{
    return internal_find_property<T, EntityTag>(_name);
}


template<typename T, class EntityTag>
void ResourceManager::set_shared(PropertyPtr<T, EntityTag>& _prop, bool _enable)
{
    if(_enable == _prop.shared()) return;

    auto sptr = std::static_pointer_cast<PropertyStorageBase>(_prop.storage());
    if (_enable) {
        if (_prop.anonymous()) {
            throw std::runtime_error("Shared properties must have a name!");
        }
        auto existing = internal_find_property<T, EntityTag>(_prop.name());
        if (existing) {
            throw std::runtime_error("A shared property with this name, type and entity type already exists.");
        }
    } else {
        set_persistent(_prop, false);
    }
    sptr->set_shared(_enable);
}

template<typename T, class EntityTag>
void ResourceManager::set_persistent(PropertyPtr<T, EntityTag>& _prop, bool _enable)
{
    if(_enable == _prop.persistent()) return;

    auto sptr = std::static_pointer_cast<PropertyStorageBase>(_prop.storage());
    if (_enable) {
        if (!_prop.shared()) {
            throw std::runtime_error("Persistent properties must be shared (set_shared).");
        }
        persistent_props_.get<EntityTag>().insert(sptr);
    } else {
        persistent_props_.get<EntityTag>().erase(sptr);
    }
    sptr->set_persistent(_enable);
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


#if 0
template<typename EntityTag>
void ResourceManager::assignPropertiesFrom(ResourceManager const& _src)
{
    auto &src_all = _src.template storage_tracker<EntityTag>();
    auto &dst_all = storage_tracker<EntityTag>();

    // If possible, re-use existing properties instead of copying
    // everything blindly.

    auto &persistent = persistent_props_.get<EntityTag>();

    // TODO OPT: this will be slow for many props (quadratic) - we could do this in nlogn.
    //           sort both sets by key (name, type), then traverse in parallel
    for (const auto &srcprop: src_all) {
        bool found = false;
        if (!srcprop->shared()) {
            // it does not make sense to copy private props,
            // noone could access them
            return;
        }
        // try to find and update existing properties in dst
        for (auto it = dst_all.begin(); it != dst_all.end(); ++it)
        {
            auto &dstprop = *it;
            if (!dstprop->shared()) {
                return;
            }
            if (dstprop->name() == srcprop->name()
                    && dstprop->internal_type_name() == srcprop->internal_type_name())
            {
                // found a correspondence!
                dstprop->assign_values_from(srcprop);
                if (srcprop->persistent() && !dstprop->persistent()) {
                    persistent.insert(dstprop->shared_from_this());
                }
                found = true;
                break;
            }
        }

        if (!found)
        {
            auto dstprop = srcprop->clone();
            dstprop->set_tracker(&storage_tracker<EntityTag>());
            if (srcprop->persistent()) {
                persistent.insert(dstprop->shared_from_this());
            }
        }
    }
}
#endif



// We define this here to avoid circular dependencies:

template <class T, typename Entity>
PropertyPtr<T, Entity>::PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def)
{
    *this = mesh->request_property<T, Entity>(_name, _def);
}


template <class T, typename EntityTag>
PropertyPtr<T, EntityTag> ResourceManager::
prop_ptr_from_storage(PropertyStorageBase *_prop)
{
    auto ps = std::static_pointer_cast<PropertyStorageT<T>>(
            _prop->shared_from_this());
    return PropertyPtr<T, EntityTag>(std::move(ps));
}



} // Namespace OpenVolumeMesh
