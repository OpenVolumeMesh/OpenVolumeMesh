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

namespace OpenVolumeMesh {

ResourceManager::ResourceManager(const ResourceManager &other)
{
  clone_persistent_properties_from(other);
}


/**
 * This turns all existing properties on the mesh into *private* properties,
 * i.e., they can not be found via the mesh API anymore.
 * However, properties in use (i.e. with a PropertyPtr)
 * are resized to fit the new n_"entity" value and continue to work.
 * Be aware that their contents might be semantically wrong for the mesh
 * after assignment.
 */

ResourceManager& ResourceManager::operator=(const ResourceManager &other)
{
    if (this == &other) return *this;

    // make existing properties anonymous (invisible to mesh API users)
    clear_all_props();

    // Resize remaining, anonymized properties:
    for_each_entity([&](auto entity_tag) {
        resize_props<decltype(entity_tag)>(other.n<decltype(entity_tag)>());
    });

    clone_persistent_properties_from(other);


    return *this;
}

void ResourceManager::clone_persistent_properties_from(ResourceManager const& other)
{
  for_each_entity([&](auto entity_tag) {
      using ET = decltype(entity_tag);
      const auto &other_props = other.persistent_props_.get<ET>();
      auto &our_props = persistent_props_.get<ET>();
      for(const auto &p: other_props) {
        auto copy = p->clone();
        copy->set_tracker(&storage_tracker<ET>());
        our_props.insert(copy->shared_from_this());
      }
    });
}


detail::Tracker<PropertyStorageBase> &
ResourceManager::storage_tracker(EntityType type) const
{
    return storage_trackers_.get(type);
}

void ResourceManager::resize_vprops(size_t _nv) {
    resize_props<Entity::Vertex>(_nv);
}

void ResourceManager::resize_eprops(size_t _ne) {
    resize_props<Entity::Edge>(_ne);
    resize_props<Entity::HalfEdge>(2 * _ne);
}

void ResourceManager::resize_fprops(size_t _nf) {
    resize_props<Entity::Face>(_nf);
    resize_props<Entity::HalfFace>(2 * _nf);
}

void ResourceManager::resize_cprops(size_t _nc) {
    resize_props<Entity::Cell>(_nc);
}

void ResourceManager::reserve_vprops(size_t _n) {
    reserve_props<Entity::Vertex>(_n);
}
void ResourceManager::reserve_eprops(size_t _n) {
    reserve_props<Entity::Edge>(_n);
    reserve_props<Entity::HalfEdge>(2 * _n);
}
void ResourceManager::reserve_fprops(size_t _n) {
    reserve_props<Entity::Face>(_n);
    reserve_props<Entity::HalfFace>(2 * _n);
}
void ResourceManager::reserve_cprops(size_t _n) {
    reserve_props<Entity::Cell>(_n);
}


void ResourceManager::vertex_deleted(const VertexHandle& _h) {
    entity_deleted(_h);
}

void ResourceManager::edge_deleted(const EdgeHandle& _h) {
    entity_deleted(_h);
    entity_deleted(_h.halfedge_handle(1));
    entity_deleted(_h.halfedge_handle(0));
}

void ResourceManager::face_deleted(const FaceHandle& _h)
{
    entity_deleted(_h);
    entity_deleted(_h.halfface_handle(1));
    entity_deleted(_h.halfface_handle(0));
}

void ResourceManager::cell_deleted(const CellHandle& _h) {
    entity_deleted(_h);
}

void ResourceManager::clear_all_props()
{
    for_each_entity([this](auto entity_tag){ clear_props<decltype(entity_tag)>();});
}


template<> size_t OVM_EXPORT ResourceManager::n<Entity::Vertex>()   const { return n_vertices(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Edge>()     const { return n_edges(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::HalfEdge>() const { return n_halfedges(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Face>()     const { return n_faces(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::HalfFace>() const { return n_halffaces(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Cell>()     const { return n_cells(); }
template<> size_t OVM_EXPORT ResourceManager::n<Entity::Mesh>()     const { return 1; }



} // Namespace OpenVolumeMesh
