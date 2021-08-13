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
    *this = other;
}

ResourceManager::ResourceManager(ResourceManager &&other)
{
    *this = std::move(other);

}
ResourceManager& ResourceManager::operator=(const ResourceManager &other)
{
    if (this != &other) {
       assignAllPropertiesFrom<false>(other);
    }
    return *this;
}

ResourceManager& ResourceManager::operator=(ResourceManager &&other)
{
    if (this != &other) {
       assignAllPropertiesFrom<true>(other);
    }
    return *this;
}

detail::Tracker<PropertyStorageBase> &
ResourceManager::storage_tracker(EntityType type) const
{
    return storage_trackers_.get(type);
}

void ResourceManager::resize_vprops(size_t _nv) {
    resize_props(storage_tracker<Entity::Vertex>(), _nv);
}

void ResourceManager::resize_eprops(size_t _ne) {
    resize_props(storage_tracker<Entity::Edge>(), _ne);
    resize_props(storage_tracker<Entity::HalfEdge>(), 2 * _ne);
}

void ResourceManager::resize_fprops(size_t _nf) {
    resize_props(storage_tracker<Entity::Face>(), _nf);
    resize_props(storage_tracker<Entity::HalfFace>(), 2 * _nf);
}

void ResourceManager::resize_cprops(size_t _nc) {
    resize_props(storage_tracker<Entity::Cell>(), _nc);
}

void ResourceManager::reserve_vprops(size_t _n) {
    reserve_props(storage_tracker<Entity::Vertex>(), _n);
}
void ResourceManager::reserve_eprops(size_t _n) {
    reserve_props(storage_tracker<Entity::Edge>(), _n);
    reserve_props(storage_tracker<Entity::HalfEdge>(), 2 * _n);
}
void ResourceManager::reserve_fprops(size_t _n) {
    reserve_props(storage_tracker<Entity::Face>(), _n);
    reserve_props(storage_tracker<Entity::HalfFace>(), 2 * _n);
}
void ResourceManager::reserve_cprops(size_t _n) {
    reserve_props(storage_tracker<Entity::Cell>(), _n);
}


void ResourceManager::vertex_deleted(const VertexHandle& _h) {
    entity_deleted<Entity::Vertex>(_h);
}

void ResourceManager::edge_deleted(const EdgeHandle& _h) {
    entity_deleted<Entity::Edge>(_h);
    entity_deleted<Entity::HalfEdge>(_h.half(1));
    entity_deleted<Entity::HalfEdge>(_h.half(0));
}

void ResourceManager::face_deleted(const FaceHandle& _h)
{
    entity_deleted<Entity::Face>(_h);
    entity_deleted<Entity::HalfFace>(_h.half(1));
    entity_deleted<Entity::HalfFace>(_h.half(0));
}

void ResourceManager::cell_deleted(const CellHandle& _h) {
    entity_deleted<Entity::Cell>(_h);
}

void ResourceManager::clear_all_props()
{
    for_each_entity([this](auto entity_tag){ clear_props<decltype(entity_tag)>();});
}


template<> size_t ResourceManager::n<Entity::Vertex>()   const { return n_vertices(); }
template<> size_t ResourceManager::n<Entity::Edge>()     const { return n_edges(); }
template<> size_t ResourceManager::n<Entity::HalfEdge>() const { return n_halfedges(); }
template<> size_t ResourceManager::n<Entity::Face>()     const { return n_faces(); }
template<> size_t ResourceManager::n<Entity::HalfFace>() const { return n_halffaces(); }
template<> size_t ResourceManager::n<Entity::Cell>()     const { return n_cells(); }
template<> size_t ResourceManager::n<Entity::Mesh>()     const { return 1; }



} // Namespace OpenVolumeMesh
