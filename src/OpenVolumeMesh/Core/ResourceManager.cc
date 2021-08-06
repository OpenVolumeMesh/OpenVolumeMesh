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
#include <OpenVolumeMesh/Core/BaseProperty.hh>

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

void ResourceManager::resize_vprops(size_t _nv) {
    resize_props(all_props_.get<Entity::Vertex>(), _nv);
}

void ResourceManager::resize_eprops(size_t _ne) {
    resize_props(all_props_.get<Entity::Edge>(), _ne);
    resize_props(all_props_.get<Entity::HalfEdge>(), 2 * _ne);
}

void ResourceManager::resize_fprops(size_t _nf) {
    resize_props(all_props_.get<Entity::Face>(), _nf);
    resize_props(all_props_.get<Entity::HalfFace>(), 2 * _nf);
}

void ResourceManager::resize_cprops(size_t _nc) {
    resize_props(all_props_.get<Entity::Cell>(), _nc);
}

void ResourceManager::reserve_vprops(size_t _n) {
    reserve_props(all_props_.get<Entity::Vertex>(), _n);
}
void ResourceManager::reserve_eprops(size_t _n) {
    reserve_props(all_props_.get<Entity::Edge>(), _n);
    reserve_props(all_props_.get<Entity::HalfEdge>(), 2 * _n);
}
void ResourceManager::reserve_fprops(size_t _n) {
    reserve_props(all_props_.get<Entity::Face>(), _n);
    reserve_props(all_props_.get<Entity::HalfFace>(), 2 * _n);
}
void ResourceManager::reserve_cprops(size_t _n) {
    reserve_props(all_props_.get<Entity::Cell>(), _n);
}


void ResourceManager::vertex_deleted(const VertexHandle& _h) {
    entity_deleted(all_props_.get<Entity::Vertex>(), _h);
}

void ResourceManager::edge_deleted(const EdgeHandle& _h) {

    entity_deleted(all_props_.get<Entity::Edge>(), _h);
    entity_deleted(all_props_.get<Entity::HalfEdge>(), HalfEdgeHandle{_h.idx()*2+1});
    entity_deleted(all_props_.get<Entity::HalfEdge>(), HalfEdgeHandle{_h.idx()*2});
}

void ResourceManager::face_deleted(const FaceHandle& _h)
{
    entity_deleted(all_props_.get<Entity::Face>(), _h);
    entity_deleted(all_props_.get<Entity::HalfFace>(), HalfFaceHandle{_h.idx()*2+1});
    entity_deleted(all_props_.get<Entity::HalfFace>(), HalfFaceHandle{_h.idx()*2});
}

void ResourceManager::cell_deleted(const CellHandle& _h) {
    entity_deleted(all_props_.get<Entity::Cell>(), _h);
}

void ResourceManager::swap_cell_properties(CellHandle _h1, CellHandle _h2)
{
    swap_property_elements<Entity::Cell>(_h1, _h2);
}

void ResourceManager::swap_face_properties(FaceHandle _h1, FaceHandle _h2)
{
    swap_property_elements<Entity::Face>(_h1, _h2);
}

void ResourceManager::swap_halfface_properties(HalfFaceHandle _h1, HalfFaceHandle _h2)
{
    swap_property_elements<Entity::HalfFace>(_h1, _h2);
}

void ResourceManager::swap_edge_properties(EdgeHandle _h1, EdgeHandle _h2)
{
    swap_property_elements<Entity::Edge>(_h1, _h2);
}

void ResourceManager::swap_halfedge_properties(HalfEdgeHandle _h1, HalfEdgeHandle _h2){

    swap_property_elements<Entity::HalfEdge>(_h1, _h2);
}

void ResourceManager::swap_vertex_properties(VertexHandle _h1, VertexHandle _h2){

    swap_property_elements<Entity::Vertex>(_h1, _h2);
}
void ResourceManager::clear_all_props()
{
    for_each_entity([this](auto entity_tag){ clear_props<decltype(entity_tag)>();});
}

void ResourceManager::delete_multiple_vertex_props(const std::vector<bool>& _tags)
{
    delete_multiple_entities<Entity::Vertex>(_tags);
}

void ResourceManager::delete_multiple_edge_props(const std::vector<bool>& _tags)
{
    delete_multiple_entities<Entity::Edge>(_tags);

    // TODO OPTI: remove_if stuff? see delete_multiple_face_props too
    // Create tags vector for halfedges
    std::vector<bool> hetags;
    for(std::vector<bool>::const_iterator t_it = _tags.begin(),
            t_end = _tags.end(); t_it != t_end; ++t_it) {
        hetags.push_back(*t_it);
        hetags.push_back(*t_it);
    }
    delete_multiple_entities<Entity::HalfEdge>(hetags);
}

void ResourceManager::delete_multiple_face_props(const std::vector<bool>& _tags)
{

    delete_multiple_entities<Entity::Face>(_tags);

    // Create tags vector for halffaces
    std::vector<bool> hftags;
    for(std::vector<bool>::const_iterator t_it = _tags.begin(),
            t_end = _tags.end(); t_it != t_end; ++t_it) {
        hftags.push_back(*t_it);
        hftags.push_back(*t_it);
    }
    delete_multiple_entities<Entity::HalfFace>(hftags);
}

void ResourceManager::delete_multiple_cell_props(const std::vector<bool>& _tags)
{
    delete_multiple_entities<Entity::Cell>(_tags);
}


template<> size_t ResourceManager::n<Entity::Vertex>()   const { return n_vertices(); }
template<> size_t ResourceManager::n<Entity::Edge>()     const { return n_edges(); }
template<> size_t ResourceManager::n<Entity::HalfEdge>() const { return n_halfedges(); }
template<> size_t ResourceManager::n<Entity::Face>()     const { return n_faces(); }
template<> size_t ResourceManager::n<Entity::HalfFace>() const { return n_halffaces(); }
template<> size_t ResourceManager::n<Entity::Cell>()     const { return n_cells(); }
template<> size_t ResourceManager::n<Entity::Mesh>()     const { return 1; }



} // Namespace OpenVolumeMesh
