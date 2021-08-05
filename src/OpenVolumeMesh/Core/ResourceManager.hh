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


#ifndef NDEBUG
#include <iostream>
#endif
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <type_traits>

#include <OpenVolumeMesh/Core/Entities.hh>
#include "OpenVolumeMesh/Config/Export.hh"
#include "PropertyStorageT.hh"
#include "PropertyHandles.hh"
#include "TypeName.hh"
#include "ForwardDeclarations.hh"

#if OVM_CXX_17
#include <optional>
#endif

namespace OpenVolumeMesh {

template<typename T>
class PerEntity
{
public:
    template<typename EntityTag>
    inline T const & get() const {return get(EntityTag::type());}
    template<typename EntityTag>
    inline T &get() {return get(EntityTag::type());}

    inline T &get(EntityType type) {
        return data_[static_cast<size_t>(type)];
    }
    inline T const &get(EntityType type) const {
        return data_[static_cast<size_t>(type)];
    }
    template<typename F>
    auto for_each(F &f) {
        std::for_each(data_.begin(), data_.end(), f);
    }
private:
    std::array<T, 7> data_;
};

// Forward declarations
class BaseProperty;

class OVM_EXPORT ResourceManager {
public:
    virtual ~ResourceManager() = default;
private:
    using Properties = std::set<PropertyStorageBase*>;
    using PersistentProperties = std::set<std::shared_ptr<PropertyStorageBase>>;


protected:
    /// Change size of stored vertex properties
    void resize_vprops(size_t _nv);

    /// Change size of stored edge properties
    void resize_eprops(size_t _ne);

    /// Change size of stored face properties
    void resize_fprops(size_t _nf);

    /// Change size of stored cell properties
    void resize_cprops(size_t _nc);

    void reserve_vprops(size_t n);
    void reserve_eprops(size_t n);
    void reserve_fprops(size_t n);
    void reserve_cprops(size_t n);

    void vertex_deleted(const VertexHandle& _h);
    void edge_deleted(const EdgeHandle& _h);
    void face_deleted(const FaceHandle& _h);
    void cell_deleted(const CellHandle& _h);

    void swap_cell_properties(CellHandle _h1, CellHandle _h2);
    void swap_face_properties(FaceHandle _h1, FaceHandle _h2);
    void swap_halfface_properties(HalfFaceHandle _h1, HalfFaceHandle _h2);
    void swap_edge_properties(EdgeHandle _h1, EdgeHandle _h2);
    void swap_halfedge_properties(HalfEdgeHandle _h1, HalfEdgeHandle _h2);
    void swap_vertex_properties(VertexHandle _h1, VertexHandle _h2);

    template <typename EntityTag>
    void swap_property_elements(HandleT<EntityTag> const &_idx_a, HandleT<EntityTag> &_idx_b);

public:
    template<typename EntityTag> void clear_props();
    inline void clear_vertex_props()   { clear_props<Entity::Vertex>();}
    inline void clear_edge_props()     { clear_props<Entity::Edge>();}
    inline void clear_halfedge_props() { clear_props<Entity::HalfEdge>();}
    inline void clear_face_props()     { clear_props<Entity::Face>();}
    inline void clear_halfface_props() { clear_props<Entity::HalfFace>();}
    inline void clear_cell_props()     { clear_props<Entity::Cell>();}
    inline void clear_mesh_props()     { clear_props<Entity::Mesh>();}

    /// Get number of vertices in mesh
    virtual size_t n_vertices() const = 0;
    /// Get number of edges in mesh
    virtual size_t n_edges() const = 0;
    /// Get number of halfedges in mesh
    virtual size_t n_halfedges() const = 0;
    /// Get number of faces in mesh
    virtual size_t n_faces() const = 0;
    /// Get number of halffaces in mesh
    virtual size_t n_halffaces() const = 0;
    /// Get number of cells in mesh
    virtual size_t n_cells() const = 0;

    /** Get or create property: if the property does not exist yet, create it.
     */
    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> request_property(const std::string& _name = std::string(), const T _def = T());

    /** Create new property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> create_property(const std::string& _name = std::string(), const T _def = T());

    /** Get existing property: if the property does not exist, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> get_property(const std::string& _name = std::string());

    template<class T> VertexPropertyT<T> request_vertex_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> EdgePropertyT<T> request_edge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> FacePropertyT<T> request_face_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> CellPropertyT<T> request_cell_property(const std::string& _name = std::string(), const T _def = T());

    template<class T> MeshPropertyT<T> request_mesh_property(const std::string& _name = std::string(), const T _def = T());


public:

    template<typename EntityTag> size_t n_props() const;

    size_t n_vertex_props() const   { return n_props<Entity::Vertex>();}
    size_t n_edge_props() const     { return n_props<Entity::Edge>();}
    size_t n_halfedge_props() const { return n_props<Entity::HalfEdge>();}
    size_t n_face_props() const     { return n_props<Entity::Face>();}
    size_t n_halfface_props() const { return n_props<Entity::HalfFace>();}
    size_t n_cell_props() const     { return n_props<Entity::Cell>();}
    size_t n_mesh_props() const     { return n_props<Entity::Mesh>();}

    template<typename T, class EntityTag>
    void set_persistent(PropertyPtr<T, EntityTag>& _prop, bool _flag = true);


    // TODO: - make custom iterator to hide underlying container
    //       - implement iteration over props for all entities

    template<typename EntityTag>
    PersistentProperties::const_iterator persistent_props_begin() const
        {return persistent_props_.get<EntityTag>().cbegin();}

    template<typename EntityTag>
    PersistentProperties::const_iterator persistent_props_end() const
        {return persistent_props_.get<EntityTag>().cend();}

public:
    template <class PropT, class EntityTag>
    bool property_exists(const std::string& _name) const
    {
        return internal_find_property<PropT, EntityTag>(_name).has_value();
    }

    template <class PropT>
    bool vertex_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::Vertex>(_name);
    }

    template <class PropT>
    bool edge_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::Edge>(_name);
    }

    template <class PropT>
    bool halfedge_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::HalfEdge>(_name);
    }

    template <class PropT>
    bool face_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::Face>(_name);
    }

    template <class PropT>
    bool halfface_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::HalfFace>( _name);
    }

    template <class PropT>
    bool cell_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::Cell>( _name);
    }

    template <class PropT>
    bool mesh_property_exists(const std::string& _name) const {
        return property_exists<PropT, Entity::Mesh>( _name);
    }

protected:

    template<typename EntityTag>
    void delete_multiple_entities(const std::vector<bool>& _tags);

    void delete_multiple_vertex_props(const std::vector<bool>& _tags);

    void delete_multiple_edge_props(const std::vector<bool>& _tags);

    void delete_multiple_face_props(const std::vector<bool>& _tags);

    void delete_multiple_cell_props(const std::vector<bool>& _tags);

private:

    template<class StdVecT>
    void resize_props(StdVecT& _vec, size_t _n);

    template<class StdVecT>
    void reserve_props(StdVecT& _vec, size_t _n);

    template<class StdVecT>
    void entity_deleted(StdVecT& _vec, const OpenVolumeMeshHandle& _h);

    template<class StdVecT>
    void remove_property(StdVecT& _vec, size_t _idx);

    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> internal_find_property(const std::string& _name) const;

    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> internal_create_property(const std::string& _name, const T _def = T());

    PerEntity<Properties> weak_props_;
    PerEntity<PersistentProperties> persistent_props_;

    template<typename Entity>
    Properties &entity_props() const;

    template<typename Entity>
    size_t n();


private:
    template<typename T, typename EntityTag>
    friend class PropertyPtr;

    template<typename EntityTag>
    void property_deleted(PropertyStorageBase *ptr) {
        weak_props_.get<EntityTag>().erase(ptr);
    }
};

}

#include "ResourceManagerT_impl.hh"

