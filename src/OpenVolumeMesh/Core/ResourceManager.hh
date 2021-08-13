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
#include <optional>

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>
#include <OpenVolumeMesh/Core/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/TypeName.hh>
#include <OpenVolumeMesh/Core/ForwardDeclarations.hh>
#include <OpenVolumeMesh/Core/detail/Tracking.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyIterator.hh>


namespace OpenVolumeMesh {
class OVM_EXPORT ResourceManager {

public:
    using Properties = std::set<PropertyStorageBase*>;

    ResourceManager() = default;
    virtual ~ResourceManager() = default;

    ResourceManager(const ResourceManager &other);
    ResourceManager(ResourceManager &&other);
    ResourceManager& operator=(const ResourceManager &other);
    ResourceManager& operator=(ResourceManager &&other);

private:
    using PersistentProperties = std::set<std::shared_ptr<PropertyStorageBase>>;
    template<class EntityTag>
    void resize_props(size_t _n);

    template<class EntityTag>
    void reserve_props(size_t _n);

    template<class Handle>
    void entity_deleted(Handle);

    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> internal_find_property(const std::string& _name) const;

    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> internal_create_property(const std::string _name, const T _def = T()) const;


    template<bool Move, typename EntityTag>
    void assignProperties(typename std::conditional<Move, ResourceManager&, const ResourceManager&>::type src);
    template<bool Move>
    void assignAllPropertiesFrom(typename std::conditional<Move, ResourceManager&, const ResourceManager&>::type src);

    PerEntity<PersistentProperties> persistent_props_;

protected:
    friend class PropertyStorageBase;

    template<typename EntityTag>
    detail::Tracker<PropertyStorageBase> & storage_tracker() const;
    detail::Tracker<PropertyStorageBase> & storage_tracker(EntityType type) const;

    mutable PerEntity<detail::Tracker<PropertyStorageBase>> storage_trackers_;

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

    template <typename Handle>
    void swap_property_elements(Handle _idx_a, Handle _idx_b);

    template <typename Handle>
    void copy_property_elements(Handle _idx_a, Handle _idx_b);

public:
    /// drop all persistent properties.
    void clear_all_props();
    /// drop persistent properties.
    template<typename EntityTag> void clear_props();

    /// Get number of entities of given kind in mesh.
    template<typename EntityTag>
    size_t n() const;


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

    /// number of tracked properties
    template<typename EntityTag> size_t n_props() const;

    /// number of persistent properties
    template<typename EntityTag> size_t n_persistent_props() const;

    /** Get or create property: if the property does not exist yet, create it.
     */
    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> request_property(const std::string& _name = std::string(), const T _def = T());

    /** Create new property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> create_property(const std::string& _name = std::string(), const T _def = T());

    /** Create private property - useful for const meshes
     */
    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> create_private_property(std::string _name = {}, const T _def = T()) const;

    /** Get existing property: if the property does not exist, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> get_property(const std::string& _name = std::string());

    template <typename T, typename EntityTag>
    bool property_exists(const std::string& _name) const
    {
        return internal_find_property<T, EntityTag>(_name).has_value();
    }

    template<typename T, class EntityTag>
    void set_persistent(PropertyPtr<T, EntityTag>& _prop, bool _flag = true);

    template<typename EntityTag>
    PropertyIterator<PersistentProperties::const_iterator>
    persistent_props_begin() const
        {return persistent_props_.get<EntityTag>().cbegin();}

    template<typename EntityTag>
    PropertyIterator<PersistentProperties::const_iterator>
    persistent_props_end() const
        {return persistent_props_.get<EntityTag>().cend();}


/// convenience functions:

    [[deprecated("Use clear_props<Entity::Vertex>() instead.")]]
    inline void clear_vertex_props()   { clear_props<Entity::Vertex>();}
    [[deprecated("Use clear_props<Entity::Edge>() instead.")]]
    inline void clear_edge_props()     { clear_props<Entity::Edge>();}
    [[deprecated("Use clear_props<Entity::HalfEdge>() instead.")]]
    inline void clear_halfedge_props() { clear_props<Entity::HalfEdge>();}
    [[deprecated("Use clear_props<Entity::Face>() instead.")]]
    inline void clear_face_props()     { clear_props<Entity::Face>();}
    [[deprecated("Use clear_props<Entity::HalfFace>() instead.")]]
    inline void clear_halfface_props() { clear_props<Entity::HalfFace>();}
    [[deprecated("Use clear_props<Entity::Cell>() instead.")]]
    inline void clear_cell_props()     { clear_props<Entity::Cell>();}
    [[deprecated("Use clear_props<Entity::Mesh>() instead.")]]
    inline void clear_mesh_props()     { clear_props<Entity::Mesh>();}

    template<class T> VertexPropertyT<T>   request_vertex_property  (const std::string& _name = std::string(), const T _def = T());
    template<class T> EdgePropertyT<T>     request_edge_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name = std::string(), const T _def = T());
    template<class T> FacePropertyT<T>     request_face_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name = std::string(), const T _def = T());
    template<class T> CellPropertyT<T>     request_cell_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> MeshPropertyT<T>     request_mesh_property    (const std::string& _name = std::string(), const T _def = T());


    //[[deprecated("Use n_props<Entity::Vertex>() instead.")]]
    size_t n_vertex_props() const   { return n_props<Entity::Vertex>();}
    //[[deprecated("Use n_props<Entity::Edge instead.")]]
    size_t n_edge_props() const     { return n_props<Entity::Edge>();}
    //[[deprecated("Use n_props<Entity::HalfEdge>() instead.")]]
    size_t n_halfedge_props() const { return n_props<Entity::HalfEdge>();}
    //[[deprecated("Use n_props<Entity::Face>() instead.")]]
    size_t n_face_props() const     { return n_props<Entity::Face>();}
    //[[deprecated("Use n_props<Entity::HalfFace>() instead.")]]
    size_t n_halfface_props() const { return n_props<Entity::HalfFace>();}
    //[[deprecated("Use n_props<Entity::Cell>() instead.")]]
    size_t n_cell_props() const     { return n_props<Entity::Cell>();}
    //[[deprecated("Use n_props<Entity::Mesh>() instead.")]]
    size_t n_mesh_props() const     { return n_props<Entity::Mesh>();}

    [[deprecated("Use persistent_props_{begin,end}<Entity::Vertex>() instead.")]]
    auto vertex_props_begin()   const {return persistent_props_begin<Entity::Vertex>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Vertex>() instead.")]]
    auto vertex_props_end()     const {return persistent_props_end  <Entity::Vertex>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Edge>() instead.")]]
    auto edge_props_begin()     const {return persistent_props_begin<Entity::Edge>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Edge>() instead.")]]
    auto edge_props_end()       const {return persistent_props_end  <Entity::Edge>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::HalfEdge>() instead.")]]
    auto halfedge_props_begin() const {return persistent_props_begin<Entity::HalfEdge>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::HalfEdge>() instead.")]]
    auto halfedge_props_end()   const {return persistent_props_end  <Entity::HalfEdge>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Face>() instead.")]]
    auto face_props_begin()     const {return persistent_props_begin<Entity::Face>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Face>() instead.")]]
    auto face_props_end()       const {return persistent_props_end  <Entity::Face>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::HalfFace>() instead.")]]
    auto halfface_props_begin() const {return persistent_props_begin<Entity::HalfFace>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::HalfFace>() instead.")]]
    auto halfface_props_end()   const {return persistent_props_end  <Entity::HalfFace>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Cell>() instead.")]]
    auto cell_props_begin()     const {return persistent_props_begin<Entity::Cell>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Cell>() instead.")]]
    auto cell_props_end()       const {return persistent_props_end  <Entity::Cell>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Mesh>() instead.")]]
    auto mesh_props_begin()     const {return persistent_props_begin<Entity::Mesh>();}
    [[deprecated("Use persistent_props_{begin,end}<Entity::Mesh>() instead.")]]
    auto mesh_props_end()       const {return persistent_props_end  <Entity::Mesh>();}



    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::Vertex>() instead.")]]
    bool vertex_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Vertex>(_name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::Edge>() instead.")]]
    bool edge_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Edge>(_name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::HalfEdge>() instead.")]]
    bool halfedge_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::HalfEdge>(_name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::Face>() instead.")]]
    bool face_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Face>(_name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::HalfFace>() instead.")]]
    bool halfface_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::HalfFace>( _name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::Cell>() instead.")]]
    bool cell_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Cell>( _name);
    }

    template <class T>
    //[[deprecated("Use propery_exists<T, Entity::Mesh>() instead.")]]
    bool mesh_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Mesh>( _name);
    }




};

}

#include <OpenVolumeMesh/Core/ResourceManagerT_impl.hh>

