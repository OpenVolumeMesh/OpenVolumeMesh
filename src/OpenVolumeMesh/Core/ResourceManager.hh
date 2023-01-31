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
#include <OpenVolumeMesh/Core/ForwardDeclarations.hh>
#include <OpenVolumeMesh/Core/detail/Tracking.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyIterator.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyStorageBase.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>


namespace OpenVolumeMesh {
class OVM_EXPORT ResourceManager {

public:
    using Properties = std::set<PropertyStorageBase*>;

    ResourceManager() = default;
    virtual ~ResourceManager() = default;

    ResourceManager(const ResourceManager &other);
    ResourceManager& operator=(const ResourceManager &other);

    ResourceManager(ResourceManager &&other) = delete;
    ResourceManager& operator=(ResourceManager &&other) = delete;

private:
    using PersistentProperties = std::set<std::shared_ptr<PropertyStorageBase>>;

    template <class T, typename EntityTag>
    static PropertyPtr<T, EntityTag> prop_ptr_from_storage(PropertyStorageBase *_prop);

    template<class EntityTag>
    void resize_props(size_t _n);

    template<class EntityTag>
    void reserve_props(size_t _n);

    template<class Handle>
    void entity_deleted(Handle);

    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> internal_find_property(const std::string& _name) const;

    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> internal_create_property(std::string _name, const T _def, bool shared) const;


    void clone_persistent_properties_from(ResourceManager const&);

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
     *  If it has an empty name, it will be a private property, otherwise shared.
     */
    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> request_property(const std::string& _name = std::string(), const T _def = T());

    /** Create new property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    [[deprecated("Use create_{shared,persistent}_property instead")]]
    std::optional<PropertyPtr<T, EntityTag>> create_property(std::string _name = std::string(), const T _def = T());

    /** Create new shared property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> create_shared_property(std::string _name, const T _def = T());

    /** Create new shared + persistent property: if the property already exists, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> create_persistent_property(std::string _name, const T _def = T());

    /** Create private property - useful for const meshes
     */
    template<typename T, typename EntityTag>
    PropertyPtr<T, EntityTag> create_private_property(std::string _name = {}, const T _def = T()) const;

    /** Get existing shared property. If the property does not exist, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<PropertyPtr<T, EntityTag>> get_property(const std::string& _name);

    /** Get a const view on an existing shared property. If the property does not exist, return no value.
     */
    template<typename T, typename EntityTag>
    std::optional<const PropertyPtr<T, EntityTag>> get_property(const std::string& _name) const;

    template <typename T, typename EntityTag>
    bool property_exists(const std::string& _name) const
    {
        return internal_find_property<T, EntityTag>(_name).has_value();
    }


    template<typename T, class EntityTag>
    void set_persistent(PropertyPtr<T, EntityTag>& _prop, bool _enable = true);

    template<typename T, class EntityTag>
    void set_shared(PropertyPtr<T, EntityTag>& _prop, bool _enable = true);

    template<typename EntityTag>
    PropertyIterator<PersistentProperties::const_iterator>
    persistent_props_begin() const
        {return persistent_props_.get<EntityTag>().cbegin();}

    template<typename EntityTag>
    PropertyIterator<PersistentProperties::const_iterator>
    persistent_props_end() const
        {return persistent_props_.get<EntityTag>().cend();}


/// convenience functions:

    inline void clear_vertex_props()   { clear_props<Entity::Vertex>();}
    inline void clear_edge_props()     { clear_props<Entity::Edge>();}
    inline void clear_halfedge_props() { clear_props<Entity::HalfEdge>();}
    inline void clear_face_props()     { clear_props<Entity::Face>();}
    inline void clear_halfface_props() { clear_props<Entity::HalfFace>();}
    inline void clear_cell_props()     { clear_props<Entity::Cell>();}
    inline void clear_mesh_props()     { clear_props<Entity::Mesh>();}

    // TODO: should we deprecate request_*_property, replace with get_or_create_*_property for more clarity?
    template<class T> VertexPropertyT<T>   request_vertex_property  (const std::string& _name = std::string(), const T _def = T());
    template<class T> EdgePropertyT<T>     request_edge_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> HalfEdgePropertyT<T> request_halfedge_property(const std::string& _name = std::string(), const T _def = T());
    template<class T> FacePropertyT<T>     request_face_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> HalfFacePropertyT<T> request_halfface_property(const std::string& _name = std::string(), const T _def = T());
    template<class T> CellPropertyT<T>     request_cell_property    (const std::string& _name = std::string(), const T _def = T());
    template<class T> MeshPropertyT<T>     request_mesh_property    (const std::string& _name = std::string(), const T _def = T());


    ////
    //// Vertex properties:
    ////

    size_t n_vertex_props()     const {return n_props<Entity::Vertex>();}
    auto vertex_props_begin()   const {return persistent_props_begin<Entity::Vertex>();}
    auto vertex_props_end()     const {return persistent_props_end  <Entity::Vertex>();}

    /** Create new shared vertex property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<VertexPropertyPtr<T>> create_shared_vertex_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::Vertex>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent vertex property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<VertexPropertyPtr<T>> create_persistent_vertex_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::Vertex>(std::move(_name), std::move(_def)); }

    /** Create private vertex property - useful for const meshes
     */
    template<typename T>
    VertexPropertyPtr<T> create_private_vertex_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::Vertex>(std::move(_name), std::move(_def)); }

    /** Get existing shared vertex property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<VertexPropertyPtr<T>> get_vertex_property(const std::string& _name) {
        return get_property<T, Entity::Vertex>(_name);
    }

    /** Get a const view on an existing shared vertex property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const VertexPropertyPtr<T>> get_vertex_property(const std::string& _name) const {
        return get_property<T, Entity::Vertex>(_name);
    }

    template <class T>
    bool vertex_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Vertex>(_name);
    }



    ////
    //// Edge properties:
    ////

    size_t n_edge_props()     const {return n_props<Entity::Edge>();}
    auto edge_props_begin()   const {return persistent_props_begin<Entity::Edge>();}
    auto edge_props_end()     const {return persistent_props_end  <Entity::Edge>();}

    /** Create new shared edge property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<EdgePropertyPtr<T>> create_shared_edge_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::Edge>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent edge property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<EdgePropertyPtr<T>> create_persistent_edge_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::Edge>(std::move(_name), std::move(_def)); }

    /** Create private edge property - useful for const meshes
     */
    template<typename T>
    EdgePropertyPtr<T> create_private_edge_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::Edge>(std::move(_name), std::move(_def)); }

    /** Get existing shared edge property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<EdgePropertyPtr<T>> get_edge_property(const std::string& _name) {
        return get_property<T, Entity::Edge>(_name);
    }

    /** Get a const view on an existing shared edge property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const EdgePropertyPtr<T>> get_edge_property(const std::string& _name) const {
        return get_property<T, Entity::Edge>(_name);
    }

    template <class T>
    bool edge_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Edge>(_name);
    }



    ////
    //// HalfEdge properties:
    ////

    size_t n_halfedge_props()     const {return n_props<Entity::HalfEdge>();}
    auto halfedge_props_begin()   const {return persistent_props_begin<Entity::HalfEdge>();}
    auto halfedge_props_end()     const {return persistent_props_end  <Entity::HalfEdge>();}

    /** Create new shared halfedge property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<HalfEdgePropertyPtr<T>> create_shared_halfedge_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::HalfEdge>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent halfedge property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<HalfEdgePropertyPtr<T>> create_persistent_halfedge_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::HalfEdge>(std::move(_name), std::move(_def)); }

    /** Create private halfedge property - useful for const meshes
     */
    template<typename T>
    HalfEdgePropertyPtr<T> create_private_halfedge_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::HalfEdge>(std::move(_name), std::move(_def)); }

    /** Get existing shared halfedge property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<HalfEdgePropertyPtr<T>> get_halfedge_property(const std::string& _name) {
        return get_property<T, Entity::HalfEdge>(_name);
    }

    /** Get a const view on an existing shared halfedge property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const HalfEdgePropertyPtr<T>> get_halfedge_property(const std::string& _name) const {
        return get_property<T, Entity::HalfEdge>(_name);
    }

    template <class T>
    bool halfedge_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::HalfEdge>(_name);
    }



    ////
    //// Face properties:
    ////

    size_t n_face_props()     const {return n_props<Entity::Face>();}
    auto face_props_begin()   const {return persistent_props_begin<Entity::Face>();}
    auto face_props_end()     const {return persistent_props_end  <Entity::Face>();}

    /** Create new shared face property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<FacePropertyPtr<T>> create_shared_face_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::Face>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent face property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<FacePropertyPtr<T>> create_persistent_face_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::Face>(std::move(_name), std::move(_def)); }

    /** Create private face property - useful for const meshes
     */
    template<typename T>
    FacePropertyPtr<T> create_private_face_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::Face>(std::move(_name), std::move(_def)); }

    /** Get existing shared face property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<FacePropertyPtr<T>> get_face_property(const std::string& _name) {
        return get_property<T, Entity::Face>(_name);
    }

    /** Get a const view on an existing shared face property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const FacePropertyPtr<T>> get_face_property(const std::string& _name) const {
        return get_property<T, Entity::Face>(_name);
    }

    template <class T>
    bool face_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Face>(_name);
    }



    ////
    //// HalfFace properties:
    ////

    size_t n_halfface_props()     const {return n_props<Entity::HalfFace>();}
    auto halfface_props_begin()   const {return persistent_props_begin<Entity::HalfFace>();}
    auto halfface_props_end()     const {return persistent_props_end  <Entity::HalfFace>();}

    /** Create new shared halfface property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<HalfFacePropertyPtr<T>> create_shared_halfface_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::HalfFace>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent halfface property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<HalfFacePropertyPtr<T>> create_persistent_halfface_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::HalfFace>(std::move(_name), std::move(_def)); }

    /** Create private halfface property - useful for const meshes
     */
    template<typename T>
    HalfFacePropertyPtr<T> create_private_halfface_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::HalfFace>(std::move(_name), std::move(_def)); }

    /** Get existing shared halfface property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<HalfFacePropertyPtr<T>> get_halfface_property(const std::string& _name) {
        return get_property<T, Entity::HalfFace>(_name);
    }

    /** Get a const view on an existing shared halfface property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const HalfFacePropertyPtr<T>> get_halfface_property(const std::string& _name) const {
        return get_property<T, Entity::HalfFace>(_name);
    }

    template <class T>
    bool halfface_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::HalfFace>(_name);
    }



    ////
    //// Cell properties:
    ////

    size_t n_cell_props()     const {return n_props<Entity::Cell>();}
    auto cell_props_begin()   const {return persistent_props_begin<Entity::Cell>();}
    auto cell_props_end()     const {return persistent_props_end  <Entity::Cell>();}

    /** Create new shared cell property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<CellPropertyPtr<T>> create_shared_cell_property(std::string _name, const T _def = T())
    { return create_shared_property<T, Entity::Cell>(std::move(_name), std::move(_def)); }

    /** Create new shared + persistent cell property: if the property already exists, return no value.
     */
    template<typename T>
    std::optional<CellPropertyPtr<T>> create_persistent_cell_property(std::string _name, const T _def = T())
    { return create_persistent_property<T, Entity::Cell>(std::move(_name), std::move(_def)); }

    /** Create private cell property - useful for const meshes
     */
    template<typename T>
    CellPropertyPtr<T> create_private_cell_property(std::string _name = {}, const T _def = T()) const
    { return create_private_property<T, Entity::Cell>(std::move(_name), std::move(_def)); }

    /** Get existing shared cell property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<CellPropertyPtr<T>> get_cell_property(const std::string& _name) {
        return get_property<T, Entity::Cell>(_name);
    }

    /** Get a const view on an existing shared cell property. If the property does not exist, return no value.
     */
    template<typename T>
    std::optional<const CellPropertyPtr<T>> get_cell_property(const std::string& _name) const {
        return get_property<T, Entity::Cell>(_name);
    }

    template <class T>
    bool cell_property_exists(const std::string& _name) const {
        return property_exists<T, Entity::Cell>(_name);
    }



};

}

#include <OpenVolumeMesh/Core/ResourceManagerT_impl.hh>

