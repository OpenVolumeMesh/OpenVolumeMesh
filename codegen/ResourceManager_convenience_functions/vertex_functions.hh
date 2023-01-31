
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


