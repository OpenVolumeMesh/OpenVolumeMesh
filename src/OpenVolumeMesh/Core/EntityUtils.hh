#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <stdexcept>
#include <algorithm>
#include <array>


namespace OpenVolumeMesh {


/// Call templated function with a dummy argument of each entity type
template<typename F, typename ... Ts>
inline void for_each_entity(F fun, Ts... args) {
    fun(Entity::Vertex(),   args...);
    fun(Entity::Edge(),     args...);
    fun(Entity::HalfEdge(), args...);
    fun(Entity::Face(),     args...);
    fun(Entity::HalfFace(), args...);
    fun(Entity::Cell(),     args...);
    fun(Entity::Mesh(),     args...);
}

/// Call templated function with a dummy argument of the dynamic entity type
template<typename F, typename ... Ts>
inline auto entitytag_dispatch(EntityType type, F fun, Ts... args)
{
    switch (type) {
    case EntityType::Vertex:   return fun(Entity::Vertex(),   std::forward<Ts>(args)...);
    case EntityType::Edge:     return fun(Entity::Edge(),     std::forward<Ts>(args)...);
    case EntityType::HalfEdge: return fun(Entity::HalfEdge(), std::forward<Ts>(args)...);
    case EntityType::Face:     return fun(Entity::Face(),     std::forward<Ts>(args)...);
    case EntityType::HalfFace: return fun(Entity::HalfFace(), std::forward<Ts>(args)...);
    case EntityType::Cell:     return fun(Entity::Cell(),     std::forward<Ts>(args)...);
    case EntityType::Mesh:     return fun(Entity::Mesh(),     std::forward<Ts>(args)...);
    }
    throw std::runtime_error("entitytag_dispatch(): unknown entity.");
}


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
    std::array<T, n_entity_types> data_;
};

} // namespace OpenVolumeMesh
