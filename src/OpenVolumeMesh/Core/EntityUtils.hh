#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <algorithm>

namespace OpenVolumeMesh {


/// Call templated function with a dummy argument of each entity type
template<typename F>
inline void for_each_entity(F fun) {
    fun(Entity::Vertex());
    fun(Entity::Edge());
    fun(Entity::HalfEdge());
    fun(Entity::Face());
    fun(Entity::HalfFace());
    fun(Entity::Cell());
    fun(Entity::Mesh());
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
