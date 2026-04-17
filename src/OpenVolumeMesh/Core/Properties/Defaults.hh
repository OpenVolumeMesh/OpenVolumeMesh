#pragma once

#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <type_traits>

namespace OpenVolumeMesh {


/***
 * We require a safe default value for every property to avoid undefined behavior.
 * Some types that are popular as property types (OM/OVM/ACG VectorT, Eigen matrix types)
 * have default constructors that leave their elements uninitialized.
 * If we were to rely on value-initialization (as performed by std::vector), any
 * swap due to a deletion/garbage collection would constitute UB if some properties were not yet filled
 * in (e.g. a color property that is only occasionally used outside the core algorithm).
 *
 * Thus we only use value-initialisation for types were we explicitly know it to be safe,
 * and allow users to specify default values for other types they want to store in properties.
 *
 * Alternatively, users may always supply a default on every call to request_property/create_property.
 *
 * We provide default_prop specializations for common types here:
*   #include <OpenVolumeMesh/Geometry/EigenTraits.hh>
 */

template<typename T, typename=void>
struct default_prop;
#ifndef _MSC_VER // cl.exe does not print which specialization is missing with this approach 
template<typename T, typename>
struct default_prop
{
        // false value, but depending on T(!), to avoid static_assert(false) always failing:
        static_assert(!std::is_same_v<T, T>,
            "Please provide a default value for this property either in your {request,create}_property() call,\n"
            "or specialize OpenVolumeMesh::default_prop<T> for this property type.\n"
            "It can look something like this:\n\n"
            "template<>\n"
            "struct OpenVolumeMesh::default_prop<YourType> {\n"
            "    inline static const YourType value = 0; // set appropriate value\n"
            "};\n"
            "Note that OpenVolumeMesh provide default_prop specializations for some common types here:\n"
            "  #include <OpenVolumeMesh/Geometry/EigenTraits.hh>\n"
            "\nSorry for the inconvenience, we're trying to avoid undefined behavior.\n"
            );
};
#endif


// We don't want this to match user-defined types whose constructor
// does not initialize the members. This is very narrow, but better to err
// on the side of caution:
template<typename T>
struct default_prop<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
    static inline constexpr const T value = {};
};

template<typename Scalar, int DIM>
struct default_prop<VectorT<Scalar, DIM>> {
    static inline constexpr const VectorT<Scalar, DIM> value = VectorT<Scalar, DIM>::zero();
};

template<>
struct default_prop<std::string> {
    // not constexpr before c++20:
    static inline const std::string value = {};
};

template<typename T>
struct default_prop<std::vector<T>> {
    // not constexpr before c++20:
    static inline const std::vector<T> value = {};
};

template<typename HandleT>
struct default_prop<HandleT, std::enable_if_t<is_handle_v<HandleT>>>
{
    inline static constexpr HandleT value = HandleT::invalid();
};

#if 0
template<typename K, typename V>
struct default_prop<std::map<K, V>> {
    static inline const std::map<K, V> value = {};
};
#endif

template<typename T>
inline T const& default_prop_v = default_prop<T>::value;

} // namespace OpenVolumeMesh
