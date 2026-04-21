
#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <unordered_map>

namespace OpenVolumeMesh {
template<typename K, typename V>
struct default_prop<std::unordered_map<K, V>> {
    static inline const std::unordered_map<K, V> value = {};
};

}
