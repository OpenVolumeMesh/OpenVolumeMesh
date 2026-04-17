#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <map>

namespace OpenVolumeMesh {
template<typename K, typename V>
struct default_prop<std::map<K, V>> {
    static inline const std::map<K, V> value = {};
};

}
