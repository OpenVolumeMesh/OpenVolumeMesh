#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <utility>

namespace OpenVolumeMesh {
template<typename T, typename U>
struct default_prop<std::pair<T, U>> {
    static inline const std::pair<T, U> value = {
        default_prop<T>::value,
        default_prop<U>::value
    };
};

}
