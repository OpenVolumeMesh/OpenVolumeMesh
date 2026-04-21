#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <tuple>

namespace OpenVolumeMesh {
template<typename... T>
struct default_prop<std::tuple<T...>> {
    static inline const std::tuple<T...> value = {
        default_prop<T>::value...
    };
};

}
