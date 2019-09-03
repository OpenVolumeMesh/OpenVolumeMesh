#include "OpenVolumeMesh/Config/Version.hh"

// C++ version at OVM build time,
// to be used for things that change ABI

#if OPENVOLUMEMESH_CXX_VERSION >= 17
    #define OVM_BUILD_CXX_17 1
#else
    #define OVM_BUILD_CXX_17 0
#endif

// C++ version when compiling (e.g. client code),
// be careful not to change ABI depending on these defines.

#if __cplusplus >= 201703L
    #define OVM_CXX_17 1
#else
    #define OVM_CXX_17 0
#endif
