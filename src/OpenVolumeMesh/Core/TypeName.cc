#include <string>
#include <typeinfo>

std::string get_type_name(std::type_info const &ti)
{
#ifdef _MSC_VER
    // MSVC's type_name only returns a friendly name with .name(),
    // get the more unique mangled name using .raw_name():
    return ti.raw_name();
#else
    // GCC and clang currently return the mangled name as .name(),
    // there is no .raw_name()
    return ti.name();
#endif
}

