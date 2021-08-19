#pragma once

namespace OpenVolumeMesh::detail {

/// while std::vector<bool>::swap(ref,ref) exists, it (currently?)
/// does not work with _GLIBCXX_DEBUG!
void swap_bool(
        std::vector<bool>::reference a,
        std::vector<bool>::reference b)
{
    bool tmp = a;
    a = b;
    b = tmp;
}

} // namespace OpenVolumeMesh::detail
