/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

// TODO: rename to Handles.hh
#pragma once

#include <algorithm>
#include <iosfwd>
#include <vector>
#include <cassert>
#include <limits>
#include <functional>
#include <type_traits>

#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/System/Deprecation.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {
// Define handle types in order to distinguish different entities by their indices

namespace detail {

template <typename Derived>
class OVM_EXPORT HandleT
{
public:
    constexpr HandleT() = default;
    explicit constexpr HandleT(int _idx) : idx_(_idx) {}

    HandleT(const HandleT& _idx) = default;
    HandleT(HandleT&& _idx) = default;
    HandleT& operator=(const HandleT& _idx) = default;
    HandleT& operator=(HandleT&& _idx) = default;

    [[deprecated]]
    HandleT& operator=(int _idx) {
        idx_ = _idx;
        return *this;
    }

    bool is_valid() const { return idx_ != -1; }

    constexpr bool operator<(const Derived& _idx) const { return (this->idx_ < _idx.idx_); }

    constexpr bool operator>(const Derived& _idx) const { return (this->idx_ > _idx.idx_); }

    constexpr bool operator==(const Derived& _h) const { return _h.idx_ == this->idx_; }

    constexpr bool operator!=(const Derived& _h) const { return _h.idx_ != this->idx_; }

    constexpr const int& idx() const { return idx_; }
    constexpr int& idx_mutable() { return idx_; }

    /// return unsigned idx - handle must be valid
    unsigned int uidx() const { assert(is_valid()); return static_cast<size_t>(idx_); }

    /// set idx
	void idx(const int& _idx) { idx_ = _idx; }

    /// make handle invalid
	void reset() { idx_ = -1; }

    static Derived from_unsigned(size_t _idx)
    {
        if (_idx <= static_cast<size_t>(std::numeric_limits<int>::max())) {
            return Derived{static_cast<int>(_idx)};
        } else {
            assert(false);
            return Derived{};
        }
    }

private:
    int idx_ = -1;
};

template <typename Derived, typename SuperHandle>
class SubHandleT : public HandleT<Derived>
{
public:
    using HandleT<Derived>::HandleT;

    constexpr int subidx() const {
        assert(this->is_valid());
        return this->idx() & 1;
    }

    constexpr SuperHandle full() const {
        assert(this->is_valid());
        return SuperHandle{this->idx() / 2 };
    }

    /// obtain opposite sub-handle.
    constexpr Derived opp() const {
        assert(this->is_valid());
        return Derived{this->idx() ^ 1};
    }

};

template <typename Derived, typename SubHandle>
class SuperHandleT : public HandleT<Derived>
{
public:
    using HandleT<Derived>::HandleT;

    constexpr SubHandle half(int subidx) const {
        assert(0 <= subidx && subidx <= 1);
        return SubHandle{2 * HandleT<Derived>::idx() + subidx};
    }

};

} // namespace detail

class VH  : public detail::HandleT<VH> {
public:
    using detail::HandleT<VH>::HandleT;
    using EntityTag = Entity::Vertex;
};

class HEH;
class EH  : public detail::SuperHandleT<EH, HEH> {
public:
    using detail::SuperHandleT<EH, HEH>::SuperHandleT;
    using EntityTag = Entity::Edge;

};
class HEH : public detail::SubHandleT<HEH, EH> {
public:
    using detail::SubHandleT<HEH, EH>::SubHandleT;
    using EntityTag = Entity::HalfEdge;
};

class HFH;
class FH  : public detail::SuperHandleT<FH, HFH> {
public:
    using detail::SuperHandleT<FH, HFH>::SuperHandleT;
    using EntityTag = Entity::Face;

};
class HFH : public detail::SubHandleT<HFH, FH> {
public:
    using detail::SubHandleT<HFH, FH>::SubHandleT;
    using EntityTag = Entity::HalfFace;
};

class CH  : public detail::HandleT<CH> {
public:
    using detail::HandleT<CH>::HandleT;
    using EntityTag = Entity::Cell;

};
class MH  : public detail::HandleT<MH> {
public:
    using detail::HandleT<MH>::HandleT;
    using EntityTag = Entity::Mesh;
};

using VertexHandle = VH;
using EdgeHandle = EH;
using HalfEdgeHandle = HEH;
using FaceHandle = FH;
using HalfFaceHandle = HFH;
using CellHandle = CH;
using MeshHandle = MH;

template<typename EntityTag>
struct handle_for_tag;


template<> struct handle_for_tag<Entity::Vertex>   { using type = VH;  };
template<> struct handle_for_tag<Entity::Edge>     { using type = EH;  };
template<> struct handle_for_tag<Entity::HalfEdge> { using type = HEH; };
template<> struct handle_for_tag<Entity::Face>     { using type = FH;  };
template<> struct handle_for_tag<Entity::HalfFace> { using type = HFH; };
template<> struct handle_for_tag<Entity::Cell>     { using type = CH;  };
template<> struct handle_for_tag<Entity::Mesh>     { using type = MH;  };


template<typename EntityTag>
using HandleT = typename handle_for_tag<EntityTag>::type;

template<typename>
struct is_handle : public std::false_type {};

template<> struct is_handle<VH>  : public std::true_type {};
template<> struct is_handle<EH>  : public std::true_type {};
template<> struct is_handle<HEH> : public std::true_type {};
template<> struct is_handle<FH>  : public std::true_type {};
template<> struct is_handle<HFH> : public std::true_type {};
template<> struct is_handle<CH>  : public std::true_type {};
template<> struct is_handle<MH>  : public std::true_type {};

template<typename Handle>
inline const bool is_handle_v = is_handle<Handle>::value;

std::ostream& operator<<(std::ostream& _ostr, VH _h);
std::ostream& operator<<(std::ostream& _ostr, EH _h);
std::ostream& operator<<(std::ostream& _ostr, HEH _h);
std::ostream& operator<<(std::ostream& _ostr, FH _h);
std::ostream& operator<<(std::ostream& _ostr, HFH _h);
std::ostream& operator<<(std::ostream& _ostr, CH _h);
std::ostream& operator<<(std::ostream& _ostr, MH _h);

std::istream& operator>>(std::istream& _istr, VH &_h);
std::istream& operator>>(std::istream& _istr, EH &_h);
std::istream& operator>>(std::istream& _istr, HEH &_h);
std::istream& operator>>(std::istream& _istr, FH &_h);
std::istream& operator>>(std::istream& _istr, HFH &_h);
std::istream& operator>>(std::istream& _istr, CH &_h);
std::istream& operator>>(std::istream& _istr, MH &_h);

} // Namespace OpenVolumeMesh

namespace std
{
    template<> struct hash<OpenVolumeMesh::VH> {
        auto operator()(OpenVolumeMesh::VH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::EH> {
        auto operator()(OpenVolumeMesh::EH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::HEH> {
        auto operator()(OpenVolumeMesh::HEH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::FH> {
        auto operator()(OpenVolumeMesh::FH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::HFH> {
        auto operator()(OpenVolumeMesh::HFH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::CH> {
        auto operator()(OpenVolumeMesh::CH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
    template<> struct hash<OpenVolumeMesh::MH> {
        auto operator()(OpenVolumeMesh::MH const& _h) const noexcept
        { return std::hash<int>{}(_h.idx()); }
    };
}

