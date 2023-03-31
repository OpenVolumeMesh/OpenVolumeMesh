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

#pragma once

#include <algorithm>
#include <iosfwd>
#include <vector>
#include <cassert>
#include <limits>
#include <functional>

#include <OpenVolumeMesh/Core/Entities.hh>
#include <functional>
#include <OpenVolumeMesh/System/Deprecation.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {


namespace detail {

class OVM_EXPORT HandleBase
{
public:
    constexpr HandleBase() = default;
    explicit constexpr HandleBase(int _idx) : idx_(_idx) {}

    constexpr HandleBase(const HandleBase& _idx) = default;
    constexpr HandleBase(HandleBase&& _idx) = default;
    HandleBase& operator=(const HandleBase& _idx) = default;
    HandleBase& operator=(HandleBase&& _idx) = default;

    [[deprecated]]
    HandleBase& operator=(int _idx) {
        idx_ = _idx;
        return *this;
    }

    bool is_valid() const { return idx_ >= 0; }

    constexpr bool operator<(const HandleBase& _idx) const { return (this->idx_ < _idx.idx_); }

    constexpr bool operator>(const HandleBase& _idx) const { return (this->idx_ > _idx.idx_); }

    constexpr bool operator==(const HandleBase& _h) const { return _h.idx_ == this->idx_; }

    constexpr bool operator!=(const HandleBase& _h) const { return _h.idx_ != this->idx_; }

    constexpr const int& idx() const { return idx_; }
    constexpr int& idx_mutable() { return idx_; }

    /// return unsigned idx - handle must be valid
    unsigned int uidx() const { assert(is_valid()); return static_cast<size_t>(idx_); }

    /// set idx; TODO: deprecate, provides no advantage over `h = VH(123);`
    void idx(const int& _idx) { idx_ = _idx; }

    /// make handle invalid
    void reset() { idx_ = -1; }
private:
    int idx_ = -1;
};

template <typename Derived>
class OVM_EXPORT HandleT : public HandleBase
{
public:
    using HandleBase::HandleBase;
    static Derived from_unsigned(size_t _idx)
    {
        if (_idx <= static_cast<size_t>(std::numeric_limits<int>::max())) {
            return Derived{static_cast<int>(_idx)};
        } else {
            assert(false);
            return Derived{};
        }
    }

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

protected:
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

protected:
    constexpr SubHandle half(int subidx) const {
        assert(0 <= subidx && subidx <= 1);
        return SubHandle{2 * HandleT<Derived>::idx() + subidx};
    }

};

} // namespace detail


class OVM_EXPORT VH  : public detail::HandleT<VH> {
public:
    using detail::HandleT<VH>::HandleT;
    using EntityTag = Entity::Vertex;
};

class HEH;
class OVM_EXPORT EH  : public detail::SuperHandleT<EH, HEH> {
public:
    using detail::SuperHandleT<EH, HEH>::SuperHandleT;
    using EntityTag = Entity::Edge;

    constexpr HEH halfedge_handle(int subidx) const;

};

class OVM_EXPORT HEH : public detail::SubHandleT<HEH, EH> {
public:
    using detail::SubHandleT<HEH, EH>::SubHandleT;
    using EntityTag = Entity::HalfEdge;

    constexpr HEH opposite_handle() const;
    constexpr EH edge_handle() const;
};


class HFH;
class OVM_EXPORT FH  : public detail::SuperHandleT<FH, HFH> {
public:
    using detail::SuperHandleT<FH, HFH>::SuperHandleT;
    using EntityTag = Entity::Face;

    constexpr HFH halfface_handle(int subidx) const;

};

class OVM_EXPORT HFH : public detail::SubHandleT<HFH, FH> {
public:
    using detail::SubHandleT<HFH, FH>::SubHandleT;
    using EntityTag = Entity::HalfFace;

    constexpr HFH opposite_handle() const;
    constexpr FH face_handle() const;
};


class OVM_EXPORT CH  : public detail::HandleT<CH> {
public:
    using detail::HandleT<CH>::HandleT;
    using EntityTag = Entity::Cell;

};
class OVM_EXPORT MH  : public detail::HandleT<MH> {
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

constexpr inline HFH FH::halfface_handle(int subidx) const {
    return half(subidx);
}

constexpr inline HEH EH::halfedge_handle(int subidx) const {
    return half(subidx);
}

constexpr inline HEH HEH::opposite_handle() const {
    return opp();
}

constexpr inline EH HEH::edge_handle() const {
    return full();
}

constexpr inline HFH HFH::opposite_handle() const {
    return opp();
}

constexpr inline FH HFH::face_handle() const {
    return full();
}

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
inline constexpr bool is_handle_v = is_handle<Handle>::value;

OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, VH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, EH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, HEH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, FH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, HFH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, CH _h);
OVM_EXPORT std::ostream& operator<<(std::ostream& _ostr, MH _h);

OVM_EXPORT std::istream& operator>>(std::istream& _istr, VH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, EH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, HEH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, FH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, HFH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, CH &_h);
OVM_EXPORT std::istream& operator>>(std::istream& _istr, MH &_h);

// Helper class that is used to decrease all handles
// exceeding a certain threshold

class VHandleCorrection {
public:
    explicit VHandleCorrection(VertexHandle _thld) : thld_(_thld) {}
    void correctValue(VertexHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 1);
    }
private:
    VertexHandle thld_;
};
class HEHandleCorrection {
public:
    explicit HEHandleCorrection(HalfEdgeHandle _thld) : thld_(_thld) {}
    void correctVecValue(std::vector<HalfEdgeHandle>& _vec) {
#if defined(__clang_major__) && (__clang_major__ >= 5)
        for(std::vector<HalfEdgeHandle>::iterator it = _vec.begin(), end = _vec.end(); it != end; ++it) {
            correctValue(*it);
        }
#else
        std::for_each(_vec.begin(), _vec.end(), std::bind(&HEHandleCorrection::correctValue, this, std::placeholders::_1));
#endif
    }
    void correctValue(HalfEdgeHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 2);
    }
private:
    HalfEdgeHandle thld_;
};
class HFHandleCorrection {
public:
    explicit HFHandleCorrection(HalfFaceHandle _thld) : thld_(_thld) {}
    void correctVecValue(std::vector<HalfFaceHandle>& _vec) {
#if defined(__clang_major__) && (__clang_major__ >= 5)
        for(std::vector<HalfFaceHandle>::iterator it = _vec.begin(), end = _vec.end(); it != end; ++it) {
            correctValue(*it);
        }
#else
        std::for_each(_vec.begin(), _vec.end(), std::bind(&HFHandleCorrection::correctValue, this, std::placeholders::_1));
#endif
    }
    void correctValue(HalfFaceHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 2);
    }
private:
    HalfFaceHandle thld_;
};
class CHandleCorrection {
public:
    explicit CHandleCorrection(CellHandle _thld) : thld_(_thld) {}
    void correctValue(CellHandle& _h) {
        if(_h > thld_) _h.idx(_h.idx() - 1);
    }
private:
    CellHandle thld_;
};

} // Namespace OpenVolumeMesh

template<> struct std::hash<OpenVolumeMesh::VH>
{
    std::size_t operator()(OpenVolumeMesh::VH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
template<> struct std::hash<OpenVolumeMesh::EH>
{
    std::size_t operator()(OpenVolumeMesh::EH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
template<> struct std::hash<OpenVolumeMesh::HEH>
{
    std::size_t operator()(OpenVolumeMesh::HEH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
template<> struct std::hash<OpenVolumeMesh::FH>
{
    std::size_t operator()(OpenVolumeMesh::FH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
template<> struct std::hash<OpenVolumeMesh::HFH>
{
    std::size_t operator()(OpenVolumeMesh::HFH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
template<> struct std::hash<OpenVolumeMesh::CH>
{
    std::size_t operator()(OpenVolumeMesh::CH const& h) const noexcept
    { return std::hash<int>{}(h.idx()); }
};
