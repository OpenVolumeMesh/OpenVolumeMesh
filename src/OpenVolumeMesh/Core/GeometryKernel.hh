#pragma once
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

#include <cassert>
#include <iostream>
#include <type_traits>

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

#if 0
template <class VecT>
class GeometryKernelT {
public:
    GeometryKernelT(ResourceManager *_resman)
        : prop_{_resman->create_persistent_property<VecT, Entity::Vertex>("ovm:position")}
    {}
    VecT const& operator[](VH vh) const {return prop_[vh];}
    VecT      & operator[](VH vh)       {return prop_[vh];}

private:
    PropertyPtr<VecT, Entity::Vertex> prop_;
};
#endif
template <class VecT>
using GeometryKernelT = PropertyPtr<VecT, Entity::Vertex>;

template <class VecT, class TopologyKernelT = TopologyKernel>
class GeometryKernel : public TopologyKernelT {
public:

    using PointT = VecT; // OVM legacy
    using Point = VecT;  // OpenMesh compatiblity
    using KernelT = TopologyKernelT;

    GeometryKernel()
        : TopologyKernelT()
        , position_{get_prop()}
    {}

    GeometryKernel(GeometryKernel&& other)
        : TopologyKernelT(std::move(other))
        , position_{get_prop()}
    {}

    GeometryKernel(GeometryKernel const& other)
        : TopologyKernelT(other)
        , position_{get_prop()}
    {}

    GeometryKernel& operator=(GeometryKernel &&other)
    {
        TopologyKernelT::operator=(std::move(other));
        position_ = get_prop();
    }

    GeometryKernel& operator=(GeometryKernel const&other)
    {
        TopologyKernelT::operator=(other);
        position_ = get_prop();
    }

    template<class OtherTopoKernel>
    void assign(const GeometryKernel<VecT, OtherTopoKernel> *other)
    {
        TopologyKernelT::assign(other);
        position_.assign(other->vertex_positions());
    }

    /// Add a geometric point to the mesh
    VertexHandle add_vertex(const VecT& _p) {
        VH vh = TopologyKernelT::add_vertex();
        position_[vh] = _p;
        return vh;
    }

    /// Set the coordinates of point _vh
    [[deprecated("Use VecT& vertex(VH) instead.")]]
    void set_vertex(VertexHandle _vh, const VecT& _p) {

        assert(_vh.idx() < (int)position_.size());

        position_[_vh] = _p;
    }

    /// Get point _vh's coordinates
    const VecT& vertex(VertexHandle _vh) const {
        return position_[_vh];
    }

    void swap_vertices(std::vector<VecT> &_other) {
        position_.swap(_other);
    }

public:

    typename PointT::value_type length(HalfEdgeHandle _heh) const {
        return vector(_heh).length();
    }

    typename PointT::value_type length(EdgeHandle _eh) const {
        return vector(_eh).length();
    }

    PointT vector(HalfEdgeHandle _heh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::halfedge(_heh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex()));
    }

    PointT vector(EdgeHandle _eh) const {

        const typename TopologyKernelT::Edge& e = TopologyKernelT::edge(_eh);
        return (vertex(e.to_vertex()) - vertex(e.from_vertex()));
    }

    PointT barycenter(EdgeHandle _eh) const {
        return PointT(0.5 * vertex(TopologyKernelT::edge(_eh).from_vertex()) +
                      0.5 * vertex(TopologyKernelT::edge(_eh).to_vertex()));
    }

    PointT barycenter(FaceHandle _fh) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        HalfFaceVertexIter hfv_it =
                TopologyKernelT::hfv_iter(TopologyKernelT::halfface_handle(_fh, 0));
        for(; hfv_it.valid(); ++hfv_it, valence += 1) {
            p += vertex(*hfv_it);
        }
        p /= valence;
        return p;
    }

    PointT barycenter(CellHandle _ch) const {
        PointT p(typename PointT::value_type(0));
        typename PointT::value_type valence = 0;
        CellVertexIter cv_it = TopologyKernelT::cv_iter(_ch);
        for(; cv_it.valid(); ++cv_it, valence += 1) {
            p += vertex(*cv_it);
        }
        p /= valence;
        return p;
    }

    /// Compute halfface normal assuming planarity (just uses first 2 edges)
    /// Note: NormalAttrib provides fast access to precomputed normals.
    PointT normal(HalfFaceHandle _hfh) const
    {
        if(TopologyKernelT::halfface(_hfh).halfedges().size() < 3) {
            std::cerr << "Warning: Degenerate face: "
                      << TopologyKernelT::face_handle(_hfh) << std::endl;
            return PointT {0.0, 0.0, 0.0};
        }

        const std::vector<HalfEdgeHandle>& halfedges = TopologyKernelT::halfface(_hfh).halfedges();
        std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();

        const PointT &p1 = vertex(TopologyKernelT::halfedge(*he_it).from_vertex());
        const PointT &p2 = vertex(TopologyKernelT::halfedge(*he_it).to_vertex());
        ++he_it;
        const PointT &p3 = vertex(TopologyKernelT::halfedge(*he_it).to_vertex());

        const PointT n = (p2 - p1).cross(p3 - p2);
        return n.normalized();
    }

    GeometryKernelT<VecT> const& vertex_positions() const & {return position_;}
    GeometryKernelT<VecT> & vertex_positions() & {return position_;}

private:
    GeometryKernelT<VecT> get_prop() {
        auto prop = this->template request_property<VecT, Entity::Vertex>("ovm:position");
        this->set_persistent(prop);
        return prop;
    }

private:
    GeometryKernelT<VecT> position_;

};

} // Namespace OpenVolumeMesh

