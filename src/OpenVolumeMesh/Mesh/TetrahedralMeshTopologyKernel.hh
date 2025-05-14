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


#ifndef NDEBUG
#include <iostream>
#endif
#include <set>

#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshIterators.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {

/**
 * \class TetrahedralMeshTopologyKernel
 *
 * \brief A data structure based on PolyhedralMesh with specializations for tetrahedra.
 *
 */

class OVM_EXPORT TetrahedralMeshTopologyKernel : public TopologyKernel {
public:

    TetrahedralMeshTopologyKernel() = default;
    ~TetrahedralMeshTopologyKernel() override = default;

    FaceHandle add_face(std::vector<HalfEdgeHandle> _halfedges, bool _topologyCheck = false) override;

    FaceHandle add_face(const std::vector<VertexHandle>& _vertices) override;

    CellHandle add_cell(std::vector<HalfFaceHandle> _halffaces, bool _topologyCheck = false) override;

    CellHandle add_cell(const std::vector<VertexHandle>& _vertices, bool _topologyCheck = false);

    CellHandle add_cell(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2, VertexHandle _vh3, bool _topologyCheck = false);

    HalfFaceHandle add_halfface(const std::vector<HalfEdgeHandle>& _halfedges, bool _topologyCheck = false);
    HalfFaceHandle add_halfface(VertexHandle _vh0, VertexHandle _vh1, VertexHandle _vh2, bool _topologyCheck = false);

    HalfEdgeHandle add_halfedge(VertexHandle _fromVertex, VertexHandle _toVertex);

    /// Get the 4 vertices of the tet ch in a specific order:
    /// 1.-3. vertices of ch's first halfface, ccw, starting with the
    /// first from_vertex of the halfface's first halfedge. 4. Then comes the 4th vertex of the tet.
    std::vector<VertexHandle> get_cell_vertices(CellHandle ch) const;

    /// Get the 4 vertices of the tet ch in a specific order, starting with vh.
    std::vector<VertexHandle> get_cell_vertices(CellHandle ch, VertexHandle vh) const;

    /// Get the 4 vertices of hfh's incident cell in a specific order:
    /// 1.-3. vertices of hfh, ccw, starting with the
    /// first from_vertex of the halfface's first halfedge. 4. Then comes the 4th vertex of the tet.
    /// Returns an empty vector, if the incident cell is invalid (hfh is boundary).
    std::vector<VertexHandle> get_cell_vertices(HalfFaceHandle hfh) const;

    /// Get the 4 vertices of hfh's incident cell in a specific order:
    /// 1. heh.to_vertex, 2. heh.to_vertex, 3. 3rd vertex of hfh, 4. 4th vertex of the tet.
    /// heh is expected to be incident to hfh
    std::vector<VertexHandle> get_cell_vertices(HalfFaceHandle hfh, HalfEdgeHandle heh) const;

    /// Get the vertex of the halfface's incident cell that is not contained in the halfface hfh.
    /// If the incident cell is invalid (hfh is boundary), returns the invalid vertex handle
    VertexHandle halfface_opposite_vertex(HalfFaceHandle hfh) const;

    /// Get the first halfface of the tet ch that does not contain the vertex vh.
    HalfFaceHandle vertex_opposite_halfface(CellHandle ch, VertexHandle vh) const;


    VertexHandle collapse_edge(HalfEdgeHandle _heh);
protected:
    void split_edge(HalfEdgeHandle _heh, VertexHandle _vh);
    void split_face(FaceHandle _fh, VertexHandle _vh);

public:


    // ======================= Specialized Iterators =============================

    friend class TetVertexIter;

    typedef class TetVertexIter TetVertexIter;

    /// Returns an iterator to iterate over the four vertices of a tetrahedron in a specific order:
    /// 1.-3. vertices of the tet's first halfface,
    /// starting with the halfface's first halfedge's from_vertex, then 4. the remaining fourth vertex of the tet.
    /// Uses get_cell_vertices(ch).
    TetVertexIter tv_iter(CellHandle _ref_h, int _max_laps = 1) const {
        return TetVertexIter(_ref_h, this, _max_laps);
    }

    std::pair<TetVertexIter,TetVertexIter> tet_vertices(CellHandle _ref_h, int _max_laps = 1) const {
        TetVertexIter begin = tv_iter(_ref_h, _max_laps);
        TetVertexIter end   = make_end_circulator(begin);
        return std::make_pair(begin, end);
    }

private:
   // void replaceHalfFace(CellHandle ch, HalfFaceHandle hf_del, HalfFaceHandle hf_ins);
   // void replaceHalfEdge(HalfFaceHandle hfh, HalfEdgeHandle he_del, HalfEdgeHandle he_ins);

    template <typename PropIterator, typename Handle>
    void swapPropertyElements(PropIterator begin, PropIterator end, Handle source, Handle destination)
    {
        PropIterator p_iter =  begin;
        for (; p_iter != end; ++p_iter)
            (*p_iter)->swap_elements(source, destination);
    }
};

} // Namespace OpenVolumeMesh

