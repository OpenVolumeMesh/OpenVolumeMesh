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
#include <set>
#include <vector>
#include <array>

#include <OpenVolumeMesh/Core/HandleIndexing.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/Iterators.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {

// provide begin() and end() for the iterator pairs provided in TopologyKernel,
// so we can use range-for, e.g. for(const auto &vh: mesh.vertices()) works.
template<class I>
typename std::enable_if<is_ovm_iterator<I>::value, I>::type
begin(const std::pair<I, I>& iterpair)
{
    return iterpair.first;
}

template<class I>
typename std::enable_if<is_ovm_iterator<I>::value, I>::type
end(const std::pair<I, I>& iterpair)
{
    return iterpair.second;
}



class OVM_EXPORT TopologyKernel : public ResourceManager {
public:

    TopologyKernel() = default;
    ~TopologyKernel() override = default;

    TopologyKernel(TopologyKernel const& other) = default;
    TopologyKernel(TopologyKernel      &&other) = delete;

    TopologyKernel& operator=(TopologyKernel const& other) = default;
    TopologyKernel& operator=(TopologyKernel && other) = delete;

    /*
     * Defines and constants
     */

    static const VertexHandle   InvalidVertexHandle;
    static const EdgeHandle     InvalidEdgeHandle;
    static const FaceHandle     InvalidFaceHandle;
    static const CellHandle     InvalidCellHandle;
    static const HalfEdgeHandle InvalidHalfEdgeHandle;
    static const HalfFaceHandle InvalidHalfFaceHandle;

    typedef OpenVolumeMeshEdge Edge;
    typedef OpenVolumeMeshFace Face;
    typedef OpenVolumeMeshCell Cell;

    // Add StatusAttrib to list of friend classes
    // since it provides a garbage collection
    // that needs access to some protected methods
    friend class StatusAttrib;

    //=====================================================================
    // Iterators
    //=====================================================================

    friend class VertexOHalfEdgeIter;
    friend class VertexVertexIter;
    friend class HalfEdgeHalfFaceIter;
    friend class VertexFaceIter;
    friend class VertexCellIter;
    friend class HalfEdgeCellIter;
    friend class CellVertexIter;
    friend class CellCellIter;
    friend class HalfFaceVertexIter;
    friend class BoundaryHalfFaceHalfFaceIter;
    friend class VertexIter;
    friend class EdgeIter;
    friend class HalfEdgeIter;
    friend class FaceIter;
    friend class HalfFaceIter;
    friend class CellIter;

    /*
     * Circulators
     */

protected:

    template <class Circulator>
    static Circulator make_end_circulator(const Circulator& _circ)
    {
        Circulator end = _circ;
        if (end.valid()) {
            end.lap(_circ.max_laps());
            end.valid(false);
        }
        return end;
    }

public:

    VertexVertexIter vv_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexVertexIter(_h, this, _max_laps);
    }

    std::pair<VertexVertexIter, VertexVertexIter> vertex_vertices(VertexHandle _h, int _max_laps = 1) const {
        VertexVertexIter begin = vv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexOHalfEdgeIter voh_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexOHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexOHalfEdgeIter, VertexOHalfEdgeIter> outgoing_halfedges(VertexHandle _h, int _max_laps = 1) const {
        VertexOHalfEdgeIter begin = voh_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexIHalfEdgeIter vih_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexIHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexIHalfEdgeIter, VertexIHalfEdgeIter> incoming_halfedges(VertexHandle _h, int _max_laps = 1) const {
        VertexIHalfEdgeIter begin = vih_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexEdgeIter ve_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexEdgeIter, VertexEdgeIter> vertex_edges(VertexHandle _h, int _max_laps = 1) const {
        VertexEdgeIter begin = ve_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexHalfFaceIter vhf_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<VertexHalfFaceIter, VertexHalfFaceIter> vertex_halffaces(VertexHandle _h, int _max_laps = 1) const {
        VertexHalfFaceIter begin = vhf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexFaceIter vf_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexFaceIter(_h, this, _max_laps);
    }

    std::pair<VertexFaceIter, VertexFaceIter> vertex_faces(VertexHandle _h, int _max_laps = 1) const {
        VertexFaceIter begin = vf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexCellIter vc_iter(VertexHandle _h, int _max_laps = 1) const {
        return VertexCellIter(_h, this, _max_laps);
    }

    std::pair<VertexCellIter, VertexCellIter> vertex_cells(VertexHandle _h, int _max_laps = 1) const {
        VertexCellIter begin = vc_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeHalfFaceIter hehf_iter(HalfEdgeHandle _h, int _max_laps = 1) const {
        return HalfEdgeHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeHalfFaceIter, HalfEdgeHalfFaceIter> halfedge_halffaces(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeHalfFaceIter begin = hehf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeFaceIter hef_iter(HalfEdgeHandle _h, int _max_laps = 1) const {
        return HalfEdgeFaceIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeFaceIter, HalfEdgeFaceIter> halfedge_faces(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeFaceIter begin = hef_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeCellIter hec_iter(HalfEdgeHandle _h, int _max_laps = 1) const {
        return HalfEdgeCellIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeCellIter, HalfEdgeCellIter> halfedge_cells(HalfEdgeHandle _h, int _max_laps = 1) const {
        HalfEdgeCellIter begin = hec_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeHalfFaceIter ehf_iter(EdgeHandle _h, int _max_laps = 1) const {
        return EdgeHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<EdgeHalfFaceIter, EdgeHalfFaceIter> edge_halffaces(EdgeHandle _h, int _max_laps = 1) const {
        EdgeHalfFaceIter begin = ehf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeFaceIter ef_iter(EdgeHandle _h, int _max_laps = 1) const {
        return EdgeFaceIter(_h, this, _max_laps);
    }

    std::pair<EdgeFaceIter, EdgeFaceIter> edge_faces(EdgeHandle _h, int _max_laps = 1) const {
        EdgeFaceIter begin = ef_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeCellIter ec_iter(EdgeHandle _h, int _max_laps = 1) const {
        return EdgeCellIter(_h, this, _max_laps);
    }

    std::pair<EdgeCellIter, EdgeCellIter> edge_cells(EdgeHandle _h, int _max_laps = 1) const {
        EdgeCellIter begin = ec_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceHalfEdgeIter hfhe_iter(HalfFaceHandle _h, int _max_laps = 1) const {
        return HalfFaceHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceHalfEdgeIter, HalfFaceHalfEdgeIter> halfface_halfedges(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceHalfEdgeIter begin = hfhe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceEdgeIter hfe_iter(HalfFaceHandle _h, int _max_laps = 1) const {
        return HalfFaceEdgeIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceEdgeIter, HalfFaceEdgeIter> halfface_edges(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceEdgeIter begin = hfe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceVertexIter fv_iter(FaceHandle _h, int _max_laps = 1) const {
        return FaceVertexIter(_h, this, _max_laps);
    }

    std::pair<FaceVertexIter, FaceVertexIter> face_vertices(FaceHandle _h, int _max_laps = 1) const {
        FaceVertexIter begin = fv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceHalfEdgeIter fhe_iter(FaceHandle _h, int _max_laps = 1) const {
        return FaceHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<FaceHalfEdgeIter, FaceHalfEdgeIter> face_halfedges(FaceHandle _h, int _max_laps = 1) const {
        FaceHalfEdgeIter begin = fhe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceEdgeIter fe_iter(FaceHandle _h, int _max_laps = 1) const {
        return FaceEdgeIter(_h, this, _max_laps);
    }

    std::pair<FaceEdgeIter, FaceEdgeIter> face_edges(FaceHandle _h, int _max_laps = 1) const {
        FaceEdgeIter begin = fe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellVertexIter cv_iter(CellHandle _h, int _max_laps = 1) const {
        return CellVertexIter(_h, this, _max_laps);
    }

    std::pair<CellVertexIter, CellVertexIter> cell_vertices(CellHandle _h, int _max_laps = 1) const {
        CellVertexIter begin = cv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellHalfEdgeIter che_iter(CellHandle _h, int _max_laps = 1) const {
        return CellHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<CellHalfEdgeIter, CellHalfEdgeIter> cell_halfedges(CellHandle _h, int _max_laps = 1) const {
        CellHalfEdgeIter begin = che_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellEdgeIter ce_iter(CellHandle _h, int _max_laps = 1) const {
        return CellEdgeIter(_h, this, _max_laps);
    }

    std::pair<CellEdgeIter, CellEdgeIter> cell_edges(CellHandle _h, int _max_laps = 1) const {
        CellEdgeIter begin = ce_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellHalfFaceIter chf_iter(CellHandle _h, int _max_laps = 1) const {
        return CellHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<CellHalfFaceIter, CellHalfFaceIter> cell_halffaces(CellHandle _h, int _max_laps = 1) const {
        CellHalfFaceIter begin = chf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellFaceIter cf_iter(CellHandle _h, int _max_laps = 1) const {
        return CellFaceIter(_h, this, _max_laps);
    }

    std::pair<CellFaceIter, CellFaceIter> cell_faces(CellHandle _h, int _max_laps = 1) const {
        CellFaceIter begin = cf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellCellIter cc_iter(CellHandle _h, int _max_laps = 1) const {
        return CellCellIter(_h, this, _max_laps);
    }

    std::pair<CellCellIter, CellCellIter> cell_cells(CellHandle _h, int _max_laps = 1) const {
        CellCellIter begin = cc_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceVertexIter hfv_iter(HalfFaceHandle _h, int _max_laps = 1) const {
        return HalfFaceVertexIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceVertexIter, HalfFaceVertexIter> halfface_vertices(HalfFaceHandle _h, int _max_laps = 1) const {
        HalfFaceVertexIter begin = hfv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    BoundaryHalfFaceHalfFaceIter bhfhf_iter(HalfFaceHandle _ref_h, int _max_laps = 1) const {
        return BoundaryHalfFaceHalfFaceIter(_ref_h, this, _max_laps);
    }

    std::pair<BoundaryHalfFaceHalfFaceIter, BoundaryHalfFaceHalfFaceIter> boundary_halfface_halffaces(HalfFaceHandle _h, int _max_laps = 1) const {
        BoundaryHalfFaceHalfFaceIter begin = bhfhf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    /*
     * Iterators
     */

    BoundaryVertexIter bv_iter() const {
        return BoundaryVertexIter(this);
    }

    BoundaryHalfEdgeIter bhe_iter() const {
        return BoundaryHalfEdgeIter(this);
    }

    BoundaryEdgeIter be_iter() const {
        return BoundaryEdgeIter(this);
    }

    BoundaryHalfFaceIter bhf_iter() const {
        return BoundaryHalfFaceIter(this);
    }

    BoundaryFaceIter bf_iter() const {
        return BoundaryFaceIter(this);
    }

    BoundaryCellIter bc_iter() const {
        return BoundaryCellIter(this);
    }

    VertexIter v_iter() const {
        return VertexIter(this);
    }

    VertexIter vertices_begin() const {
        return VertexIter(this, VertexHandle(0));
    }

    VertexIter vertices_end() const {
        return VertexIter(this, VertexHandle((int)n_vertices()));
    }

    std::pair<VertexIter, VertexIter> vertices() const {
        return std::make_pair(vertices_begin(), vertices_end());
    }

    EdgeIter e_iter() const {
        return EdgeIter(this);
    }

    EdgeIter edges_begin() const {
        return EdgeIter(this, EdgeHandle(0));
    }

    EdgeIter edges_end() const {
        return EdgeIter(this, EdgeHandle((int)edges_.size()));
    }

    std::pair<EdgeIter, EdgeIter> edges() const {
        return std::make_pair(edges_begin(), edges_end());
    }

    HalfEdgeIter he_iter() const {
        return HalfEdgeIter(this);
    }

    HalfEdgeIter halfedges_begin() const {
        return HalfEdgeIter(this, HalfEdgeHandle(0));
    }

    HalfEdgeIter halfedges_end() const {
        return HalfEdgeIter(this, HalfEdgeHandle((int)edges_.size() * 2));
    }

    std::pair<HalfEdgeIter, HalfEdgeIter> halfedges() const {
        return std::make_pair(halfedges_begin(), halfedges_end());
    }

    FaceIter f_iter() const {
        return FaceIter(this);
    }

    FaceIter faces_begin() const {
        return FaceIter(this, FaceHandle(0));
    }

    FaceIter faces_end() const {
        return FaceIter(this, FaceHandle((int)faces_.size()));
    }

    std::pair<FaceIter, FaceIter> faces() const {
        return std::make_pair(faces_begin(), faces_end());
    }

    HalfFaceIter hf_iter() const {
        return HalfFaceIter(this);
    }

    HalfFaceIter halffaces_begin() const {
        return HalfFaceIter(this, HalfFaceHandle(0));
    }

    HalfFaceIter halffaces_end() const {
        return HalfFaceIter(this, HalfFaceHandle((int)faces_.size() * 2));
    }

    std::pair<HalfFaceIter, HalfFaceIter> halffaces() const {
        return std::make_pair(halffaces_begin(), halffaces_end());
    }

    CellIter c_iter() const {
        return CellIter(this);
    }

    CellIter cells_begin() const {
        return CellIter(this, CellHandle(0));
    }

    CellIter cells_end() const {
        return CellIter(this, CellHandle((int)cells_.size()));
    }

    std::pair<CellIter, CellIter> cells() const {
        return std::make_pair(cells_begin(), cells_end());
    }

    /*
     * Convenience functions
     */

    std::array<VertexHandle, 2> halfedge_vertices(HalfEdgeHandle _h) const {
        return {from_vertex_handle(_h), to_vertex_handle(_h)};
    }

    std::array<VertexHandle, 2> edge_vertices(EdgeHandle _h) const {
        return halfedge_vertices(halfedge_handle(_h, 0));
    }

    std::array<HalfEdgeHandle, 2> edge_halfedges(EdgeHandle _h) const {
        return {
            halfedge_handle(_h, 0),
            halfedge_handle(_h, 1)};
    }

    std::array<HalfFaceHandle,2> face_halffaces(FaceHandle _h) const {
        return {
            halfface_handle(_h, 0),
            halfface_handle(_h, 1)};
    }

    std::array<CellHandle, 2> face_cells(FaceHandle _h) const {
        return {
        incident_cell(halfface_handle(_h, 0)),
        incident_cell(halfface_handle(_h, 1))};
    }

    /*
     * Virtual functions with implementation
     */

    /// Get number of vertices in mesh
    size_t n_vertices()   const override { return n_vertices_; }
    /// Get number of edges in mesh
    size_t n_edges()      const override { return edges_.size(); }
    /// Get number of halfedges in mesh
    size_t n_halfedges()  const override { return edges_.size() * 2u; }
    /// Get number of faces in mesh
    size_t n_faces()      const override { return faces_.size(); }
    /// Get number of halffaces in mesh
    size_t n_halffaces()  const override { return faces_.size() * 2u; }
    /// Get number of cells in mesh
    size_t n_cells()      const override { return cells_.size(); }

    /// Get number of undeleted vertices in mesh
    size_t n_logical_vertices()   const { return n_vertices_ - n_deleted_vertices_; }
    /// Get number of undeleted edges in mesh
    size_t n_logical_edges()      const { return edges_.size() - n_deleted_edges_; }
    /// Get number of undeleted halfedges in mesh
    size_t n_logical_halfedges()  const { return n_logical_edges() * 2u; }
    /// Get number of undeleted faces in mesh
    size_t n_logical_faces()      const { return faces_.size() - n_deleted_faces_; }
    /// Get number of undeleted halffaces in mesh
    size_t n_logical_halffaces()  const { return n_faces() * 2u; }
    /// Get number of undeleted cells in mesh
    size_t n_logical_cells()      const { return cells_.size() - n_deleted_cells_; }

    int genus() const {

        int g = (1 - (int)(n_vertices() -
                           n_edges() +
                           n_faces() -
                           n_cells()));

        if(g % 2 == 0) return (g / 2);

        // An error occured
        // The mesh might not be manifold
        return  -1;
    }

private:

    // Cache total vertex number
    size_t n_vertices_ = 0u;

public:
    void add_n_vertices(size_t n);
    void reserve_vertices(size_t n);
    void reserve_edges(size_t n);
    void reserve_faces(size_t n);
    void reserve_cells(size_t n);

    /// Add abstract vertex
    virtual VertexHandle add_vertex();

    //=======================================================================

    /// Add edge
    virtual EdgeHandle add_edge(VertexHandle _fromVertex, VertexHandle _toHandle, bool _allowDuplicates = false);

    /// Add face via incident edges
    ///
    /// \return Handle of the new face, InvalidFaceHandle if \a _halfedges
    ///         are not connected and \a _topologyCheck is \a true.
    ///
    /// \warning If _halfedges are not connected or not in the right order and \a _topologyCheck is \a false,
    ///          the behavior is undefined.
    virtual FaceHandle add_face(std::vector<HalfEdgeHandle> _halfedges, bool _topologyCheck = false);

    /// Add face via incident vertices
    virtual FaceHandle add_face(const std::vector<VertexHandle>& _vertices);

    /// Add cell via incident halffaces
    ///
    /// \return Handle of the new cell, InvalidCellHandle if \a _topologyCheck is \a true and
    ///         \a _halffaces are not connected.
    ///
    /// \warning If _halffaces are not connected and \a _topologyCheck is \a false,
    ///          the behavior is undefined.
    virtual CellHandle add_cell(std::vector<HalfFaceHandle> _halffaces, bool _topologyCheck = false);

    /// Set the vertices of an edge
    void set_edge(EdgeHandle _eh, VertexHandle _fromVertex, VertexHandle _toVertex);

    /// Set the half-edges of a face
    void set_face(FaceHandle _fh, const std::vector<HalfEdgeHandle>& _hes);

    /// Set the half-faces of a cell
    void set_cell(CellHandle _ch, const std::vector<HalfFaceHandle>& _hfs);

    /// Recompute cyclic ordering of (half)faces incident to an edge (used by iterators)
    void reorder_incident_halffaces(EdgeHandle _eh);

    /*
     * Non-virtual functions
     */

    /// Get edge with handle _edgeHandle
    const Edge& edge(EdgeHandle _edgeHandle) const;

    /// Get face with handle _faceHandle
    const Face& face(FaceHandle _faceHandle) const;

    /// Get cell with handle _cellHandle
    const Cell& cell(CellHandle _cellHandle) const;

    /// Get edge with handle _edgeHandle
    Edge& edge(EdgeHandle _edgeHandle);

    /// Get face with handle _faceHandle
    Face& face(FaceHandle _faceHandle);

    /// Get cell with handle _cellHandle
    Cell& cell(CellHandle _cellHandle);

    /// Get edge that corresponds to halfedge with handle _halfEdgeHandle
    Edge halfedge(HalfEdgeHandle _halfEdgeHandle) const;

    /// Get face that corresponds to halfface with handle _halfFaceHandle
    Face halfface(HalfFaceHandle _halfFaceHandle) const;

    /// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
    Edge opposite_halfedge(HalfEdgeHandle _halfEdgeHandle) const;

    /// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
    Face opposite_halfface(HalfFaceHandle _halfFaceHandle) const;

    /// Get halfedge from vertex _vh1 to _vh2
    HalfEdgeHandle find_halfedge(VertexHandle _vh1, VertexHandle _vh2) const;
    [[deprecated("please use find_halfedge instead")]]
    HalfEdgeHandle halfedge     (VertexHandle _vh1, VertexHandle _vh2) const; // deprecated, use the one above instead!!!

    /// Get halfedge from vertex _vh1 to _vh2 but restricted to halfedges of cell _ch
    HalfEdgeHandle find_halfedge_in_cell(VertexHandle _vh1, VertexHandle _vh2, CellHandle _ch) const;

    /// Get half-face from list of incident vertices (in connected order)
    ///
    /// \note Only the first three vertices are checked
    HalfFaceHandle find_halfface(const std::vector<VertexHandle>& _vs) const;
    [[deprecated("please use find_halfface instead")]]
    HalfFaceHandle halfface     (const std::vector<VertexHandle>& _vs) const; // deprecated, use the one above instead!!!

    /// Get half-face from list of incident vertices (in connected order)
    ///
    /// \note Only the first three vertices are checked
    HalfFaceHandle find_halfface_in_cell(const std::vector<VertexHandle>& _vs, CellHandle _ch) const;

  /// Get half-face from list of incident vertices (in connected order)
    ///
    /// \note All vertices are checked
    HalfFaceHandle find_halfface_extensive(const std::vector<VertexHandle>& _vs) const;
    [[deprecated("please use find_halfface_extensive instead")]]
    HalfFaceHandle halfface_extensive     (const std::vector<VertexHandle>& _vs) const; // deprecated, use the one above instead!!!

    /// Get half-face from list of incident half-edges
    ///
    /// \note Only the first two half-edges are checked
    HalfFaceHandle find_halfface(const std::vector<HalfEdgeHandle>& _hes) const;
    [[deprecated("please use find_halfface instead")]]
    HalfFaceHandle halfface     (const std::vector<HalfEdgeHandle>& _hes) const; // deprecated, use the one above instead!!!

    /// Get next halfedge within a halfface
    HalfEdgeHandle next_halfedge_in_halfface(HalfEdgeHandle _heh, HalfFaceHandle _hfh) const;

    /// Get previous halfedge within a halfface
    HalfEdgeHandle prev_halfedge_in_halfface(HalfEdgeHandle _heh, HalfFaceHandle _hfh) const;

    /// Get the vertex the halfedge starts from
    VertexHandle from_vertex_handle(HalfEdgeHandle _h) const {
        return halfedge(_h).from_vertex();
    }

    /// Get the vertex the halfedge points to
    VertexHandle to_vertex_handle(HalfEdgeHandle _h) const {
        return halfedge(_h).to_vertex();
    }

    /// Get valence of vertex (number of incident edges)
    inline size_t valence(VertexHandle _vh) const {
        assert(is_valid(_vh));
        assert(has_vertex_bottom_up_incidences());

        return outgoing_hes_per_vertex_[_vh].size();
    }

    /// Get valence of edge (number of incident faces)
    inline size_t valence(EdgeHandle _eh) const {
        assert(is_valid(_eh));
        assert(has_edge_bottom_up_incidences());

        return incident_hfs_per_he_[halfedge_handle(_eh, 0)].size();
    }

    /// Get valence of face (number of incident edges)
    inline size_t valence(FaceHandle _fh) const {
        assert(is_valid(_fh));

        return face(_fh).halfedges().size();
    }

    /// Get valence of cell (number of incident faces)
    inline size_t valence(CellHandle _ch) const {
        assert(is_valid(_ch));

        return cell(_ch).halffaces().size();
    }

    /// Get vertices of a halfface
    std::vector<VertexHandle> get_halfface_vertices(HalfFaceHandle hfh) const;
    /// Get vertices of a halfface orderd to start from vh
    std::vector<VertexHandle> get_halfface_vertices(HalfFaceHandle hfh, VertexHandle vh) const;
    /// Get vertices of a halfface orderd to start from from_vertex_handle(heh)
    std::vector<VertexHandle> get_halfface_vertices(HalfFaceHandle hfh, HalfEdgeHandle heh) const;

    /// check whether face _fh and edge _eh are incident
    bool is_incident( FaceHandle _fh, EdgeHandle _eh) const;

    //=====================================================================
    // Delete entities
    //=====================================================================

public:

    virtual VertexIter delete_vertex(VertexHandle _h);

    virtual EdgeIter delete_edge(EdgeHandle _h);

    virtual FaceIter delete_face(FaceHandle _h);

    virtual CellIter delete_cell(CellHandle _h);

    virtual void collect_garbage();


    virtual bool is_deleted(VertexHandle _h)   const { return vertex_deleted_[_h]; }
    virtual bool is_deleted(EdgeHandle _h)     const { return edge_deleted_[_h]; }
    virtual bool is_deleted(HalfEdgeHandle _h) const { return edge_deleted_[_h.edge_handle()]; }
    virtual bool is_deleted(FaceHandle _h)     const { return face_deleted_[_h]; }
    virtual bool is_deleted(HalfFaceHandle _h) const { return face_deleted_[_h.face_handle()]; }
    virtual bool is_deleted(CellHandle _h)     const { return cell_deleted_[_h]; }

private:

    template <class ContainerT>
    void get_incident_edges(const ContainerT& _vs, std::set<EdgeHandle>& _es) const;

    template <class ContainerT>
    void get_incident_faces(const ContainerT& _es, std::set<FaceHandle>& _fs) const;

    template <class ContainerT>
    void get_incident_cells(const ContainerT& _fs, std::set<CellHandle>& _cs) const;

    VertexIter delete_vertex_core(VertexHandle _h);

    EdgeIter delete_edge_core(EdgeHandle _h);

    FaceIter delete_face_core(FaceHandle _h);

    CellIter delete_cell_core(CellHandle _h);

public:

    /// Exchanges the indices of two cells while keeping the mesh otherwise unaffected.
    virtual void swap_cell_indices(CellHandle _h1, CellHandle _h2);

    /// Exchanges the indices of two faces while keeping the mesh otherwise unaffected.
    virtual void swap_face_indices(FaceHandle _h1, FaceHandle _h2);

    /// Exchanges the indices of two edges while keeping the mesh otherwise unaffected.
    virtual void swap_edge_indices(EdgeHandle _h1, EdgeHandle _h2);

    /// Exchanges the indices of two vertices while keeping the mesh otherwise unaffected.
    virtual void swap_vertex_indices(VertexHandle _h1, VertexHandle _h2);

protected:

    class EdgeCorrector {
    public:
        explicit EdgeCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Edge& _edge) {
            _edge.set_from_vertex(VertexHandle(newIndices_[_edge.from_vertex().uidx()]));
            _edge.set_to_vertex(VertexHandle(newIndices_[_edge.to_vertex().uidx()]));
        }
    private:
        const std::vector<int>& newIndices_;
    };

    class FaceCorrector {
    public:
        explicit FaceCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Face& _face) {
            std::vector<HalfEdgeHandle> hes = _face.halfedges();
            for(std::vector<HalfEdgeHandle>::iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {

                EdgeHandle eh = edge_handle(*he_it);
                unsigned char opp = he_it->idx() == halfedge_handle(eh, 1).idx();
                *he_it = halfedge_handle(EdgeHandle(newIndices_[eh.uidx()]), opp);
            }
            _face.set_halfedges(hes);
        }
    private:
        const std::vector<int>& newIndices_;
    };

    class CellCorrector {
    public:
        explicit CellCorrector(const std::vector<int>& _newIndices) :
            newIndices_(_newIndices) {}

        void operator()(Cell& _cell) {
            std::vector<HalfFaceHandle> hfs = _cell.halffaces();
            for(std::vector<HalfFaceHandle>::iterator hf_it = hfs.begin(),
                    hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

                FaceHandle fh = face_handle(*hf_it);
                unsigned char opp = hf_it->idx() == halfface_handle(fh, 1).idx();
                *hf_it = halfface_handle(FaceHandle(newIndices_[fh.uidx()]), opp);
            }
            _cell.set_halffaces(hfs);
        }
    private:
        const std::vector<int>& newIndices_;
    };

public:

    /// Clear whole mesh
    virtual void clear(bool _clearProps = true) {

        edges_.clear();
        faces_.clear();
        cells_.clear();
        vertex_deleted_.clear();
        edge_deleted_.clear();
        face_deleted_.clear();
        cell_deleted_.clear();
        n_deleted_vertices_ = 0;
        n_deleted_edges_ = 0;
        n_deleted_faces_ = 0;
        n_deleted_cells_ = 0;
        outgoing_hes_per_vertex_.clear();
        incident_hfs_per_he_.clear();
        incident_cell_per_hf_.clear();
        n_vertices_ = 0;

        if(_clearProps) {
            clear_all_props();
        } else {
            // Resize props
            resize_vprops(0u);
            resize_eprops(0u);
            resize_fprops(0u);
            resize_cprops(0u);
        }
    }

    //=====================================================================
    // Bottom-up Incidences
    //=====================================================================

public:

    void enable_bottom_up_incidences(bool _enable = true)
    {
        enable_vertex_bottom_up_incidences(_enable);
        enable_edge_bottom_up_incidences(_enable);
        enable_face_bottom_up_incidences(_enable);
    }

    void enable_vertex_bottom_up_incidences(bool _enable = true) {

        if(_enable && !v_bottom_up_) {
            // Vertex bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_vertex_bottom_up_incidences();
        }

        if(!_enable) {
            outgoing_hes_per_vertex_.clear();
        }

        v_bottom_up_ = _enable;
    }

    void enable_edge_bottom_up_incidences(bool _enable = true) {

        if(_enable && !e_bottom_up_) {
            // Edge bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_edge_bottom_up_incidences();

            if(f_bottom_up_) {
                for (const auto &eh: edges()) {
                    reorder_incident_halffaces(eh);
                }
            }
        }

        if(!_enable) {
            incident_hfs_per_he_.clear();
        }

        e_bottom_up_ = _enable;
    }

    void enable_face_bottom_up_incidences(bool _enable = true) {

        bool updateOrder = false;
        if(_enable && !f_bottom_up_) {
            // Face bottom-up incidences have to be
            // recomputed for the whole mesh
            compute_face_bottom_up_incidences();

            updateOrder = true;
        }

        if(!_enable) {
            incident_cell_per_hf_.clear();
        }

        f_bottom_up_ = _enable;

        if(updateOrder) {
            if(e_bottom_up_) {
                for (const auto &eh: edges()) {
                    reorder_incident_halffaces(eh);
                }
            }
        }
    }

    bool has_full_bottom_up_incidences() const {
        return (has_vertex_bottom_up_incidences() &&
                has_edge_bottom_up_incidences() &&
                has_face_bottom_up_incidences());
    }

    bool has_vertex_bottom_up_incidences() const { return v_bottom_up_; }

    bool has_edge_bottom_up_incidences() const { return e_bottom_up_; }

    bool has_face_bottom_up_incidences() const { return f_bottom_up_; }


    void enable_deferred_deletion(bool _enable = true);
    bool deferred_deletion_enabled() const { return deferred_deletion_; }


    void enable_fast_deletion(bool _enable = true) { fast_deletion_ = _enable; }
    bool fast_deletion_enabled() const { return fast_deletion_; }


protected:

    void compute_vertex_bottom_up_incidences();

    void compute_edge_bottom_up_incidences();

    void compute_face_bottom_up_incidences();

    // Outgoing halfedges per vertex
    VertexVector<std::vector<HalfEdgeHandle> > outgoing_hes_per_vertex_;

    // Incident halffaces per (directed) halfedge
    HalfEdgeVector<std::vector<HalfFaceHandle> > incident_hfs_per_he_;

    // Incident cell (at most one) per halfface
    HalfFaceVector<CellHandle> incident_cell_per_hf_;

private:
    bool v_bottom_up_ = true;

    bool e_bottom_up_ = true;

    bool f_bottom_up_ = true;

    bool deferred_deletion_ = true;

    bool fast_deletion_ = true;

    //=====================================================================
    // Connectivity
    //=====================================================================

public:

    /// \brief Get halfface that is adjacent (w.r.t. a common halfedge) within the same cell.
    ///        It correctly handles self-adjacent cells where the halfedge orientation matters.
    ///        For legacy reasons the halfedge orientation can be arbitrary if there are no self-adjacencies in the cell.
    ///
    /// \return Handle of the adjacent half-face if \a _halfFaceHandle is not
    ///         at a boundary, \a InvalidHalfFaceHandle otherwise.
    ///
    /// \warning The mesh must have face bottom-up incidences.

    // TODO: We might also want a "common halfedge" API that may assume he \in hf once the deprecated function is eliminated for good
    HalfFaceHandle adjacent_halfface_in_cell(HalfFaceHandle _halfFaceHandle, HalfEdgeHandle _halfEdgeHandle) const;

    /// Get cell that is incident to the given halfface
    CellHandle incident_cell(HalfFaceHandle _halfFaceHandle) const;

    bool is_boundary(HalfFaceHandle _halfFaceHandle) const
    {
        assert(is_valid(_halfFaceHandle));
        assert(has_face_bottom_up_incidences());
        return incident_cell_per_hf_[_halfFaceHandle] == InvalidCellHandle;
    }

    bool is_boundary(FaceHandle _faceHandle) const {
        assert(is_valid(_faceHandle));
        assert(has_face_bottom_up_incidences());
        return  is_boundary(halfface_handle(_faceHandle, 0)) ||
                is_boundary(halfface_handle(_faceHandle, 1));
    }

    bool is_boundary(EdgeHandle _edgeHandle) const {
        assert(is_valid(_edgeHandle));
        assert(has_edge_bottom_up_incidences());

        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(_edgeHandle, 0));
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(HalfEdgeHandle _halfedgeHandle) const {
        assert(is_valid(_halfedgeHandle));
        assert(has_edge_bottom_up_incidences());

        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(_halfedgeHandle);
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(VertexHandle _vertexHandle) const {
        assert(is_valid(_vertexHandle));
        assert(has_vertex_bottom_up_incidences());

        for(VertexOHalfEdgeIter voh_it = voh_iter(_vertexHandle); voh_it.valid(); ++voh_it) {
            if(is_boundary(*voh_it)) return true;
        }
        return false;
    }

    bool is_boundary(CellHandle _cellHandle) const {
        assert(is_valid(_cellHandle));

        for(CellFaceIter cf_it = cf_iter(_cellHandle); cf_it.valid(); ++cf_it) {
            if(is_boundary(*cf_it)) return true;
        }
        return false;
    }

    size_t n_vertices_in_cell(CellHandle _ch) const {
        assert(is_valid(_ch));

        std::set<VertexHandle> vhs;
        for(const auto hfh: cell_halffaces(_ch)) {
            for(const auto heh: halfface_halfedges(hfh)) {
                vhs.insert(to_vertex_handle(heh));
            }
        }
        return vhs.size();
    }

    //=========================================================================

    /*
     * Non-virtual functions
     */

    Edge opposite_halfedge(const Edge& _edge) const {
        return Edge(_edge.to_vertex(), _edge.from_vertex());
    }

    Face opposite_halfface(const Face& _face) const {
        std::vector<HalfEdgeHandle> opp_halfedges;
        opp_halfedges.resize(_face.halfedges().size());
        std::transform(_face.halfedges().rbegin(),
                       _face.halfedges().rend(),
                       opp_halfedges.begin(),
                       opposite_halfedge_handle);
        return Face(opp_halfedges);
    }

    /*
     * Static functions
     */

    /// Conversion function
    static inline HalfEdgeHandle halfedge_handle(EdgeHandle _h, const unsigned char _subIdx) {
        // Is handle in range?
        assert(_h.is_valid());
        assert(_subIdx < 2);
        // if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfEdgeHandle;
        return HalfEdgeHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Conversion function
    static inline HalfFaceHandle halfface_handle(FaceHandle _h, const unsigned char _subIdx) {
        // Is handle in range?
        assert(_h.is_valid());
        assert(_subIdx < 2);
        // if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfFaceHandle;
        return HalfFaceHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Handle conversion
    static inline EdgeHandle edge_handle(HalfEdgeHandle _h) {
        assert(_h.is_valid());
        return EdgeHandle(_h.idx() / 2);
    }

    static inline FaceHandle face_handle(HalfFaceHandle _h) {
        assert(_h.is_valid());
        return FaceHandle(_h.idx() / 2);
    }

    static inline HalfEdgeHandle opposite_halfedge_handle(HalfEdgeHandle _h) {
        assert(_h.is_valid());
        return HalfEdgeHandle(_h.idx() ^ 1);
    }

    static inline HalfFaceHandle opposite_halfface_handle(HalfFaceHandle _h) {
        assert(_h.is_valid());
        return HalfFaceHandle(_h.idx() ^ 1);
    }

    bool inline needs_garbage_collection() const {
        return n_deleted_vertices_ > 0 || n_deleted_edges_ > 0 || n_deleted_faces_ > 0 || n_deleted_cells_ > 0;
    }

    /// test is_valid and perform index range check
    template<typename Handle>
    bool is_valid(Handle _h) const {
        static_assert(is_handle_v<Handle>);
        return _h.is_valid() && _h.uidx() < n<typename Handle::EntityTag>();
    }

private:

    // top-down incidences:
    EdgeVector<Edge> edges_;
    FaceVector<Face> faces_;
    CellVector<Cell> cells_;

    // deferred deletion:
    VertexVector<bool> vertex_deleted_;
    EdgeVector<bool> edge_deleted_;
    FaceVector<bool> face_deleted_;
    CellVector<bool> cell_deleted_;

    // number of elements deleted, but not yet garbage collected
    size_t n_deleted_vertices_ = 0;
    size_t n_deleted_edges_ = 0;
    size_t n_deleted_faces_ = 0;
    size_t n_deleted_cells_ = 0;
};

}

