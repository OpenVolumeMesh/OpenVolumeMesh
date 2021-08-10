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

#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/Iterators.hh>
#include <OpenVolumeMesh/Core/detail/VertexHalfEdgeIncidence.hh>
#include <OpenVolumeMesh/Core/detail/HalfEdgeHalfFaceIncidence.hh>
#include <OpenVolumeMesh/Core/detail/HalfFaceCellIncidence.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {

class OVM_EXPORT TopologyKernel
        : public ResourceManager
        , public VertexHalfEdgeIncidence<TopologyKernel>
        , public HalfEdgeHalfFaceIncidence<TopologyKernel>
        , public HalfFaceCellIncidence<TopologyKernel>
{
public:

    TopologyKernel();
    ~TopologyKernel() override = default;

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

    VertexVertexIter vv_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexVertexIter(_h, this, _max_laps);
    }

    std::pair<VertexVertexIter, VertexVertexIter> vertex_vertices(const VertexHandle& _h, int _max_laps = 1) const {
        VertexVertexIter begin = vv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexOHalfEdgeIter voh_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexOHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexOHalfEdgeIter, VertexOHalfEdgeIter> outgoing_halfedges(const VertexHandle& _h, int _max_laps = 1) const {
        VertexOHalfEdgeIter begin = voh_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexIHalfEdgeIter vih_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexIHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexIHalfEdgeIter, VertexIHalfEdgeIter> incoming_halfedges(const VertexHandle& _h, int _max_laps = 1) const {
        VertexIHalfEdgeIter begin = vih_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexEdgeIter ve_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexEdgeIter(_h, this, _max_laps);
    }

    std::pair<VertexEdgeIter, VertexEdgeIter> vertex_edges(const VertexHandle& _h, int _max_laps = 1) const {
        VertexEdgeIter begin = ve_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexHalfFaceIter vhf_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<VertexHalfFaceIter, VertexHalfFaceIter> vertex_halffaces(const VertexHandle& _h, int _max_laps = 1) const {
        VertexHalfFaceIter begin = vhf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexFaceIter vf_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexFaceIter(_h, this, _max_laps);
    }

    std::pair<VertexFaceIter, VertexFaceIter> vertex_faces(const VertexHandle& _h, int _max_laps = 1) const {
        VertexFaceIter begin = vf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    VertexCellIter vc_iter(const VertexHandle& _h, int _max_laps = 1) const {
        return VertexCellIter(_h, this, _max_laps);
    }

    std::pair<VertexCellIter, VertexCellIter> vertex_cells(const VertexHandle& _h, int _max_laps = 1) const {
        VertexCellIter begin = vc_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeHalfFaceIter hehf_iter(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        return HalfEdgeHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeHalfFaceIter, HalfEdgeHalfFaceIter> halfedge_halffaces(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        HalfEdgeHalfFaceIter begin = hehf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeFaceIter hef_iter(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        return HalfEdgeFaceIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeFaceIter, HalfEdgeFaceIter> halfedge_faces(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        HalfEdgeFaceIter begin = hef_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfEdgeCellIter hec_iter(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        return HalfEdgeCellIter(_h, this, _max_laps);
    }

    std::pair<HalfEdgeCellIter, HalfEdgeCellIter> halfedge_cells(const HalfEdgeHandle& _h, int _max_laps = 1) const {
        HalfEdgeCellIter begin = hec_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeHalfFaceIter ehf_iter(const EdgeHandle& _h, int _max_laps = 1) const {
        return EdgeHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<EdgeHalfFaceIter, EdgeHalfFaceIter> edge_halffaces(const EdgeHandle& _h, int _max_laps = 1) const {
        EdgeHalfFaceIter begin = ehf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeFaceIter ef_iter(const EdgeHandle& _h, int _max_laps = 1) const {
        return EdgeFaceIter(_h, this, _max_laps);
    }

    std::pair<EdgeFaceIter, EdgeFaceIter> edge_faces(const EdgeHandle& _h, int _max_laps = 1) const {
        EdgeFaceIter begin = ef_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    EdgeCellIter ec_iter(const EdgeHandle& _h, int _max_laps = 1) const {
        return EdgeCellIter(_h, this, _max_laps);
    }

    std::pair<EdgeCellIter, EdgeCellIter> edge_cells(const EdgeHandle& _h, int _max_laps = 1) const {
        EdgeCellIter begin = ec_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceHalfEdgeIter hfhe_iter(const HalfFaceHandle& _h, int _max_laps = 1) const {
        return HalfFaceHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceHalfEdgeIter, HalfFaceHalfEdgeIter> halfface_halfedges(const HalfFaceHandle& _h, int _max_laps = 1) const {
        HalfFaceHalfEdgeIter begin = hfhe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceEdgeIter hfe_iter(const HalfFaceHandle& _h, int _max_laps = 1) const {
        return HalfFaceEdgeIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceEdgeIter, HalfFaceEdgeIter> halfface_edges(const HalfFaceHandle& _h, int _max_laps = 1) const {
        HalfFaceEdgeIter begin = hfe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceVertexIter fv_iter(const FaceHandle& _h, int _max_laps = 1) const {
        return FaceVertexIter(_h, this, _max_laps);
    }

    std::pair<FaceVertexIter, FaceVertexIter> face_vertices(const FaceHandle& _h, int _max_laps = 1) const {
        FaceVertexIter begin = fv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceHalfEdgeIter fhe_iter(const FaceHandle& _h, int _max_laps = 1) const {
        return FaceHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<FaceHalfEdgeIter, FaceHalfEdgeIter> face_halfedges(const FaceHandle& _h, int _max_laps = 1) const {
        FaceHalfEdgeIter begin = fhe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    FaceEdgeIter fe_iter(const FaceHandle& _h, int _max_laps = 1) const {
        return FaceEdgeIter(_h, this, _max_laps);
    }

    std::pair<FaceEdgeIter, FaceEdgeIter> face_edges(const FaceHandle& _h, int _max_laps = 1) const {
        FaceEdgeIter begin = fe_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellVertexIter cv_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellVertexIter(_h, this, _max_laps);
    }

    std::pair<CellVertexIter, CellVertexIter> cell_vertices(const CellHandle& _h, int _max_laps = 1) const {
        CellVertexIter begin = cv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellHalfEdgeIter che_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellHalfEdgeIter(_h, this, _max_laps);
    }

    std::pair<CellHalfEdgeIter, CellHalfEdgeIter> cell_halfedges(const CellHandle& _h, int _max_laps = 1) const {
        CellHalfEdgeIter begin = che_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellEdgeIter ce_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellEdgeIter(_h, this, _max_laps);
    }

    std::pair<CellEdgeIter, CellEdgeIter> cell_edges(const CellHandle& _h, int _max_laps = 1) const {
        CellEdgeIter begin = ce_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellHalfFaceIter chf_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellHalfFaceIter(_h, this, _max_laps);
    }

    std::pair<CellHalfFaceIter, CellHalfFaceIter> cell_halffaces(const CellHandle& _h, int _max_laps = 1) const {
        CellHalfFaceIter begin = chf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellFaceIter cf_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellFaceIter(_h, this, _max_laps);
    }

    std::pair<CellFaceIter, CellFaceIter> cell_faces(const CellHandle& _h, int _max_laps = 1) const {
        CellFaceIter begin = cf_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    CellCellIter cc_iter(const CellHandle& _h, int _max_laps = 1) const {
        return CellCellIter(_h, this, _max_laps);
    }

    std::pair<CellCellIter, CellCellIter> cell_cells(const CellHandle& _h, int _max_laps = 1) const {
        CellCellIter begin = cc_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    HalfFaceVertexIter hfv_iter(const HalfFaceHandle& _h, int _max_laps = 1) const {
        return HalfFaceVertexIter(_h, this, _max_laps);
    }

    std::pair<HalfFaceVertexIter, HalfFaceVertexIter> halfface_vertices(const HalfFaceHandle& _h, int _max_laps = 1) const {
        HalfFaceVertexIter begin = hfv_iter(_h, _max_laps);
        return std::make_pair(begin, make_end_circulator(begin));
    }

    BoundaryHalfFaceHalfFaceIter bhfhf_iter(const HalfFaceHandle& _ref_h, int _max_laps = 1) const {
        return BoundaryHalfFaceHalfFaceIter(_ref_h, this, _max_laps);
    }

    std::pair<BoundaryHalfFaceHalfFaceIter, BoundaryHalfFaceHalfFaceIter> boundary_halfface_halffaces(const HalfFaceHandle& _h, int _max_laps = 1) const {
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

    std::vector<VertexHandle> halfedge_vertices(const HalfEdgeHandle& _h) const {
        std::vector<VertexHandle> res(2);
        res[0] = from_vertex_handle(_h);
        res[1] = to_vertex_handle(_h);
        return res;
    }

    std::vector<VertexHandle> edge_vertices(const EdgeHandle& _h) const {
        return halfedge_vertices(halfedge_handle(_h, 0));
    }

    std::vector<HalfEdgeHandle> edge_halfedges(const EdgeHandle& _h) const {
        std::vector<HalfEdgeHandle> res(2);
        res[0] = halfedge_handle(_h, 0);
        res[1] = halfedge_handle(_h, 1);
        return res;
    }

    std::vector<HalfFaceHandle> face_halffaces(const FaceHandle& _h) const {
        std::vector<HalfFaceHandle> res(2);
        res[0] = halfface_handle(_h, 0);
        res[1] = halfface_handle(_h, 1);
        return res;
    }

    std::vector<CellHandle> face_cells(const FaceHandle& _h) const {
        std::vector<CellHandle> res(2);
        res[0] = incident_cell(halfface_handle(_h, 0));
        res[1] = incident_cell(halfface_handle(_h, 1));
        return res;
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
    void reserve_vertices(size_t n);
    void reserve_edges(size_t n);
    void reserve_faces(size_t n);
    void reserve_cells(size_t n);

    /// Add abstract vertex
    virtual VertexHandle add_vertex();

    //=======================================================================

    /// Add edge
    virtual EdgeHandle add_edge(const VertexHandle& _fromVertex, const VertexHandle& _toHandle, bool _allowDuplicates = false);

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
    void set_edge(const EdgeHandle& _eh, const VertexHandle& _fromVertex, const VertexHandle& _toVertex);

    /// Set the half-edges of a face
    void set_face(const FaceHandle& _fh, const std::vector<HalfEdgeHandle>& _hes);

    /// Set the half-faces of a cell
    void set_cell(const CellHandle& _ch, const std::vector<HalfFaceHandle>& _hfs);

    /*
     * Non-virtual functions
     */

    /// Get edge with handle _edgeHandle
    const Edge& edge(const EdgeHandle& _edgeHandle) const;

    /// Get face with handle _faceHandle
    const Face& face(const FaceHandle& _faceHandle) const;

    /// Get cell with handle _cellHandle
    const Cell& cell(const CellHandle& _cellHandle) const;

protected:
    /// Get edge with handle _edgeHandle
    Edge& edge_mutable(const EdgeHandle& _edgeHandle);

    /// Get face with handle _faceHandle
    Face& face_mutable(const FaceHandle& _faceHandle);

    /// Get cell with handle _cellHandle
    Cell& cell_mutable(const CellHandle& _cellHandle);

public:

    /// Get edge that corresponds to halfedge with handle _halfEdgeHandle
    Edge halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get face that corresponds to halfface with handle _halfFaceHandle
    Face halfface(const HalfFaceHandle& _halfFaceHandle) const;

    /// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
    Edge opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
    Face opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const;

    /// Get halfedge from vertex _vh1 to _vh2
    HalfEdgeHandle halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const;

    /// Get half-face from list of incident vertices (in connected order)
    ///
    /// \note Only the first three vertices are checked
    HalfFaceHandle halfface(const std::vector<VertexHandle>& _vs) const;

    /// Get half-face from list of incident vertices (in connected order)
    ///
    /// \note All vertices are checked
    HalfFaceHandle halfface_extensive(const std::vector<VertexHandle>& _vs) const;

    /// Get half-face from list of incident half-edges
    ///
    /// \note Only the first two half-edges are checked
    HalfFaceHandle halfface(const std::vector<HalfEdgeHandle>& _hes) const;

    /// Get next halfedge within a halfface
    HalfEdgeHandle next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get previous halfedge within a halfface
    HalfEdgeHandle prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const;

    /// Get the vertex the halfedge starts from
    VertexHandle from_vertex_handle(const HalfEdgeHandle& _h) const {
        return halfedge(_h).from_vertex();
    }

    /// Get the vertex the halfedge points to
    VertexHandle to_vertex_handle(const HalfEdgeHandle& _h) const {
        return halfedge(_h).to_vertex();
    }

    /// Get valence of vertex (number of incident edges)
    inline size_t valence(const VertexHandle& _vh) const {
        return VertexHalfEdgeIncidence::incident(_vh).size();
    }

    /// Get valence of edge (number of incident faces)
    inline size_t valence(const EdgeHandle& _eh) const {
        return HalfEdgeHalfFaceIncidence::incident(halfedge_handle(_eh, 0)).size();
    }

    /// Get valence of face (number of incident edges)
    inline size_t valence(const FaceHandle& _fh) const {
        assert(_fh.is_valid() && _fh.uidx() < faces_.size());

        return face(_fh).halfedges().size();
    }

    /// Get valence of cell (number of incident faces)
    inline size_t valence(const CellHandle& _ch) const {
        assert(_ch.is_valid() && _ch.uidx() < cells_.size());

        return cell(_ch).halffaces().size();
    }

    //=====================================================================
    // Delete entities
    //=====================================================================

public:

    // TODO: no need to make these virtual:
    virtual void delete_vertex(const VertexHandle& _h);
    virtual void delete_edge(const EdgeHandle& _h);
    virtual void delete_face(const FaceHandle& _h);
    virtual void delete_cell(const CellHandle& _h);
    virtual void collect_garbage();


    virtual bool is_deleted(const VertexHandle& _h)   const { return vertex_deleted_[_h.uidx()]; }
    virtual bool is_deleted(const EdgeHandle& _h)     const { return edge_deleted_[_h.uidx()];   }
    virtual bool is_deleted(const HalfEdgeHandle& _h) const { return edge_deleted_[_h.uidx()/2]; }
    virtual bool is_deleted(const FaceHandle& _h)     const { return face_deleted_[_h.uidx()];   }
    virtual bool is_deleted(const HalfFaceHandle& _h) const { return face_deleted_[_h.uidx()/2]; }
    virtual bool is_deleted(const CellHandle& _h)     const { return cell_deleted_[_h.uidx()];   }

private:

    template <class ContainerT>
    void get_incident_edges(const ContainerT& _vs, std::set<EdgeHandle>& _es) const;

    template <class ContainerT>
    void get_incident_faces(const ContainerT& _es, std::set<FaceHandle>& _fs) const;

    template <class ContainerT>
    void get_incident_cells(const ContainerT& _fs, std::set<CellHandle>& _cs) const;

    void delete_vertex_core(const VertexHandle& _h);
    void delete_edge_core(const EdgeHandle& _h);
    void delete_face_core(const FaceHandle& _h);
    void delete_cell_core(const CellHandle& _h);

public:

    /// Exchanges the indices of two cells while keeping the mesh otherwise unaffected.
    virtual void swap_cell_indices(CellHandle _h1, CellHandle _h2);

    /// Exchanges the indices of two faces while keeping the mesh otherwise unaffected.
    virtual void swap_face_indices(FaceHandle _h1, FaceHandle _h2);

    /// Exchanges the indices of two edges while keeping the mesh otherwise unaffected.
    virtual void swap_edge_indices(EdgeHandle _h1, EdgeHandle _h2);

    /// Exchanges the indices of two vertices while keeping the mesh otherwise unaffected.
    virtual void swap_vertex_indices(VertexHandle _h1, VertexHandle _h2);

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
        n_vertices_ = 0;

        if(_clearProps) {
            // Delete all property data
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

    void enable_vertex_bottom_up_incidences(bool _enable = true)
    {
        VertexHalfEdgeIncidence::set_enabled(_enable);
    }

    void enable_edge_bottom_up_incidences(bool _enable = true)
    {
        HalfEdgeHalfFaceIncidence::set_enabled(_enable);
    }

    void enable_face_bottom_up_incidences(bool _enable = true)
    {
        HalfFaceCellIncidence::set_enabled(_enable);
    }

    bool has_full_bottom_up_incidences() const {
        return (has_vertex_bottom_up_incidences() &&
                has_edge_bottom_up_incidences() &&
                has_face_bottom_up_incidences());
    }

    bool has_vertex_bottom_up_incidences() const { return VertexHalfEdgeIncidence::enabled(); }

    bool has_edge_bottom_up_incidences() const { return HalfEdgeHalfFaceIncidence::enabled(); }

    bool has_face_bottom_up_incidences() const { return HalfFaceCellIncidence::enabled(); }


    void enable_deferred_deletion(bool _enable = true);
    bool deferred_deletion_enabled() const { return deferred_deletion_; }


    [[deprecated("Fast deletion is always enabled now.")]]
    void enable_fast_deletion(bool = true) {}
    [[deprecated("Fast deletion is always enabled now.")]]
    bool fast_deletion_enabled() const { return true; }


private:

    bool deferred_deletion_ = true;

    //=====================================================================
    // Connectivity
    //=====================================================================

public:


    VertexHalfEdgeIncidence::Incidences const& outgoing_hes_per_vertex(VertexHandle vh) const;
    HalfEdgeHalfFaceIncidence::Incidences const& incident_hfs_per_he(HalfEdgeHandle heh) const;

    // TODO: this should probably be "common edge", the halfedge direction is ignored.
    //       but we probably also want a "common halfedge" API that may assume he \in hf
    /// \brief Get halfface that is adjacent (w.r.t. a common halfedge) within the same cell
    ///
    /// \param without_he_hf_incidences: for use when creating he-hf incidences
    ///
    /// \return Handle of the adjacent half-face if \a _halfFaceHandle is not
    ///         at a boundary, \a InvalidHalfFaceHandle otherwise.
    ///
    /// \warning The mesh must have face bottom-up incidences.
    HalfFaceHandle adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const;

    /// Get cell that is incident to the given halfface
    CellHandle incident_cell(const HalfFaceHandle& _halfFaceHandle) const;

    bool is_boundary(const HalfFaceHandle& _hfh) const {
        return !HalfFaceCellIncidence::incident(_hfh).is_valid();
    }

    bool is_boundary(const FaceHandle& _faceHandle) const {
        assert(_faceHandle.is_valid() && _faceHandle.uidx() < faces_.size());
        assert(has_face_bottom_up_incidences());
        return  is_boundary(halfface_handle(_faceHandle, 0)) ||
                is_boundary(halfface_handle(_faceHandle, 1));
    }

    bool is_boundary(const EdgeHandle& _edgeHandle) const {
        assert(has_edge_bottom_up_incidences());
        assert(_edgeHandle.is_valid() && _edgeHandle.uidx() < edges_.size());

        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(_edgeHandle, 0));
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const HalfEdgeHandle& _halfedgeHandle) const {
        assert(has_edge_bottom_up_incidences());
        assert(_halfedgeHandle.is_valid() && _halfedgeHandle.uidx() < edges_.size() * 2u);

        for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(_halfedgeHandle);
                hehf_it.valid(); ++hehf_it) {
            if(is_boundary(face_handle(*hehf_it))) {
                return true;
            }
        }
        return false;
    }

    bool is_boundary(const VertexHandle& _vertexHandle) const {
        assert(has_vertex_bottom_up_incidences());
        assert(_vertexHandle.is_valid() && _vertexHandle.uidx() < n_vertices());

        for(VertexOHalfEdgeIter voh_it = voh_iter(_vertexHandle); voh_it.valid(); ++voh_it) {
            if(is_boundary(*voh_it)) return true;
        }
        return false;
    }

    bool is_boundary(const CellHandle& _cellHandle) const {
        assert(_cellHandle.is_valid() && (size_t)_cellHandle.idx() < n_cells());

        for(CellFaceIter cf_it = cf_iter(_cellHandle); cf_it.valid(); ++cf_it) {
            if(is_boundary(*cf_it)) return true;
        }
        return false;
    }

    size_t n_vertices_in_cell(const CellHandle& _ch) const {
        assert(_ch.is_valid() && _ch.uidx() < cells_.size());

        std::set<VertexHandle> vhs;
        std::vector<HalfFaceHandle> hfs = cell(_ch).halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin();
                hf_it != hfs.end(); ++hf_it) {
            std::vector<HalfEdgeHandle> hes = halfface(*hf_it).halfedges();
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
                he_it != hes.end(); ++he_it) {
                vhs.insert(halfedge(*he_it).to_vertex());
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
    static inline HalfEdgeHandle halfedge_handle(const EdgeHandle& _h, const unsigned char _subIdx) {
        // Is handle in range?
        assert(_h.is_valid());
        assert(_subIdx < 2);
        // if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfEdgeHandle;
        return HalfEdgeHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Conversion function
    static inline HalfFaceHandle halfface_handle(const FaceHandle& _h, const unsigned char _subIdx) {
        // Is handle in range?
        assert(_h.is_valid());
        assert(_subIdx < 2);
        // if(_h.idx() < 0 || _subIdx > 1) return InvalidHalfFaceHandle;
        return HalfFaceHandle((2 * _h.idx()) + (_subIdx ? 1 : 0));
    }

    /// Handle conversion
    static inline EdgeHandle edge_handle(const HalfEdgeHandle& _h) {
        // Is handle in range?
        assert(_h.is_valid());
        // if(_h.idx() < 0) return InvalidEdgeHandle;
        return EdgeHandle((int)(_h.idx() / 2));
    }

    static inline FaceHandle face_handle(const HalfFaceHandle& _h) {
        // Is handle in range?
        assert(_h.is_valid());
        // if(_h.idx() < 0) return InvalidFaceHandle;
        return FaceHandle((int)(_h.idx() / 2));
    }

    static inline HalfEdgeHandle opposite_halfedge_handle(const HalfEdgeHandle& _h) {
        assert(_h.is_valid());
        return HalfEdgeHandle(_h.idx() ^ 1);
    }

    static inline HalfFaceHandle opposite_halfface_handle(const HalfFaceHandle& _h) {
        assert(_h.is_valid());
        return HalfFaceHandle(_h.idx() ^ 1);
    }

    bool inline needs_garbage_collection() const {
        return n_deleted_vertices_ > 0 || n_deleted_edges_ > 0 || n_deleted_faces_ > 0 || n_deleted_cells_ > 0;
    }

public:
    template<typename Entity>
    bool is_valid(HandleT<Entity> _h) const {
        return _h.uidx() < n<Entity>();
    }

protected:

    // List of edges
    std::vector<Edge> edges_;

    // List of faces
    std::vector<Face> faces_;

    // List of cells
    std::vector<Cell> cells_;

    std::vector<bool> vertex_deleted_;
    std::vector<bool> edge_deleted_;
    std::vector<bool> face_deleted_;
    std::vector<bool> cell_deleted_;

    // number of elements deleted, but not yet garbage collected
    size_t n_deleted_vertices_ = 0;
    size_t n_deleted_edges_ = 0;
    size_t n_deleted_faces_ = 0;
    size_t n_deleted_cells_ = 0;

};

}

