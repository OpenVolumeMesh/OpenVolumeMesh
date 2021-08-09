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

#include <algorithm>

#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

// Initialize constants
const VertexHandle      TopologyKernel::InvalidVertexHandle   = VertexHandle(-1);
const EdgeHandle        TopologyKernel::InvalidEdgeHandle     = EdgeHandle(-1);
const HalfEdgeHandle    TopologyKernel::InvalidHalfEdgeHandle = HalfEdgeHandle(-1);
const FaceHandle        TopologyKernel::InvalidFaceHandle     = FaceHandle(-1);
const HalfFaceHandle    TopologyKernel::InvalidHalfFaceHandle = HalfFaceHandle(-1);
const CellHandle        TopologyKernel::InvalidCellHandle     = CellHandle(-1);


void TopologyKernel::reserve_vertices(size_t n)
{
    ResourceManager::reserve_vprops(n);
    vertex_deleted_.reserve(n);
}

void TopologyKernel::reserve_edges(size_t n)
{
    ResourceManager::reserve_eprops(n);
    edges_.reserve(n);
    edge_deleted_.reserve(n);
}

void TopologyKernel::reserve_faces(size_t n)
{
    ResourceManager::reserve_fprops(n);
    faces_.reserve(n);
    face_deleted_.reserve(n);
}

void TopologyKernel::reserve_cells(size_t n)
{
    ResourceManager::reserve_cprops(n);
    cells_.reserve(n);
    cell_deleted_.reserve(n);
}

VertexHandle TopologyKernel::add_vertex() {

    ++n_vertices_;
    vertex_deleted_.push_back(false);

    // Create item for vertex bottom-up incidences
    if(has_vertex_bottom_up_incidences()) {
        outgoing_hes_per_vertex_.resize(n_vertices_);
    }

    // Resize vertex props
    resize_vprops(n_vertices_);

    // Return 0-indexed handle
    return VertexHandle((int)(n_vertices_ - 1));
}

//========================================================================================

/// Add edge
EdgeHandle TopologyKernel::add_edge(const VertexHandle& _fromVertex,
                                    const VertexHandle& _toVertex,
                                    bool _allowDuplicates) {

    // If the conditions are not fulfilled, assert will fail (instead
	// of returning an invalid handle)
    assert(_fromVertex.is_valid() && (size_t)_fromVertex.idx() < n_vertices() && !is_deleted(_fromVertex));
    assert(_toVertex.is_valid() && (size_t)_toVertex.idx() < n_vertices() && !is_deleted(_toVertex));

    // Test if edge does not exist, yet
    if(!_allowDuplicates) {
        auto heh = halfedge(_fromVertex, _toVertex);
        if (heh.is_valid()) {
            return edge_handle(heh);
        }
    }

    // Store edge
    edges_.emplace_back(_fromVertex, _toVertex);
    edge_deleted_.push_back(false);

    resize_eprops(n_edges());

    EdgeHandle eh((int)edges_.size()-1);

    VertexHalfEdgeIncidence::add_edge(eh, edges_.back());
    HalfEdgeHalfFaceIncidence::resize(n_halfedges());

    // Get handle of recently created edge
    return eh;
}

//========================================================================================

/// Add face via incident edges
FaceHandle TopologyKernel::add_face(std::vector<HalfEdgeHandle> _halfedges, bool _topologyCheck) {

#ifndef NDEBUG
    // Assert that halfedges are valid
    for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin(),
            end = _halfedges.end(); it != end; ++it) {
        is_valid(*it);
    }
#endif

    // Perform topology check
    if(_topologyCheck) {
        for (size_t i = 0; i + 1< _halfedges.size(); ++i) {
            if (to_vertex_handle(_halfedges[i]) != from_vertex_handle(_halfedges[i+1])) {
                return InvalidFaceHandle;
            }
        }
        if (to_vertex_handle(_halfedges.back()) != from_vertex_handle(_halfedges.front())) {
            return InvalidFaceHandle;
        }
        // The halfedges are now guaranteed to be connected
    }

    // Create face
    faces_.emplace_back(std::move(_halfedges));
    face_deleted_.push_back(false);

    // Get added face's handle
    FaceHandle fh((int)faces_.size() - 1);

    // Resize props
    resize_fprops(n_faces());

    HalfEdgeHalfFaceIncidence::add_face(fh, faces_.back());
    HalfFaceCellIncidence::resize(n_halffaces());

    return fh;
}

//========================================================================================

// TODO: add iterator-based interface + range adapter

/// Add face via incident vertices
/// Define the _vertices in counter-clockwise order (from the "outside")
FaceHandle TopologyKernel::add_face(const std::vector<VertexHandle>& _vertices) {

#ifndef NDEBUG
    // Assert that all vertices have valid indices
    for(std::vector<VertexHandle>::const_iterator it = _vertices.begin(),
            end = _vertices.end(); it != end; ++it)
        assert(it->is_valid() && (size_t)it->idx() < n_vertices() && !is_deleted(*it));
#endif

    // Add edge for each pair of vertices
    std::vector<HalfEdgeHandle> halfedges;
    std::vector<VertexHandle>::const_iterator it = _vertices.begin();
    std::vector<VertexHandle>::const_iterator end = _vertices.end();
    for(; (it+1) != end; ++it) {
        EdgeHandle e_idx = add_edge(*it, *(it+1));

        // Swap halfedge if edge already existed and
        // has been initially defined in reverse orientation
        char swap = edge(e_idx).to_vertex() == *it;

        halfedges.push_back(halfedge_handle(e_idx, swap));
    }
    EdgeHandle e_idx = add_edge(*it, *_vertices.begin());
    char swap = edge(e_idx).to_vertex() == *it;
    halfedges.push_back(halfedge_handle(e_idx, swap));

    // Add face
#ifndef NDEBUG
    return add_face(halfedges, true);
#else
    return add_face(halfedges, false);
#endif
}

//========================================================================================

void TopologyKernel::reorder_incident_halffaces(const EdgeHandle& _eh) {

}

//========================================================================================

/// Add cell via incident halffaces
CellHandle TopologyKernel::add_cell(std::vector<HalfFaceHandle> _halffaces, bool _topologyCheck)
{
    assert(_halffaces.size() > 0);

#ifndef NDEBUG
    // Assert that halffaces have valid indices
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin(),
            end = _halffaces.end(); it != end; ++it)
        assert(it->is_valid() && ((size_t)it->idx() < faces_.size() * 2u) && !is_deleted(*it));
#endif

    // Perform topology check
    if(_topologyCheck) {

        /*
         * We test the following necessary properties for a closed 2-manifold cell:
         *   - each halfedge may only be used once (this implies a halfface may only be used once)
         *   - if a halfedge is used, its opposite halfface must be used too
         */

        // collect a vector of all used halfedges
        std::vector<HalfEdgeHandle> incidentHalfedges;
        size_t guess_n_halfedges = _halffaces.size() * valence(face_handle(_halffaces[0]));
        // actually, use double the guess, better to allocate a bit more than
        // risk reallocation.
        incidentHalfedges.reserve(2 * guess_n_halfedges);

        for (const auto &hfh: _halffaces) {
            const auto &hes = face(face_handle(hfh)).halfedges();
            if ((hfh.idx() & 1) == 0) { // first halfface
                std::copy(hes.begin(), hes.end(),
                        std::back_inserter(incidentHalfedges));
            } else {
                std::transform(hes.rbegin(),
                        hes.rend(),
                        std::back_inserter(incidentHalfedges),
                        opposite_halfedge_handle);
            }
        }
        std::sort(incidentHalfedges.begin(), incidentHalfedges.end());
        auto duplicate = std::adjacent_find(incidentHalfedges.begin(), incidentHalfedges.end());
        if (duplicate != incidentHalfedges.end()) {
#ifndef NDEBUG
            std::cerr << "add_cell(): Halfedge #" << duplicate->idx() << " is contained in more than 1 halfface." << std::endl;
#endif
            return InvalidCellHandle;
        }
        size_t n_halfedges = incidentHalfedges.size();
        auto e_end = std::unique(incidentHalfedges.begin(), incidentHalfedges.end(),
                [](HalfEdgeHandle a, HalfEdgeHandle b) {return a.idx()/2 == b.idx()/2;});
        auto n_edges = static_cast<size_t>(std::distance(incidentHalfedges.begin(), e_end));

        if(n_halfedges != 2u * n_edges) {
#ifndef NDEBUG
            std::cerr << "add_cell(): The specified half-faces are not connected!" << std::endl;
#endif
            return InvalidCellHandle;
        }

        // The halffaces are now guaranteed to form a two-manifold
    }

    // Create new cell
    cells_.emplace_back(std::move(_halffaces));
    cell_deleted_.push_back(false);

    // Resize props
    resize_cprops(n_cells());

    CellHandle ch((int)cells_.size()-1);

    // Update face bottom-up incidences
    if(has_face_bottom_up_incidences()) {
        HalfFaceCellIncidence::add_cell(ch, cells_.back());
    }

    if (has_edge_bottom_up_incidences() && has_face_bottom_up_incidences()) {
        for (const auto &hfh: cells_.back().halffaces()) {
            // both of each incident edge's halfedges are included exactly once,
            // so we just pick the ones with subidx 0.
            for (const auto &heh: halfface_halfedges(hfh)) {
                if ((heh.idx() & 1) == 0) {
                    HalfEdgeHalfFaceIncidence::reorder_halffaces(edge_handle(heh));
                }
            }
        }
    }

    return ch;
}

//========================================================================================

/// Set the vertices of an edge
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_edge(const EdgeHandle& _eh, const VertexHandle& _fromVertex, const VertexHandle& _toVertex) {

    assert(_fromVertex.is_valid() && (size_t)_fromVertex.idx() < n_vertices() && !is_deleted(_fromVertex));
    assert(_toVertex.is_valid() && (size_t)_toVertex.idx() < n_vertices() && !is_deleted(_toVertex));

    Edge& e = edge(_eh);

    VertexHalfEdgeIncidence::delete_edge(_eh, e);
    e.set_from_vertex(_fromVertex);
    e.set_to_vertex(_toVertex);
    VertexHalfEdgeIncidence::add_edge(_eh, e);
}

//========================================================================================

/// Set the half-edges of a face
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_face(const FaceHandle& _fh, const std::vector<HalfEdgeHandle>& _hes) {

    Face& f = face(_fh);

    HalfEdgeHalfFaceIncidence::delete_face(_fh, f);
    f.set_halfedges(_hes);
    HalfEdgeHalfFaceIncidence::add_face(_fh, f);
    HalfEdgeHalfFaceIncidence::invalidate_order(_fh);
}

//========================================================================================

/// Set the half-faces of a cell
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_cell(const CellHandle& _ch, const std::vector<HalfFaceHandle>& _hfs) {

    Cell& c = cell(_ch);

    HalfFaceCellIncidence::delete_cell(_ch, c);
    c.set_halffaces(_hfs);
    HalfFaceCellIncidence::add_cell(_ch, c);

    HalfEdgeHalfFaceIncidence::invalidate_order(_ch);
}

//========================================================================================

/**
 * \brief Delete vertex from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to vertex _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the vertex to be deleted
 */
void TopologyKernel::delete_vertex(const VertexHandle& _h) {

    assert(!is_deleted(_h));

    std::vector<VertexHandle> vs;
    vs.push_back(_h);

    std::set<EdgeHandle> incidentEdges_s;
    get_incident_edges(vs, incidentEdges_s);

    std::set<FaceHandle> incidentFaces_s;
    get_incident_faces(incidentEdges_s, incidentFaces_s);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(incidentFaces_s, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete faces
    for(std::set<FaceHandle>::const_reverse_iterator f_it = incidentFaces_s.rbegin(),
            f_end = incidentFaces_s.rend(); f_it != f_end; ++f_it) {
        delete_face_core(*f_it);
    }

    // Delete edges
    for(std::set<EdgeHandle>::const_reverse_iterator e_it = incidentEdges_s.rbegin(),
            e_end = incidentEdges_s.rend(); e_it != e_end; ++e_it) {
        delete_edge_core(*e_it);
    }

    // Delete vertex
    delete_vertex_core(_h);
}

//========================================================================================

/**
 * \brief Delete edge from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to edge _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the edge to be deleted
 */
void TopologyKernel::delete_edge(const EdgeHandle& _h) {

    assert(!is_deleted(_h));

    std::vector<EdgeHandle> es;
    es.push_back(_h);

    std::set<FaceHandle> incidentFaces_s;
    get_incident_faces(es, incidentFaces_s);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(incidentFaces_s, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete faces
    for(std::set<FaceHandle>::const_reverse_iterator f_it = incidentFaces_s.rbegin(),
            f_end = incidentFaces_s.rend(); f_it != f_end; ++f_it) {
        delete_face_core(*f_it);
    }

    // Delete edge
    delete_edge_core(_h);
}

//========================================================================================

/**
 * \brief Delete face from mesh
 *
 * Get all incident higher-dimensional entities and delete the complete
 * subtree of the mesh incident to face _h.
 * In this function all incident entities are gathered
 * and deleted using the delete_*_core functions
 * that do the actual deletion including the update
 * of the bottom-up incidences, etc.
 *
 * @param _h The handle to the face to be deleted
 */
void TopologyKernel::delete_face(const FaceHandle& _h) {

    assert(!is_deleted(_h));

    std::vector<FaceHandle> fs;
    fs.push_back(_h);

    std::set<CellHandle> incidentCells_s;
    get_incident_cells(fs, incidentCells_s);

    // Delete cells
    for(std::set<CellHandle>::const_reverse_iterator c_it = incidentCells_s.rbegin(),
            c_end = incidentCells_s.rend(); c_it != c_end; ++c_it) {
        delete_cell_core(*c_it);
    }

    // Delete face
    delete_face_core(_h);
}

//========================================================================================

/**
 * \brief Delete cell from mesh
 *
 * Since there's no higher dimensional incident
 * entity to a cell, we can safely delete it from the
 * mesh.
 *
 * @param _h The handle to the cell to be deleted
 */
void TopologyKernel::delete_cell(const CellHandle& _h) {

    assert(!is_deleted(_h));
    delete_cell_core(_h);
}

/**
 * \brief Delete all entities that are marked as deleted
 */
void TopologyKernel::collect_garbage()
{
    if (!deferred_deletion_enabled() || !needs_garbage_collection())
        return; // nothing todo

    deferred_deletion_ = false;

    for (int i = (int)n_cells(); i > 0; --i) {
        if (is_deleted(CellHandle(i - 1))) {
            cell_deleted_[i - 1] = false;
            delete_cell_core(CellHandle(i - 1));
        }
    }
    n_deleted_cells_ = 0;

    for (int i = (int)n_faces(); i > 0; --i) {
        if (is_deleted(FaceHandle(i - 1))) {
            face_deleted_[i - 1] = false;
            delete_face_core(FaceHandle(i - 1));
        }
    }
    n_deleted_faces_ = 0;

    for (int i = (int)n_edges(); i > 0; --i) {
        if (is_deleted(EdgeHandle(i - 1))) {
            edge_deleted_[i - 1] = false;
            delete_edge_core(EdgeHandle(i - 1));
        }
    }
    n_deleted_edges_ = 0;

    for (int i = (int)n_vertices(); i > 0; --i) {
        if (is_deleted(VertexHandle(i - 1))) {
            vertex_deleted_[i - 1] = false;
            delete_vertex_core(VertexHandle(i - 1));
        }
    }
    n_deleted_vertices_ = 0;

    deferred_deletion_ = true;

}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_edges(const ContainerT& _vs,
                                        std::set<EdgeHandle>& _es) const {

    _es.clear();

    if(has_vertex_bottom_up_incidences()) {

        for(typename ContainerT::const_iterator v_it = _vs.begin(),
                v_end = _vs.end(); v_it != v_end; ++v_it) {

            auto const &inc_hes = VertexHalfEdgeIncidence::incident(*v_it);

            for(std::vector<HalfEdgeHandle>::const_iterator he_it = inc_hes.begin(),
                    he_end = inc_hes.end(); he_it != he_end; ++he_it) {

                _es.insert(edge_handle(*he_it));
            }
        }
    } else {

        for(typename ContainerT::const_iterator v_it = _vs.begin(),
                v_end = _vs.end(); v_it != v_end; ++v_it) {

            for(EdgeIter e_it = edges_begin(), e_end = edges_end(); e_it != e_end; ++e_it) {

                const Edge& e = edge(*e_it);

                if(e.from_vertex() == *v_it || e.to_vertex() == *v_it) {
                    _es.insert(*e_it);
                }
            }
        }
    }
}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_faces(const ContainerT& _es,
                                        std::set<FaceHandle>& _fs) const {

    _fs.clear();

    if(has_edge_bottom_up_incidences()) {

        for(typename ContainerT::const_iterator e_it = _es.begin(),
                e_end = _es.end(); e_it != e_end; ++e_it) {

            for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(halfedge_handle(*e_it, 0));
                    hehf_it.valid(); ++hehf_it) {

                const FaceHandle fh = face_handle(*hehf_it);

                _fs.insert(fh);
            }
        }
    } else {

        for(typename ContainerT::const_iterator e_it = _es.begin(),
                e_end = _es.end(); e_it != e_end; ++e_it) {

            for(FaceIter f_it = faces_begin(),
                    f_end = faces_end(); f_it != f_end; ++f_it) {

                const std::vector<HalfEdgeHandle>& hes = face(*f_it).halfedges();

                for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                        he_end = hes.end(); he_it != he_end; ++he_it) {

                    if(edge_handle(*he_it) == *e_it) {
                        _fs.insert(*f_it);
                        break;
                    }
                }
            }
        }
    }
}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_cells(const ContainerT& _fs,
                                        std::set<CellHandle>& _cs) const {

    _cs.clear();

    if(has_face_bottom_up_incidences()) {

        for(typename ContainerT::const_iterator f_it = _fs.begin(),
            f_end = _fs.end(); f_it != f_end; ++f_it) {

            const HalfFaceHandle hfh0 = halfface_handle(*f_it, 0);
            const HalfFaceHandle hfh1 = halfface_handle(*f_it, 1);

            const CellHandle c0 = incident_cell(hfh0);
            const CellHandle c1 = incident_cell(hfh1);

            if(c0.is_valid()) _cs.insert(c0);
            if(c1.is_valid()) _cs.insert(c1);
        }
    } else {

        for(typename ContainerT::const_iterator f_it = _fs.begin(),
            f_end = _fs.end(); f_it != f_end; ++f_it) {

            for(CellIter c_it = cells_begin(), c_end = cells_end();
                c_it != c_end; ++c_it) {

                const std::vector<HalfFaceHandle>& hfs = cell(*c_it).halffaces();

                for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                        hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

                    if(face_handle(*hf_it) == *f_it) {
                        _cs.insert(*c_it);
                        break;
                    }
                }
            }
        }
    }
}

//========================================================================================

/**
 * \brief Delete vertex from mesh
 *
 * There may be no more edges, faces, or cells incident to this edge.
 *
 * @param _h A vertex's handle
 */
void TopologyKernel::delete_vertex_core(const VertexHandle& _h) {

    VertexHandle h = _h;
    assert(h.is_valid() && (size_t)h.idx() < n_vertices());

    if (deferred_deletion_enabled())
    {
        ++n_deleted_vertices_;
        vertex_deleted_[h.idx()] = true;
        return;
    }

    VertexHandle last_undeleted_vertex = VertexHandle((int)n_vertices()-1);
    assert(!vertex_deleted_[last_undeleted_vertex.idx()]);
    swap_vertex_indices(h, last_undeleted_vertex);
    h = last_undeleted_vertex;

    if(has_vertex_bottom_up_incidences()) {
        assert((size_t)h.idx() < outgoing_hes_per_vertex_.size());
        outgoing_hes_per_vertex_.erase(outgoing_hes_per_vertex_.begin() + h.idx());
    }

    --n_vertices_;
    vertex_deleted_.erase(vertex_deleted_.begin() + h.idx());

    ResourceManager::vertex_deleted(h);

    return;
}

//========================================================================================

/**
 * \brief Delete edge from mesh
 *
 * There may be no more faces or cells incident to this edge.
 *
 * @param _h An edge's handle
 */
void TopologyKernel::delete_edge_core(const EdgeHandle& _h) {

    EdgeHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < edges_.size());

    if (!deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        EdgeHandle last_edge = EdgeHandle((int)edges_.size()-1);
        assert(!edge_deleted_[last_edge.idx()]);
        swap_edge_indices(h, last_edge);
        h = last_edge;
    }

    VertexHalfEdgeIncidence::delete_edge(h, edge(h));

    if (deferred_deletion_enabled())
    {
        ++n_deleted_edges_;
        edge_deleted_[h.idx()] = true;
        return;
    }

    if(has_edge_bottom_up_incidences()) {
        assert((size_t)halfedge_handle(h, 1).idx() < incident_hfs_per_he_.size());

        incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 1).idx());
        incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 0).idx());
    }

    // 5)
    edges_.erase(edges_.begin() + h.idx());
    edge_deleted_.erase(edge_deleted_.begin() + h.idx());


    // 6)

    ResourceManager::edge_deleted(h);
}

//========================================================================================

/**
 * \brief Delete face from mesh
 *
 * There may not be any cells incident to this face.
 *
 * @param _h An face's handle
 */
void TopologyKernel::delete_face_core(const FaceHandle& _h) {

    FaceHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < faces_.size());


    if (!deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        FaceHandle last_face = FaceHandle((int)faces_.size()-1);
        assert(!face_deleted_[last_face.idx()]);
        swap_face_indices(h, last_face);
        h = last_face;
    }

    // 1)
    if(has_edge_bottom_up_incidences()) {

        const std::vector<HalfEdgeHandle>& hes = face(h).halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

            assert((size_t)std::max(he_it->idx(), opposite_halfedge_handle(*he_it).idx()) < incident_hfs_per_he_.size());

            incident_hfs_per_he_[he_it->idx()].erase(
                    std::remove(incident_hfs_per_he_[he_it->idx()].begin(),
                                incident_hfs_per_he_[he_it->idx()].end(),
                                halfface_handle(h, 0)), incident_hfs_per_he_[he_it->idx()].end());


            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].erase(
                    std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].begin(),
                                incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end(),
                                halfface_handle(h, 1)), incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end());

            reorder_incident_halffaces(edge_handle(*he_it));
        }
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_faces_;
        face_deleted_[h.idx()] = true;
        return;
    }
    // 3)
    if(has_face_bottom_up_incidences()) {
        assert((size_t)halfface_handle(h, 1).idx() < incident_cell_per_hf_.size());

        incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 1).idx());
        incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 0).idx());
    }

    faces_.erase(faces_.begin() + h.idx());
    face_deleted_.erase(face_deleted_.begin() + h.idx());

    ResourceManager::face_deleted(h);
}

//========================================================================================

/**
 * \brief Delete cell from mesh
 *
 * After performing this operation, all cells
 * following cell _h in the array will be accessible
 * through their old handle decreased by one.
 * These steps are performed:
 *
 * 1) Delete links in BU: HF -> C
 * 2) Decrease all entries > c in BU: HF -> C
 * 3) Delete cell from storage array
 * 4) Delete property item
 *
 * @param _h A cell handle
 */
void TopologyKernel::delete_cell_core(const CellHandle& _h) {

    CellHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < cells_.size());


    if (!deferred_deletion_enabled()) // for fast deletion swap handle with last not deleted cell
    {
        CellHandle last_undeleted_cell = CellHandle((int)cells_.size()-1);
        assert(!cell_deleted_[last_undeleted_cell.idx()]);
        swap_cell_indices(h, last_undeleted_cell);
        h = last_undeleted_cell;
    }


    // 1)
    if(has_face_bottom_up_incidences()) {
        const std::vector<HalfFaceHandle>& hfs = cell(h).halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {
            assert((size_t)hf_it->idx() < incident_cell_per_hf_.size());
            if (incident_cell_per_hf_[hf_it->idx()] == h)
                incident_cell_per_hf_[hf_it->idx()] = InvalidCellHandle;
        }
        std::set<EdgeHandle> edges;
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {
          const auto& hf = halfface(*hf_it);
          for (const auto&  heh : hf.halfedges())
            edges.insert(edge_handle(heh));
        }
        for (auto eh : edges)
          reorder_incident_halffaces(eh);
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_cells_;
        cell_deleted_[h.idx()] = true;
        return;
    }

    // 3)
    cells_.erase(cells_.begin() + h.idx());
    cell_deleted_.erase(cell_deleted_.begin() + h.idx());

    // 4)
    ResourceManager::cell_deleted(h);
}

void TopologyKernel::swap_cell_indices(CellHandle _h1, CellHandle _h2)
{
    assert(_h1.idx() >= 0 && _h1.idx() < (int)cells_.size());
    assert(_h2.idx() >= 0 && _h2.idx() < (int)cells_.size());

    if (_h1 == _h2)
        return;

    int id1 = _h1.idx();
    int id2 = _h2.idx();

    Cell c1 = cells_[id1];
    Cell c2 = cells_[id2];

    // correct pointers to those cells
    std::vector<HalfFaceHandle> hfhs1 = c1.halffaces();
    for (unsigned int i = 0; i < hfhs1.size(); ++i)
    {
        HalfFaceHandle hfh = hfhs1[i];
        if (incident_cell_per_hf_[hfh.idx()] == _h1)
            incident_cell_per_hf_[hfh.idx()] = _h2;
    }

    std::vector<HalfFaceHandle> hfhs2 = c2.halffaces();
    for (unsigned int i = 0; i < hfhs2.size(); ++i)
    {
        HalfFaceHandle hfh = hfhs2[i];
        if (incident_cell_per_hf_[hfh.idx()] == _h2)
            incident_cell_per_hf_[hfh.idx()] = _h1;
    }

    // swap vector entries
    std::swap(cells_[id1], cells_[id2]);
    bool tmp = cell_deleted_[id1];
    cell_deleted_[id1] = cell_deleted_[id2];
    cell_deleted_[id2] = tmp;
    swap_cell_properties(_h1, _h2);
}

void TopologyKernel::swap_face_indices(FaceHandle _h1, FaceHandle _h2)
{
    assert(_h1.idx() >= 0 && _h1.idx() < (int)faces_.size());
    assert(_h2.idx() >= 0 && _h2.idx() < (int)faces_.size());

    if (_h1 == _h2)
        return;


    std::vector<unsigned int> ids;
    ids.push_back(_h1.idx());
    ids.push_back(_h2.idx());

    unsigned int id1 = _h1.idx();
    unsigned int id2 = _h2.idx();

    // correct pointers to those faces

    // correct cells that contain a swapped faces
    if (has_face_bottom_up_incidences())
    {
        std::set<unsigned int> processed_cells; // to ensure ids are only swapped once (in the case that the two swapped faces belong to a common cell)
        for (unsigned int i = 0; i < 2; ++i) // For both swapped faces
        {
            unsigned int id = ids[i];
            for (unsigned int j = 0; j < 2; ++j) // for both halffaces
            {
                HalfFaceHandle hfh = HalfFaceHandle(2*id+j);
                CellHandle ch = incident_cell_per_hf_[hfh.idx()];
                if (!ch.is_valid())
                    continue;



                if (processed_cells.find(ch.idx()) == processed_cells.end())
                {

                    Cell& c = cells_[ch.idx()];

                    // replace old halffaces with new halffaces where the ids are swapped

                    std::vector<HalfFaceHandle> new_halffaces;
                    for (unsigned int k = 0; k < c.halffaces().size(); ++k)
                        if (c.halffaces()[k].idx()/2 == (int)id1) // if halfface belongs to swapped face
                            new_halffaces.push_back(HalfFaceHandle(2 * id2 + (c.halffaces()[k].idx() % 2)));
                        else if (c.halffaces()[k].idx()/2 == (int)id2) // if halfface belongs to swapped face
                            new_halffaces.push_back(HalfFaceHandle(2 * id1 + (c.halffaces()[k].idx() % 2)));
                        else
                            new_halffaces.push_back(c.halffaces()[k]);
                    c.set_halffaces(new_halffaces);

                    processed_cells.insert(ch.idx());
                }
            }
        }
    }
    else
    {
        // serach for all cells that contain a swapped face
        for (unsigned int i = 0; i < cells_.size(); ++i)
        {
            Cell& c = cells_[i];

            // check if c contains a swapped face
            bool contains_swapped_face = false;
            for (unsigned int k = 0; k < c.halffaces().size(); ++k)
            {
                if (c.halffaces()[k].idx()/2 == (int)id1)
                    contains_swapped_face = true;
                if (c.halffaces()[k].idx()/2 == (int)id2)
                    contains_swapped_face = true;
                if (contains_swapped_face)
                    break;
            }

            if (contains_swapped_face)
            {
            // replace old halffaces with new halffaces where the ids are swapped
                std::vector<HalfFaceHandle> new_halffaces;
                for (unsigned int k = 0; k < c.halffaces().size(); ++k)
                    if (c.halffaces()[k].idx()/2 == (int)id1) // if halfface belongs to swapped face
                        new_halffaces.push_back(HalfFaceHandle(2 * id2 + (c.halffaces()[k].idx() % 2)));
                    else if (c.halffaces()[k].idx()/2 == (int)id2) // if halfface belongs to swapped face
                        new_halffaces.push_back(HalfFaceHandle(2 * id1 + (c.halffaces()[k].idx() % 2)));
                    else
                        new_halffaces.push_back(c.halffaces()[k]);
                c.set_halffaces(new_halffaces);
            }
        }
    }

    // correct bottom up indices

    if (has_edge_bottom_up_incidences())
    {
        std::set<HalfEdgeHandle> processed_halfedges; // to ensure ids are only swapped once (in the case that a halfedge is incident to both swapped faces)
        for (unsigned int i = 0; i < 2; ++i) // For both swapped faces
        {
            unsigned int id = ids[i];
            for (unsigned int j = 0; j < 2; ++j) // for both halffaces
            {
                HalfFaceHandle hfh = HalfFaceHandle(2*id+j);
                Face hf = halfface(hfh);

                for (unsigned int k = 0; k < hf.halfedges().size(); ++k)
                {
                    HalfEdgeHandle heh = hf.halfedges()[k];

                    if (processed_halfedges.find(heh) != processed_halfedges.end())
                        continue;

                    std::vector<HalfFaceHandle>& incident_halffaces = incident_hfs_per_he_[heh.idx()];
                    for (unsigned int l = 0; l < incident_halffaces.size(); ++l)
                    {
                        HalfFaceHandle& hfh2 = incident_halffaces[l];

                        if (hfh2.idx()/2 == (int)id1) // if halfface belongs to swapped face
                            hfh2 = HalfFaceHandle(2 * id2 + (hfh2.idx() % 2));
                        else if (hfh2.idx()/2 == (int)id2) // if halfface belongs to swapped face
                            hfh2 = HalfFaceHandle(2 * id1 + (hfh2.idx() % 2));
                    }

                    processed_halfedges.insert(heh);
                }
            }
        }
    }

    // swap vector entries
    std::swap(faces_[ids[0]], faces_[ids[1]]);
    bool tmp = face_deleted_[ids[0]];
    face_deleted_[ids[0]] = face_deleted_[ids[1]];
    face_deleted_[ids[1]] = tmp;
    std::swap(incident_cell_per_hf_[2*ids[0]+0], incident_cell_per_hf_[2*ids[1]+0]);
    std::swap(incident_cell_per_hf_[2*ids[0]+1], incident_cell_per_hf_[2*ids[1]+1]);
    swap_face_properties(_h1, _h2);
    swap_halfface_properties(halfface_handle(_h1, 0), halfface_handle(_h2, 0));
    swap_halfface_properties(halfface_handle(_h1, 1), halfface_handle(_h2, 1));

}

void TopologyKernel::swap_edge_indices(EdgeHandle _h1, EdgeHandle _h2)
{
    assert(_h1.idx() >= 0 && _h1.idx() < (int)edges_.size());
    assert(_h2.idx() >= 0 && _h2.idx() < (int)edges_.size());

    if (_h1 == _h2)
        return;

    std::vector<unsigned int> ids;
    ids.push_back(_h1.idx());
    ids.push_back(_h2.idx());


    // correct pointers to those edges

    if (has_edge_bottom_up_incidences())
    {
        std::set<unsigned int> processed_faces; // to ensure ids are only swapped once (in the case that the two swapped edges belong to a common face)

        for (unsigned int i = 0; i < 2; ++i) // For both swapped edges
        {
            HalfEdgeHandle heh = HalfEdgeHandle(2*ids[i]);


            std::vector<HalfFaceHandle>& incident_halffaces = incident_hfs_per_he_[heh.idx()];
            for (unsigned int j = 0; j < incident_halffaces.size(); ++j) // for each incident halfface
            {
                HalfFaceHandle hfh = incident_halffaces[j];

                unsigned int f_id = hfh.idx() / 2;

                if (processed_faces.find(f_id) == processed_faces.end())
                {

                    Face& f = faces_[f_id];

                    // replace old incident halfedges with new incident halfedges where the ids are swapped
                    std::vector<HalfEdgeHandle> new_halfedges;
                    for (unsigned int k = 0; k < f.halfedges().size(); ++k)
                    {
                        HalfEdgeHandle heh2 = f.halfedges()[k];
                        if (heh2.idx() / 2 == (int)ids[0])
                            new_halfedges.push_back(HalfEdgeHandle(2*ids[1] + (heh2.idx() % 2)));
                        else if (heh2.idx() / 2 == (int)ids[1])
                            new_halfedges.push_back(HalfEdgeHandle(2*ids[0] + (heh2.idx() % 2)));
                        else
                            new_halfedges.push_back(heh2);
                    }
                    f.set_halfedges(new_halfedges);

                    processed_faces.insert(f_id);
                }
            }
        }
    }
    else
    {
        // search for all faces that contain one of the swapped edges
        for (unsigned int i = 0; i < faces_.size(); ++i)
        {
            Face& f = faces_[i];

            // check if f contains a swapped edge
            bool contains_swapped_edge = false;
            for (unsigned int k = 0; k < f.halfedges().size(); ++k)
            {
                if (f.halfedges()[k].idx()/2 == (int)ids[0])
                    contains_swapped_edge = true;
                if (f.halfedges()[k].idx()/2 == (int)ids[1])
                    contains_swapped_edge = true;
                if (contains_swapped_edge)
                    break;
            }

            if (contains_swapped_edge)
            {
                // replace old incident halfedges with new incident halfedges where the ids are swapped
                std::vector<HalfEdgeHandle> new_halfedges;
                for (unsigned int k = 0; k < f.halfedges().size(); ++k)
                {
                    HalfEdgeHandle heh2 = f.halfedges()[k];
                    if (heh2.idx() / 2 == (int)ids[0])
                        new_halfedges.push_back(HalfEdgeHandle(2*ids[1] + (heh2.idx() % 2)));
                    else if (heh2.idx() / 2 == (int)ids[1])
                        new_halfedges.push_back(HalfEdgeHandle(2*ids[0] + (heh2.idx() % 2)));
                    else
                        new_halfedges.push_back(heh2);
                }
                f.set_halfedges(new_halfedges);
            }
        }
    }

    // correct bottom up incidences

    if (has_vertex_bottom_up_incidences())
    {
        std::set<VertexHandle> processed_vertices;
        for (unsigned int i = 0; i < 2; ++i) // For both swapped edges
        {
            Edge e = edge(EdgeHandle(ids[i]));
            std::vector<VertexHandle> vhs;
            vhs.push_back(e.from_vertex());
            vhs.push_back(e.to_vertex());

            for (unsigned int j = 0; j < 2; ++j) // for both incident vertices
            {
                if (processed_vertices.find(vhs[j]) != processed_vertices.end())
                    continue;

                std::vector<HalfEdgeHandle>& outgoing_hes = outgoing_hes_per_vertex_[vhs[j].idx()];
                for (unsigned int k = 0; k < outgoing_hes.size(); ++k)
                {
                    HalfEdgeHandle& heh = outgoing_hes[k];
                    if (heh.idx() / 2 == (int)ids[0])
                        heh = HalfEdgeHandle(2 * ids[1] + (heh.idx() % 2));
                    else if (heh.idx() / 2 == (int)ids[1])
                        heh = HalfEdgeHandle(2 * ids[0] + (heh.idx() % 2));
                }
                processed_vertices.insert(vhs[j]);
            }

        }
    }

    // swap vector entries
    std::swap(edges_[ids[0]], edges_[ids[1]]);
    bool tmp = edge_deleted_[ids[0]];
    edge_deleted_[ids[0]] = edge_deleted_[ids[1]];
    edge_deleted_[ids[1]] = tmp;
    std::swap(incident_hfs_per_he_[2*ids[0]+0], incident_hfs_per_he_[2*ids[1]+0]);
    std::swap(incident_hfs_per_he_[2*ids[0]+1], incident_hfs_per_he_[2*ids[1]+1]);
    swap_edge_properties(_h1, _h2);
    swap_halfedge_properties(halfedge_handle(_h1, 0), halfedge_handle(_h2, 0));
    swap_halfedge_properties(halfedge_handle(_h1, 1), halfedge_handle(_h2, 1));
}

void TopologyKernel::swap_vertex_indices(VertexHandle _h1, VertexHandle _h2)
{
    assert(_h1.idx() >= 0 && _h1.idx() < (int)n_vertices_);
    assert(_h2.idx() >= 0 && _h2.idx() < (int)n_vertices_);

    if (_h1 == _h2)
        return;

    auto swap_edge_indices = [=](Edge &e) {
        if (e.from_vertex() == _h1)
            e.set_from_vertex(_h2);
        else if (e.from_vertex() == _h2)
            e.set_from_vertex(_h1);
        if (e.to_vertex() == _h1)
            e.set_to_vertex(_h2);
        else if (e.to_vertex() == _h2)
            e.set_to_vertex(_h1);
    };

    if (has_vertex_bottom_up_incidences())
    {
        std::set<EdgeHandle> processed_edges; // to ensure ids are only swapped once (in the case that the two swapped vertices are connected by an edge)
        for (const auto vh: {_h1, _h2})
        {
            for (const auto &heh: outgoing_halfedges(vh))
            {
                EdgeHandle eh = edge_handle(heh);
                if (processed_edges.find(eh) == processed_edges.end()) {
                    continue;
                }
                swap_edge_indices(edge(eh));
                processed_edges.insert(eh);
            }
        }

    }
    else
    {
        // search for all edges containing a swapped vertex
        for (const auto &eh: edges()) {
            swap_edge_indices(edge(eh));
        }
    }

    // swap vector entries
    std::swap(vertex_deleted_[_h1.uidx()], vertex_deleted_[_h2.uidx()]);
    std::swap(outgoing_hes_per_vertex_[_h1.uidx()], outgoing_hes_per_vertex_[_h2.uidx()]);
    swap_vertex_properties(_h1, _h2);
}

//========================================================================================

void TopologyKernel::delete_multiple_vertices(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_vertices());

    std::vector<int> newIndices(n_vertices(), -1);
    int curIdx = 0;

    std::vector<int>::iterator idx_it = newIndices.begin();
    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it) {

        if(!(*t_it)) {
            // Not marked as deleted
            *idx_it = curIdx;
            ++curIdx;
        } else {
            --n_vertices_;
        }
    }

    // Delete properties accordingly
    delete_multiple_vertex_props(_tag);

    EdgeCorrector corrector(newIndices);
    std::for_each(edges_.begin(), edges_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_edges(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_edges());

    std::vector<int> newIndices(n_edges(), -1);
    int curIdx = 0;

    std::vector<Edge> newEdges;

    std::vector<int>::iterator idx_it = newIndices.begin();
    std::vector<Edge>::const_iterator e_it = edges_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it, ++e_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newEdges.push_back(*e_it);

            *idx_it = curIdx;
            ++curIdx;
        }
    }

    // Swap edges
    edges_.swap(newEdges);

    // Delete properties accordingly
    delete_multiple_edge_props(_tag);

    FaceCorrector corrector(newIndices);
    std::for_each(faces_.begin(), faces_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_faces(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_faces());

    std::vector<int> newIndices(n_faces(), -1);
    int curIdx = 0;

    std::vector<Face> newFaces;

    std::vector<int>::iterator idx_it = newIndices.begin();
    std::vector<Face>::const_iterator f_it = faces_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++idx_it, ++f_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newFaces.push_back(*f_it);

            *idx_it = curIdx;
            ++curIdx;
        }
    }

    // Swap faces
    faces_.swap(newFaces);

    // Delete properties accordingly
    delete_multiple_face_props(_tag);

    CellCorrector corrector(newIndices);
    std::for_each(cells_.begin(), cells_.end(), corrector);
}

//========================================================================================

void TopologyKernel::delete_multiple_cells(const std::vector<bool>& _tag) {

    assert(_tag.size() == n_cells());

    std::vector<Cell> newCells;

    std::vector<Cell>::const_iterator c_it = cells_.begin();

    for(std::vector<bool>::const_iterator t_it = _tag.begin(),
            t_end = _tag.end(); t_it != t_end; ++t_it, ++c_it) {

        if(!(*t_it)) {
            // Not marked as deleted

            newCells.push_back(*c_it);
        }
    }

    // Swap cells
    cells_.swap(newCells);

    // Delete properties accordingly
    delete_multiple_cell_props(_tag);
}

//========================================================================================

// cppcheck-suppress unusedFunction ; public interface
CellIter TopologyKernel::delete_cell_range(const CellIter& _first, const CellIter& _last) {

    assert(_first >= cells_begin());
    assert(_last <= cells_end());

    std::vector<Cell>::iterator it = cells_.erase(cells_.begin() + _first->idx(), cells_.begin() + _last->idx());

    // Re-compute face bottom-up incidences if necessary
    if(has_face_bottom_up_incidences()) {
        f_bottom_up_ = false;
        enable_face_bottom_up_incidences(true);
    }

    return CellIter(this, CellHandle((int)(it - cells_.begin())));
}

void TopologyKernel::enable_deferred_deletion(bool _enable)
{
    if (deferred_deletion_ && !_enable)
        collect_garbage();

    deferred_deletion_ = _enable;
}

void TopologyKernel::compute_vertex_bottom_up_incidences()
{

}

void TopologyKernel::compute_edge_bottom_up_incidences()
{

}

void TopologyKernel::compute_face_bottom_up_incidences()
{

}

//========================================================================================

/// Get edge with handle _edgeHandle
const OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) const {

    // Test if edge is valid
    assert(_edgeHandle.is_valid() && (size_t)_edgeHandle.idx() < edges_.size());

    return edges_[_edgeHandle.idx()];
}

//========================================================================================

/// Get face with handle _faceHandle
const OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) const {

    // Test if face is valid
    assert(_faceHandle.is_valid() && (size_t)_faceHandle.idx() < faces_.size());

    return faces_[_faceHandle.idx()];
}

//========================================================================================

/// Get cell with handle _cellHandle
const OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) const {

    // Test if cell is valid
    assert(_cellHandle.is_valid() && (size_t)_cellHandle.idx() < cells_.size());

    return cells_[_cellHandle.idx()];
}

//========================================================================================

/// Get edge with handle _edgeHandle
OpenVolumeMeshEdge& TopologyKernel::edge(const EdgeHandle& _edgeHandle) {

    // Test if edge is valid
    assert(_edgeHandle.is_valid() && (size_t)_edgeHandle.idx() < edges_.size());

    return edges_[_edgeHandle.idx()];
}

//========================================================================================

/// Get face with handle _faceHandle
OpenVolumeMeshFace& TopologyKernel::face(const FaceHandle& _faceHandle) {

    // Test if face is valid
    assert((size_t)_faceHandle.idx() < faces_.size());
    assert(_faceHandle.idx() >= 0);

    return faces_[_faceHandle.idx()];
}

//========================================================================================

/// Get cell with handle _cellHandle
OpenVolumeMeshCell& TopologyKernel::cell(const CellHandle& _cellHandle) {

    // Test if cell is valid
    assert((size_t)_cellHandle.idx() < cells_.size());
    assert(_cellHandle.idx() >= 0);

    return cells_[_cellHandle.idx()];
}

//========================================================================================

/// Get edge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert((size_t)_halfEdgeHandle.idx() < (edges_.size() * 2));
    assert(_halfEdgeHandle.idx() >= 0);

    // In case the handle is even, just return the corresponding edge
    /// Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle.idx() % 2 == 0)
        return edges_[(int)(_halfEdgeHandle.idx() / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle.idx() / 2)]);
}

//========================================================================================

/// Get face that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert((size_t)_halfFaceHandle.idx() < (faces_.size() * 2));
    assert(_halfFaceHandle.idx() >= 0);

    // In case the handle is even, just return the corresponding face
    // Otherwise return the opposite halfface via opposite()
    if(_halfFaceHandle.idx() % 2 == 0)
        return faces_[(int)(_halfFaceHandle.idx() / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle.idx() / 2)]);
}

//========================================================================================

/// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::opposite_halfedge(const HalfEdgeHandle& _halfEdgeHandle) const {

    // Is handle in range?
    assert(_halfEdgeHandle.idx() >= 0);
    assert((size_t)_halfEdgeHandle.idx() < (edges_.size() * 2));

    // In case the handle is not even, just return the corresponding edge
    // Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle.idx() % 2 != 0)
        return edges_[(int)(_halfEdgeHandle.idx() / 2)];
    else
        return opposite_halfedge(edges_[(int)(_halfEdgeHandle.idx() / 2)]);
}

//========================================================================================

/// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::opposite_halfface(const HalfFaceHandle& _halfFaceHandle) const {

    // Is handle in range?
    assert(_halfFaceHandle.idx() >= 0);
    assert((size_t)_halfFaceHandle.idx() < (faces_.size() * 2));

    // In case the handle is not even, just return the corresponding face
    // Otherwise return the opposite via the first face's opposite() function
    if(_halfFaceHandle.idx() % 2 != 0)
        return faces_[(int)(_halfFaceHandle.idx() / 2)];
    else
        return opposite_halfface(faces_[(int)(_halfFaceHandle.idx() / 2)]);
}

//========================================================================================

HalfEdgeHandle TopologyKernel::halfedge(const VertexHandle& _vh1, const VertexHandle& _vh2) const {

    assert(_vh1.idx() < (int)n_vertices());
    assert(_vh2.idx() < (int)n_vertices());

    for(VertexOHalfEdgeIter voh_it = voh_iter(_vh1); voh_it.valid(); ++voh_it) {
        if(halfedge(*voh_it).to_vertex() == _vh2) {
            return *voh_it;
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

HalfFaceHandle TopologyKernel::halfface(const std::vector<VertexHandle>& _vs) const {

    assert(_vs.size() > 2);

    VertexHandle v0 = _vs[0], v1 = _vs[1], v2 = _vs[2];

    assert(v0.is_valid() && v1.is_valid() && v2.is_valid());

    HalfEdgeHandle he0 = halfedge(v0, v1);
    if(!he0.is_valid()) return InvalidHalfFaceHandle;
    HalfEdgeHandle he1 = halfedge(v1, v2);
    if(!he1.is_valid()) return InvalidHalfFaceHandle;

    std::vector<HalfEdgeHandle> hes;
    hes.push_back(he0);
    hes.push_back(he1);

    return halfface(hes);
}

HalfFaceHandle TopologyKernel::halfface_extensive(const std::vector<VertexHandle>& _vs) const
{
  //TODO: schöner machen

  assert(_vs.size() > 2);

  VertexHandle v0 = _vs[0];
  VertexHandle v1 = _vs[1];

  assert(v0.is_valid() && v1.is_valid());

  HalfEdgeHandle he0 = halfedge(v0, v1);
  if(!he0.is_valid()) return InvalidHalfFaceHandle;

  for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(he0); hehf_it.valid(); ++hehf_it)
  {
    std::vector<HalfEdgeHandle> hes = halfface(*hehf_it).halfedges();

    if (hes.size() != _vs.size())
      continue;

    int offset = 0;
    for (unsigned int i = 0; i < hes.size(); ++i)
      if (hes[i] == he0)
        offset = i;

    bool all_vertices_found = true;

    for (unsigned int i = 0; i < hes.size(); ++i)
    {
      HalfEdgeHandle heh = hes[(i+offset)%hes.size()];
      if (halfedge(heh).from_vertex() != _vs[i])
      {
        all_vertices_found = false;
        break;
      }
    }

    if (all_vertices_found)
      return *hehf_it;
  }

  return InvalidHalfFaceHandle;
}

//========================================================================================

HalfFaceHandle TopologyKernel::halfface(const std::vector<HalfEdgeHandle>& _hes) const {

    assert(_hes.size() >= 2);

    HalfEdgeHandle he0 = _hes[0], he1 = _hes[1];

    assert(he0.is_valid() && he1.is_valid());

    for(HalfEdgeHalfFaceIter hehf_it = hehf_iter(he0); hehf_it.valid(); ++hehf_it) {

        std::vector<HalfEdgeHandle> hes = halfface(*hehf_it).halfedges();
        if(std::find(hes.begin(), hes.end(), he1) != hes.end()) {
            return *hehf_it;
        }
    }

    return InvalidHalfFaceHandle;
}

//========================================================================================

HalfEdgeHandle TopologyKernel::next_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert(_heh.is_valid() && (size_t)_heh.idx() < edges_.size() * 2u);
    assert(_hfh.is_valid() && (size_t)_hfh.idx() < faces_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if((it + 1) != hes.end()) return *(it + 1);
            else return *hes.begin();
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

HalfEdgeHandle TopologyKernel::prev_halfedge_in_halfface(const HalfEdgeHandle& _heh, const HalfFaceHandle& _hfh) const {

    assert(_heh.is_valid() && (size_t)_heh.idx() < edges_.size() * 2u);
    assert(_hfh.is_valid() && (size_t)_hfh.idx() < faces_.size() * 2u);

    std::vector<HalfEdgeHandle> hes = halfface(_hfh).halfedges();

    for(std::vector<HalfEdgeHandle>::const_iterator it = hes.begin();
            it != hes.end(); ++it) {
        if(*it == _heh) {
            if(it != hes.begin()) return *(it - 1);
            else return *(hes.end() - 1);
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

HalfFaceHandle
TopologyKernel::adjacent_halfface_in_cell(const HalfFaceHandle& _halfFaceHandle, const HalfEdgeHandle& _halfEdgeHandle) const
{

    assert(_halfFaceHandle.is_valid() && (size_t)_halfFaceHandle.idx() < faces_.size() * 2u);
    assert(_halfEdgeHandle.is_valid() && (size_t)_halfEdgeHandle.idx() < edges_.size() * 2u);
    assert(has_face_bottom_up_incidences());
    assert((size_t)_halfFaceHandle.idx() < incident_cell_per_hf_.size());

    if(incident_cell_per_hf_[_halfFaceHandle.idx()] == InvalidCellHandle) {
        // Specified halfface is on the outside of the complex
        return InvalidHalfFaceHandle;
    }

    // Make sure that _halfFaceHandle is incident to _halfEdgeHandle
    bool skipped = false;
    HalfFaceHandle idx = InvalidHalfFaceHandle;

    const auto eh = edge_handle(_halfEdgeHandle);
    const auto fh = face_handle(_halfFaceHandle);
    const auto ch = incident_cell(_halfFaceHandle);

    // if we have bottom up incidences, looping through the incident faces of the halfedge should be faster.
    if (has_edge_bottom_up_incidences()) {
        for (const auto hfh: incident_hfs_per_he_[_halfEdgeHandle.idx()]) {
            if (face_handle(hfh) == fh) {
                assert(!skipped);
                if (idx.is_valid())
                    return idx;
                skipped = true;
                continue;
            }
            if (incident_cell(hfh) == ch) {
                assert(!idx.is_valid());
                if (skipped)
                    return hfh;
                idx = hfh;
            }
            const auto opp_hfh = opposite_halfface_handle(hfh);
            if (incident_cell(opp_hfh) == ch) {
                assert(!idx.is_valid());
                if (skipped)
                    return opp_hfh;
                idx = opp_hfh;
            }
        }
        return InvalidHalfFaceHandle;
    }


    for(const auto &hfh: cell(ch).halffaces()) {
        if(hfh == _halfFaceHandle) {
            assert(!skipped);
            skipped = true;
            if (idx.is_valid()) {
                return idx;
            }
        } else {
            for (const auto &heh: face(face_handle(hfh)).halfedges()) {
                if(edge_handle(heh) == eh) {
                    assert(!idx.is_valid());
                    if (skipped) {
                        return hfh;
                    } else {
                        idx = hfh;
                        continue;
                    }
                }
            }
        }
    }
    return InvalidHalfFaceHandle;
}

//========================================================================================

CellHandle TopologyKernel::incident_cell(const HalfFaceHandle& _halfFaceHandle) const {

    assert(has_face_bottom_up_incidences());
    assert((size_t)_halfFaceHandle.idx() < incident_cell_per_hf_.size() && _halfFaceHandle.idx() >= 0);

    return incident_cell_per_hf_[_halfFaceHandle.idx()];
}

//========================================================================================
//========================================================================================


} // Namespace OpenVolumeMesh
