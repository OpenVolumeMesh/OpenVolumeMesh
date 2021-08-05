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

#include "TopologyKernel.hh"

namespace OpenVolumeMesh {

// Initialize constants
const VertexHandle      TopologyKernel::InvalidVertexHandle   = VertexHandle(-1);
const EdgeHandle        TopologyKernel::InvalidEdgeHandle     = EdgeHandle(-1);
const HalfEdgeHandle    TopologyKernel::InvalidHalfEdgeHandle = HalfEdgeHandle(-1);
const FaceHandle        TopologyKernel::InvalidFaceHandle     = FaceHandle(-1);
const HalfFaceHandle    TopologyKernel::InvalidHalfFaceHandle = HalfFaceHandle(-1);
const CellHandle        TopologyKernel::InvalidCellHandle     = CellHandle(-1);

//========================================================================================
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
    if(v_bottom_up_) {
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
        if(v_bottom_up_) {

            assert((size_t)_fromVertex.idx() < outgoing_hes_per_vertex_.size());
            std::vector<HalfEdgeHandle>& ohes = outgoing_hes_per_vertex_[_fromVertex.idx()];
            for(std::vector<HalfEdgeHandle>::const_iterator he_it = ohes.begin(),
                    he_end = ohes.end(); he_it != he_end; ++he_it) {
                if(halfedge(*he_it).to_vertex() == _toVertex) {
                    return edge_handle(*he_it);
                }
            }
        } else {
            for(int i = 0; i < (int)edges_.size(); ++i) {
                if(edge(EdgeHandle(i)).from_vertex() == _fromVertex && edge(EdgeHandle(i)).to_vertex() == _toVertex) {
                    return EdgeHandle(i);
                } else if(edge(EdgeHandle(i)).from_vertex() == _toVertex && edge(EdgeHandle(i)).to_vertex() == _fromVertex) {
                    return EdgeHandle(i);
                }
            }
        }
    }

    // Store edge
    edges_.emplace_back(_fromVertex, _toVertex);
    edge_deleted_.push_back(false);

    resize_eprops(n_edges());

    EdgeHandle eh((int)edges_.size()-1);

    // Update vertex bottom-up incidences
    if(v_bottom_up_) {
        assert((size_t)_fromVertex.idx() < outgoing_hes_per_vertex_.size());
        assert((size_t)_toVertex.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[_fromVertex.idx()].push_back(halfedge_handle(eh, 0));
        outgoing_hes_per_vertex_[_toVertex.idx()].push_back(halfedge_handle(eh, 1));
    }

    // Create item for edge bottom-up incidences
    if(e_bottom_up_) {
        incident_hfs_per_he_.resize(n_halfedges());
    }

    // Get handle of recently created edge
    return eh;
}

//========================================================================================

/// Add face via incident edges
FaceHandle TopologyKernel::add_face(std::vector<HalfEdgeHandle> _halfedges, bool _topologyCheck) {

#ifndef NDEBUG
    // Assert that halfedges are valid
    for(std::vector<HalfEdgeHandle>::const_iterator it = _halfedges.begin(),
            end = _halfedges.end(); it != end; ++it)
        assert(it->is_valid() && (size_t)it->idx() < edges_.size() * 2u && !is_deleted(*it));
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

    // Update edge bottom-up incidences
    if(e_bottom_up_) {

        const auto &face_halfedges = faces_[fh.idx()].halfedges();
        for (const auto heh: face_halfedges) {
            auto opp = opposite_halfedge_handle(heh);

            assert((size_t)heh.idx() < incident_hfs_per_he_.size());
            assert((size_t)opp.idx() < incident_hfs_per_he_.size());

            incident_hfs_per_he_[heh.idx()].push_back(halfface_handle(fh, 0));
            incident_hfs_per_he_[opp.idx()].push_back(halfface_handle(fh, 1));
            // we added a single face, so the configuration is not manifold.
            // no need to call reorder_incident_halffaces(edge_handle(heh))
        }
    }

    // Create item for face bottom-up incidences
    if(f_bottom_up_) {
        incident_cell_per_hf_.resize(n_halffaces(), InvalidCellHandle);
    }

    // Return handle of recently created face
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

    /* Put halffaces in clockwise order via the
     * same cell property which now exists.
     * Note, this only works for manifold configurations though.
     * Proceed as follows: Pick one starting halfface. Assuming
     * that all halfface normals point into the incident cell,
     * we find the adjacent halfface within the incident cell
     * along the considered halfedge. We set the found halfface
     * to be the one to be processed next. If we reach an outside
     * region, we try to go back from the starting halfface in reverse
     * order. If the complex is properly connected (the pairwise
     * intersection of two adjacent 3-dimensional cells is always
     * a 2-dimensional entity, namely a facet), such an ordering
     * always exists and will be found. If not, a correct order
     * can not be given and, as a result, the related iterators
     * will address the related entities in an arbitrary fashion.
     */

    HalfEdgeHandle heh = halfedge_handle(_eh, 0);
    assert((size_t)heh.idx() < incident_hfs_per_he_.size());
    auto &incident_hfs = incident_hfs_per_he_[heh.idx()];

    const size_t n_hfs = incident_hfs.size();

    if(n_hfs < 2)
        return;

    std::vector<HalfFaceHandle> new_halffaces;
    new_halffaces.reserve(n_hfs);

    // Start with one incident halfface and go into the first direction
    auto start_hf = incident_hfs.front();
    auto cur_hf = start_hf;

    do {
        new_halffaces.push_back(cur_hf);
        if (new_halffaces.size() > incident_hfs.size()) {
            //std::cerr << "reorder_incident_halffaces(" << _eh.idx() << "): weird topology, aborting." << std::endl;
            return;
        };

        if (incident_cell(cur_hf) == InvalidCellHandle
                || is_deleted(incident_cell(cur_hf)))
            break;

        cur_hf = adjacent_halfface_in_cell(cur_hf, heh);
        if(cur_hf == InvalidHalfFaceHandle) {
            return;
        }
        cur_hf = opposite_halfface_handle(cur_hf);

    } while (cur_hf != start_hf);

    // First direction has terminated
    // If new_halffaces has the same size as old (unordered)
    // vector of incident halffaces, we are done here
    // If not, try the other way round
    // (this must be a boundary edge)
    if(new_halffaces.size() != incident_hfs.size()) {

        cur_hf = start_hf;

        while(true) {
            cur_hf = opposite_halfface_handle(cur_hf);

            if (incident_cell(cur_hf) == InvalidCellHandle
                    || is_deleted(incident_cell(cur_hf))) {
                // reached the other boundary
                break;
            }

            cur_hf = adjacent_halfface_in_cell(cur_hf, heh);
            if (cur_hf == InvalidHalfFaceHandle) {
                return;
            }

            // TODO PERF: just move everything we already have to the end *once* and fill backwards
            new_halffaces.insert(new_halffaces.begin(), cur_hf);
            if(new_halffaces.size() > incident_hfs.size()) {
                //std::cerr << "reorder_incident_halffaces(" << _eh.idx() << ") #2: weird topology, aborting" << std::endl;
                return;
            }
        }
    }

    // Everything worked just fine, set the new ordered vector
    if(new_halffaces.size() == incident_hfs.size()) {
        incident_hfs = std::move(new_halffaces);
        // update incident halffaces of the opposite halfedge:
        std::transform(incident_hfs.rbegin(), incident_hfs.rend(),
                incident_hfs_per_he_[opposite_halfedge_handle(heh).idx()].begin(),
                opposite_halfface_handle);
    }
#if 0
    else {
        std::cerr << "reorder_incident_halffaces: found " << new_halffaces.size() << " of " << incident_hfs.size()
            << " incident halffaces, likely the edge has more than one boundary! Currently not supported, not reordering." << std::endl;
        // TODO FIXME: we should support this case.
    }
#endif

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
         * Test if all halffaces are connected and form a two-manifold
         * => Cell is closed
         *
         * This test is simple: The number of involved half-edges has to be
         * exactly twice the number of involved edges.
         */

        size_t guess_n_halfedges = _halffaces.size() * valence(face_handle(_halffaces[0]));
        std::vector<HalfEdgeHandle> incidentHalfedges;
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
        auto he_end = std::unique(incidentHalfedges.begin(), incidentHalfedges.end());
        auto n_halfedges = std::distance(incidentHalfedges.begin(), he_end);
        auto e_end = std::unique(incidentHalfedges.begin(), he_end,
                [](HalfEdgeHandle a, HalfEdgeHandle b) {return a.idx()/2 == b.idx()/2;});
        auto n_edges = std::distance(incidentHalfedges.begin(), e_end);

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
    if(f_bottom_up_) {

        const auto &cell_halffaces = cells_[ch.idx()].halffaces();
        std::set<EdgeHandle> cell_edges;
        for(const auto &hfh: cell_halffaces) {
            assert((size_t)hfh.idx() < incident_cell_per_hf_.size());

#ifndef NDEBUG
            if(_topologyCheck) {
                if(incident_cell_per_hf_[hfh.idx()] != InvalidCellHandle) {
                    // Shouldn't this situation be dealt with before adding the
                    // cell and return InvalidCellHandle in this case?
                    // Mike: Not if the user intends to add non-manifold
                    // configurations. Although, in this case, he should be
                    // warned about it.
                    std::cerr << "add_cell(): One of the specified half-faces is already incident to another cell!" << std::endl;
                }
            }
#endif

            // Overwrite incident cell for current half-face
            incident_cell_per_hf_[hfh.idx()] = ch;

            // Collect all edges of cell
            for(const auto &eh: face_edges(face_handle(hfh))) {
                cell_edges.insert(eh);
            }
        }

        if(e_bottom_up_) {

            // Try to reorder all half-faces w.r.t.
            // their incident half-edges such that all
            // half-faces are in cyclic order around
            // a half-edge
            for (const auto &eh: cell_edges) {
                reorder_incident_halffaces(eh);
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

    // Update bottom-up entries
    if(has_vertex_bottom_up_incidences()) {

        const VertexHandle& fv = e.from_vertex();
        const VertexHandle& tv = e.to_vertex();

        const HalfEdgeHandle heh0 = halfedge_handle(_eh, 0);
        const HalfEdgeHandle heh1 = halfedge_handle(_eh, 1);

        std::vector<HalfEdgeHandle>::iterator h_end =
        		std::remove(outgoing_hes_per_vertex_[fv.idx()].begin(), outgoing_hes_per_vertex_[fv.idx()].end(), heh0);
        outgoing_hes_per_vertex_[fv.idx()].resize(h_end - outgoing_hes_per_vertex_[fv.idx()].begin());

        h_end = std::remove(outgoing_hes_per_vertex_[tv.idx()].begin(), outgoing_hes_per_vertex_[tv.idx()].end(), heh1);
        outgoing_hes_per_vertex_[tv.idx()].resize(h_end - outgoing_hes_per_vertex_[tv.idx()].begin());

        outgoing_hes_per_vertex_[_fromVertex.idx()].push_back(heh0);
        outgoing_hes_per_vertex_[_toVertex.idx()].push_back(heh1);
    }

    e.set_from_vertex(_fromVertex);
    e.set_to_vertex(_toVertex);
}

//========================================================================================

/// Set the half-edges of a face
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_face(const FaceHandle& _fh, const std::vector<HalfEdgeHandle>& _hes) {

    Face& f = face(_fh);

    if(has_edge_bottom_up_incidences()) {

        const HalfFaceHandle hf0 = halfface_handle(_fh, 0);
        const HalfFaceHandle hf1 = halfface_handle(_fh, 1);

        const std::vector<HalfEdgeHandle>& hes = f.halfedges();

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

        	std::vector<HalfFaceHandle>::iterator h_end =
        			std::remove(incident_hfs_per_he_[he_it->idx()].begin(),
                        		incident_hfs_per_he_[he_it->idx()].end(), hf0);
            incident_hfs_per_he_[he_it->idx()].resize(h_end - incident_hfs_per_he_[he_it->idx()].begin());

            h_end =  std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].begin(),
                        		 incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].end(), hf1);
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].resize(h_end - incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].begin());
        }

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = _hes.begin(),
                he_end = _hes.end(); he_it != he_end; ++he_it) {

            incident_hfs_per_he_[he_it->idx()].push_back(hf0);
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it).idx()].push_back(hf1);
        }

        // TODO: Reorder incident half-faces
    }

    f.set_halfedges(_hes);
}

//========================================================================================

/// Set the half-faces of a cell
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_cell(const CellHandle& _ch, const std::vector<HalfFaceHandle>& _hfs) {

    Cell& c = cell(_ch);

    if(has_face_bottom_up_incidences()) {
        for (const auto &hfh: c.halffaces()) {
            incident_cell_per_hf_[hfh.idx()] = InvalidCellHandle;
        }
        for (const auto &hfh: _hfs) {
            incident_cell_per_hf_[hfh.idx()] = _ch;
        }
    }

    c.set_halffaces(_hfs);
    // TODO: possibly reorder incident halffaces of the edges of the "new" cell
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
VertexIter TopologyKernel::delete_vertex(const VertexHandle& _h) {

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
    return delete_vertex_core(_h);
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
EdgeIter TopologyKernel::delete_edge(const EdgeHandle& _h) {

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
    return delete_edge_core(_h);
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
FaceIter TopologyKernel::delete_face(const FaceHandle& _h) {

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
    return delete_face_core(_h);
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
CellIter TopologyKernel::delete_cell(const CellHandle& _h) {

    assert(!is_deleted(_h));
    return delete_cell_core(_h);
}

/**
 * \brief Delete all entities that are marked as deleted
 */
void TopologyKernel::collect_garbage()
{
    if (!deferred_deletion_enabled() || !needs_garbage_collection())
        return; // nothing todo

    deferred_deletion = false;

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

    deferred_deletion = true;

}

//========================================================================================

template <class ContainerT>
void TopologyKernel::get_incident_edges(const ContainerT& _vs,
                                        std::set<EdgeHandle>& _es) const {

    _es.clear();

    if(v_bottom_up_) {

        for(typename ContainerT::const_iterator v_it = _vs.begin(),
                v_end = _vs.end(); v_it != v_end; ++v_it) {

            const std::vector<HalfEdgeHandle>& inc_hes = outgoing_hes_per_vertex_[v_it->idx()];

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

    if(e_bottom_up_) {

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

    if(f_bottom_up_) {

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
 * After performing this operation, all vertices
 * following vertex _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the vertex links
 * in all edges. These steps are performed:
 *
 * 1) Decrease all vertex handles > _h in incident edges
 * 2) Delete entry in bottom-up list: V -> HE
 * 3) Delete vertex itself (not necessary here since
 *    a vertex is only represented by a number)
 * 4) Delete property entry
 *
 * @param _h A vertex's handle
 */
VertexIter TopologyKernel::delete_vertex_core(const VertexHandle& _h) {

    VertexHandle h = _h;
    assert(h.is_valid() && (size_t)h.idx() < n_vertices());

    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last not deleted vertex
    {
        VertexHandle last_undeleted_vertex = VertexHandle((int)n_vertices()-1);
        assert(!vertex_deleted_[last_undeleted_vertex.idx()]);
        swap_vertex_indices(h, last_undeleted_vertex);
        h = last_undeleted_vertex;
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_vertices_;
        vertex_deleted_[h.idx()] = true;
//        deleted_vertices_.push_back(h);

        // Iterator to next element in vertex list
//        return (vertices_begin() + h.idx()+1);
        return VertexIter(this, VertexHandle(h.idx()+1));
    }
    else
    {
        // 1)
        if(v_bottom_up_) {

            // Decrease all vertex handles >= _h in all edge definitions
            for(int i = h.idx(), end = (int)n_vertices(); i < end; ++i) {
                const std::vector<HalfEdgeHandle>& hes = outgoing_hes_per_vertex_[i];
                for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                    he_end = hes.end(); he_it != he_end; ++he_it) {

                    Edge& e = edge(edge_handle(*he_it));
                    if(e.from_vertex().idx() == i) {
                        e.set_from_vertex(VertexHandle(i-1));
                    }
                    if(e.to_vertex().idx() == i) {
                        e.set_to_vertex(VertexHandle(i-1));
                    }
                }
            }

        } else {

            // Iterate over all edges
            for(EdgeIter e_it = edges_begin(), e_end = edges_end();
                e_it != e_end; ++e_it) {

                // Decrease all vertex handles in edge definitions that are greater than _h
                if(edge(*e_it).from_vertex() > h) {
                    edge(*e_it).set_from_vertex(VertexHandle(edge(*e_it).from_vertex().idx() - 1));
                }
                if(edge(*e_it).to_vertex() > h) {
                    edge(*e_it).set_to_vertex(VertexHandle(edge(*e_it).to_vertex().idx() - 1));
                }
            }
        }

        // 2)

        if(v_bottom_up_) {
            assert((size_t)h.idx() < outgoing_hes_per_vertex_.size());
            outgoing_hes_per_vertex_.erase(outgoing_hes_per_vertex_.begin() + h.idx());
        }


        // 3)

        --n_vertices_;
        vertex_deleted_.erase(vertex_deleted_.begin() + h.idx());

        // 4)

        vertex_deleted(h);

        // Iterator to next element in vertex list
//        return (vertices_begin() + h.idx());
        return VertexIter(this, h);

    }
}

//========================================================================================

/**
 * \brief Delete edge from mesh
 *
 * After performing this operation, all edges
 * following edge _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the edge links
 * in all faces. These steps are performed:
 *
 * 1) Delete bottom-up links from incident vertices
 * 2) Decrease all half-edge handles > _h in incident faces
 * 3) Delete entry in bottom-up list: HE -> HF
 * 4) Decrease all half-edge handles > 2*_h.idx() in
 *    vertex bottom-up list
 * 5) Delete edge itself
 * 6) Delete property entry
 *
 * @param _h An edge's handle
 */
EdgeIter TopologyKernel::delete_edge_core(const EdgeHandle& _h) {

    EdgeHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < edges_.size());

    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        EdgeHandle last_edge = EdgeHandle((int)edges_.size()-1);
        assert(!edge_deleted_[last_edge.idx()]);
        swap_edge_indices(h, last_edge);
        h = last_edge;
    }


    // 1)
    if(v_bottom_up_) {

        VertexHandle v0 = edge(h).from_vertex();
        VertexHandle v1 = edge(h).to_vertex();
        assert(v0.is_valid() && (size_t)v0.idx() < outgoing_hes_per_vertex_.size());
        assert(v1.is_valid() && (size_t)v1.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[v0.idx()].erase(
                std::remove(outgoing_hes_per_vertex_[v0.idx()].begin(),
                            outgoing_hes_per_vertex_[v0.idx()].end(),
                            halfedge_handle(h, 0)),
                            outgoing_hes_per_vertex_[v0.idx()].end());

        outgoing_hes_per_vertex_[v1.idx()].erase(
                std::remove(outgoing_hes_per_vertex_[v1.idx()].begin(),
                            outgoing_hes_per_vertex_[v1.idx()].end(),
                            halfedge_handle(h, 1)),
                            outgoing_hes_per_vertex_[v1.idx()].end());
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_edges_;
        edge_deleted_[h.idx()] = true;
//        deleted_edges_.push_back(h);

        // Return iterator to next element in list
//        return (edges_begin() + h.idx()+1);
        return EdgeIter(this, EdgeHandle(h.idx()+1));
    }
    else
    {

        if (!fast_deletion_enabled())
        {
            // 2)
            if(e_bottom_up_) {

                assert((size_t)halfedge_handle(h, 0).idx() < incident_hfs_per_he_.size());

                // Decrease all half-edge handles > he and
                // delete all half-edge handles == he in face definitions
                // Get all faces that need updates
                std::set<FaceHandle> update_faces;
                for(std::vector<std::vector<HalfFaceHandle> >::const_iterator iit =
                    (incident_hfs_per_he_.begin() + halfedge_handle(h, 0).idx()),
                    iit_end = incident_hfs_per_he_.end(); iit != iit_end; ++iit) {
                    for(std::vector<HalfFaceHandle>::const_iterator it = iit->begin(),
                        end = iit->end(); it != end; ++it) {
                        update_faces.insert(face_handle(*it));
                    }
                }
                // Update respective handles
                HEHandleCorrection cor(halfedge_handle(h, 1));
                for(std::set<FaceHandle>::iterator f_it = update_faces.begin(),
                    f_end = update_faces.end(); f_it != f_end; ++f_it) {

                    std::vector<HalfEdgeHandle> hes = face(*f_it).halfedges();

                    // Delete current half-edge from face's half-edge list
                    hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(h, 0)), hes.end());
                    hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(h, 1)), hes.end());

                    cor.correctVecValue(hes);

                    face(*f_it).set_halfedges(hes);
                }
            } else {

                // Iterate over all faces
                for(FaceIter f_it = faces_begin(), f_end = faces_end();
                    f_it != f_end; ++f_it) {

                    // Get face's half-edges
                    std::vector<HalfEdgeHandle> hes = face(*f_it).halfedges();

                    // Delete current half-edge from face's half-edge list
                    hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(h, 0)), hes.end());
                    hes.erase(std::remove(hes.begin(), hes.end(), halfedge_handle(h, 1)), hes.end());

                    // Decrease all half-edge handles greater than _h in face
                    HEHandleCorrection cor(halfedge_handle(h, 1));
                    cor.correctVecValue(hes);
                    face(*f_it).set_halfedges(hes);
                }
            }
        }

        // 3)

        if(e_bottom_up_) {
            assert((size_t)halfedge_handle(h, 1).idx() < incident_hfs_per_he_.size());

            incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 1).idx());
            incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 0).idx());
        }

        if (!fast_deletion_enabled())
        {
            // 4)
            if(v_bottom_up_) {
                HEHandleCorrection cor(halfedge_handle(h, 1));
                for (auto &ohehs: outgoing_hes_per_vertex_) {
                    cor.correctVecValue(ohehs);
                }
            }
        }


        // 5)
        edges_.erase(edges_.begin() + h.idx());
        edge_deleted_.erase(edge_deleted_.begin() + h.idx());


        // 6)

        edge_deleted(h);

        // Return iterator to next element in list
//        return (edges_begin() + h.idx());
        return EdgeIter(this, h);

    }
}

//========================================================================================

/**
 * \brief Delete face from mesh
 *
 * After performing this operation, all faces
 * following face _h in the array will be accessible
 * through their old handle decreased by one.
 * This function directly fixes the face links
 * in all cells. These steps are performed:
 *
 * 1) Delete bottom-up links from incident edges
 * 2) Decrease all half-face handles > _h in incident cells
 * 3) Delete entry in bottom-up list: HF -> C
 * 4) Decrease all half-face handles > 2*_h.idx() in
 *    half-edge bottom-up list
 * 5) Delete face itself
 * 6) Delete property entry
 *
 * @param _h An face's handle
 */
FaceIter TopologyKernel::delete_face_core(const FaceHandle& _h) {

    FaceHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < faces_.size());


    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        FaceHandle last_face = FaceHandle((int)faces_.size()-1);
        assert(!face_deleted_[last_face.idx()]);
        swap_face_indices(h, last_face);
        h = last_face;
    }

    // 1)
    if(e_bottom_up_) {

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
//        deleted_faces_.push_back(h);

        // Return iterator to next element in list
//        return (faces_begin() + h.idx()+1);
        return FaceIter(this, FaceHandle(h.idx()+1));
    }
    else
    {

        if (!fast_deletion_enabled())
        {
            // 2)
            if(f_bottom_up_) {

                // Decrease all half-face handles > _h in all cells
                // and delete all half-face handles == _h
                std::set<CellHandle> update_cells;
                for(std::vector<CellHandle>::const_iterator c_it = (incident_cell_per_hf_.begin() + halfface_handle(h, 0).idx()),
                    c_end = incident_cell_per_hf_.end(); c_it != c_end; ++c_it) {
                    if(!c_it->is_valid()) continue;
                    update_cells.insert(*c_it);
                }
                for(std::set<CellHandle>::const_iterator c_it = update_cells.begin(),
                    c_end = update_cells.end(); c_it != c_end; ++c_it) {

                    std::vector<HalfFaceHandle> hfs = cell(*c_it).halffaces();

                    // Delete current half-faces from cell's half-face list
                    hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(h, 0)), hfs.end());
                    hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(h, 1)), hfs.end());

                    HFHandleCorrection cor(halfface_handle(h, 1));
                    cor.correctVecValue(hfs);

                    cell(*c_it).set_halffaces(hfs);
                }

            } else {

                // Iterate over all cells
                for(CellIter c_it = cells_begin(), c_end = cells_end(); c_it != c_end; ++c_it) {

                    std::vector<HalfFaceHandle> hfs = cell(*c_it).halffaces();

                    // Delete current half-faces from cell's half-face list
                    hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(h, 0)), hfs.end());
                    hfs.erase(std::remove(hfs.begin(), hfs.end(), halfface_handle(h, 1)), hfs.end());

                    HFHandleCorrection cor(halfface_handle(h, 1));
                    cor.correctVecValue(hfs);
                    cell(*c_it).set_halffaces(hfs);
                }
            }
        }


        // 3)
        if(f_bottom_up_) {
            assert((size_t)halfface_handle(h, 1).idx() < incident_cell_per_hf_.size());

            incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 1).idx());
            incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 0).idx());
        }


        if (!fast_deletion_enabled())
        {
            // 4)
            if(e_bottom_up_) {
                HFHandleCorrection cor(halfface_handle(h, 1));
                for (auto &ihfs: incident_hfs_per_he_) {
                    cor.correctVecValue(ihfs);
                }
            }
        }

        // 5)
        faces_.erase(faces_.begin() + h.idx());
        face_deleted_.erase(face_deleted_.begin() + h.idx());

        // 6)
        face_deleted(h);

        // Return iterator to next element in list
//        return (faces_begin() + h.idx());
        return FaceIter(this, h);
    }

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
CellIter TopologyKernel::delete_cell_core(const CellHandle& _h) {

    CellHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < cells_.size());


    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last not deleted cell
    {
        CellHandle last_undeleted_cell = CellHandle((int)cells_.size()-1);
        assert(!cell_deleted_[last_undeleted_cell.idx()]);
        swap_cell_indices(h, last_undeleted_cell);
        h = last_undeleted_cell;
    }


    // 1)
    if(f_bottom_up_) {
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
//        deleted_cells_.push_back(h);
//        deleted_cells_set.insert(h);

//        return (cells_begin() + h.idx()+1);
        return CellIter(this, CellHandle(h.idx()+1));
    }
    else
    {
        // 2)
        if (!fast_deletion_enabled())
        {
            if(f_bottom_up_) {
                CHandleCorrection cor(h);
                cor.correctVecValue(incident_cell_per_hf_);
            }
        }

        // 3)
        cells_.erase(cells_.begin() + h.idx());
        cell_deleted_.erase(cell_deleted_.begin() + h.idx());

        // 4)
        cell_deleted(h);

        // return handle to original position
//        return (cells_begin() + h.idx()+1);
        return CellIter(this, h);
    }
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

    std::vector<unsigned int> ids;
    ids.push_back(_h1.idx());
    ids.push_back(_h2.idx());


    // correct pointers to those vertices

    if (has_vertex_bottom_up_incidences())
    {
        for (unsigned int i = 0; i < 2; ++i) // For both swapped vertices
        {
            std::set<unsigned int> processed_edges; // to ensure ids are only swapped once (in the case that the two swapped vertices are connected by an edge)
            std::vector<HalfEdgeHandle>& outgoing_hes = outgoing_hes_per_vertex_[ids[i]];
            for (unsigned int k = 0; k < outgoing_hes.size(); ++k) // for each outgoing halfedge
            {
                unsigned int e_id = outgoing_hes[k].idx() / 2;

                if (processed_edges.find(e_id) == processed_edges.end())
                {
                    Edge& e = edges_[e_id];
                    if (e.from_vertex().idx() == (int)ids[0])
                        e.set_from_vertex(VertexHandle(ids[1]));
                    else if (e.from_vertex().idx() == (int)ids[1])
                        e.set_from_vertex(VertexHandle(ids[0]));

                    if (e.to_vertex().idx() == (int)ids[0])
                        e.set_to_vertex(VertexHandle(ids[1]));
                    else if (e.to_vertex().idx() == (int)ids[1])
                        e.set_to_vertex(VertexHandle(ids[0]));

                    processed_edges.insert(e_id);
                }
            }
        }

    }
    else
    {
        // search for all edges containing a swapped vertex

        for (unsigned int i = 0; i < edges_.size(); ++i)
        {
            Edge& e = edges_[i];
            if (e.from_vertex().idx() == (int)ids[0])
                e.set_from_vertex(VertexHandle(ids[1]));
            else if (e.from_vertex().idx() == (int)ids[1])
                e.set_from_vertex(VertexHandle(ids[0]));

            if (e.to_vertex().idx() == (int)ids[0])
                e.set_to_vertex(VertexHandle(ids[1]));
            else if (e.to_vertex().idx() == (int)ids[1])
                e.set_to_vertex(VertexHandle(ids[0]));
        }
    }

    // swap vector entries
    bool tmp = vertex_deleted_[ids[0]];
    vertex_deleted_[ids[0]] = vertex_deleted_[ids[1]];
    vertex_deleted_[ids[1]] = tmp;
    std::swap(outgoing_hes_per_vertex_[ids[0]], outgoing_hes_per_vertex_[ids[1]]);
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
    if(f_bottom_up_) {
        f_bottom_up_ = false;
        enable_face_bottom_up_incidences(true);
    }

    return CellIter(this, CellHandle((int)(it - cells_.begin())));
}

void TopologyKernel::enable_deferred_deletion(bool _enable)
{
    if (deferred_deletion && !_enable)
        collect_garbage();

    deferred_deletion = _enable;
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

void TopologyKernel::compute_vertex_bottom_up_incidences() {

    // Clear incidences
    outgoing_hes_per_vertex_.clear();
    outgoing_hes_per_vertex_.resize(n_vertices());
    std::vector<int> n_edges_per_vertex(n_vertices(), 0);
    for (const auto &eh: edges()) {
        ++n_edges_per_vertex[edges_[eh.idx()].from_vertex().idx()];
        ++n_edges_per_vertex[edges_[eh.idx()].to_vertex().idx()];
    }
    for (const auto &vh: vertices()) {
        outgoing_hes_per_vertex_[vh.idx()].reserve(n_edges_per_vertex[vh.idx()]);
    }

    // Store outgoing halfedges per vertex
    int n_edges = (int)edges_.size();
    for(int i = 0; i < n_edges; ++i) {
        if (edge_deleted_[i])
            continue;

        VertexHandle from = edges_[i].from_vertex();
        // If this condition is not fulfilled, it is out of caller's control and
        // definitely our bug, therefore an assert
        assert((size_t)from.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[from.idx()].push_back(halfedge_handle(EdgeHandle(i), 0));

        VertexHandle to = edges_[i].to_vertex();
        assert((size_t)to.idx() < outgoing_hes_per_vertex_.size());

        // Store opposite halfedge handle
        outgoing_hes_per_vertex_[to.idx()].push_back(halfedge_handle(EdgeHandle(i), 1));
    }
}

//========================================================================================

void TopologyKernel::compute_edge_bottom_up_incidences() {

    // Clear
    incident_hfs_per_he_.resize(n_halfedges());

    std::vector<int> n_faces_per_edge(n_edges(), 0);
    for (const auto &fh: faces()) {
        for (const auto &heh: face(fh).halfedges()) {
            ++n_faces_per_edge[edge_handle(heh).idx()];
        }
    }
    for (const auto &eh: edges()) {
        incident_hfs_per_he_[halfedge_handle(eh, 0).idx()].clear();
        incident_hfs_per_he_[halfedge_handle(eh, 0).idx()].reserve(n_faces_per_edge[eh.idx()]);
        incident_hfs_per_he_[halfedge_handle(eh, 1).idx()].clear();
        incident_hfs_per_he_[halfedge_handle(eh, 1).idx()].reserve(n_faces_per_edge[eh.idx()]);
    }
    // Store incident halffaces per halfedge
    for (const auto &fh: faces()) {
        for(const auto &heh: faces_[fh.idx()].halfedges()) {
            auto opp = opposite_halfedge_handle(heh);
            incident_hfs_per_he_[heh.idx()].push_back(halfface_handle(fh, 0));
            incident_hfs_per_he_[opp.idx()].push_back(halfface_handle(fh, 1));
        }
    }
}

//========================================================================================

void TopologyKernel::compute_face_bottom_up_incidences() {

    // Clear
    incident_cell_per_hf_.clear();
    incident_cell_per_hf_.resize(faces_.size() * 2u, InvalidCellHandle);

    int n_cells = (int)cells_.size();
    for(int i = 0; i < n_cells; ++i) {
        if (cell_deleted_[i])
            continue;

        std::vector<HalfFaceHandle> halffaces = cells_[i].halffaces();

        // Go over all halffaces
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = halffaces.begin();
                hf_it != halffaces.end(); ++hf_it) {

            if(incident_cell_per_hf_[hf_it->idx()] == InvalidCellHandle) {

                incident_cell_per_hf_[hf_it->idx()] = CellHandle(i);

            } else {

#ifndef NDEBUG
                std::cerr << "compute_face_bottom_up_incidences(): Detected non-three-manifold configuration!" << std::endl;
                std::cerr << "Connectivity probably won't work." << std::endl;
#endif
                continue;
            }
        }
    }
}

} // Namespace OpenVolumeMesh
