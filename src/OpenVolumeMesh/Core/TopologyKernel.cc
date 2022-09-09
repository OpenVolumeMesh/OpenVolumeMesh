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

#include <queue>

#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Core/detail/swap_bool.hh>

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

void TopologyKernel::add_n_vertices(size_t n)
{
    resize_vprops(n_vertices_ + n);
    n_vertices_ += n;
    vertex_deleted_.resize(n_vertices_, false);
    if(has_vertex_bottom_up_incidences()) {
        outgoing_hes_per_vertex_.resize(n_vertices_);
    }
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
EdgeHandle TopologyKernel::add_edge(VertexHandle _fromVertex,
                                    VertexHandle _toVertex,
                                    bool _allowDuplicates) {

    // If the conditions are not fulfilled, assert will fail (instead
	// of returning an invalid handle)
    assert(_fromVertex.is_valid() && (size_t)_fromVertex.idx() < n_vertices() && !is_deleted(_fromVertex));
    assert(_toVertex.is_valid() && (size_t)_toVertex.idx() < n_vertices() && !is_deleted(_toVertex));

    // Test if edge does not exist, yet
    if(!_allowDuplicates) {
        if(has_vertex_bottom_up_incidences()) {

            assert((size_t)_fromVertex.idx() < outgoing_hes_per_vertex_.size());
            std::vector<HalfEdgeHandle>& ohes = outgoing_hes_per_vertex_[_fromVertex];
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
    if(has_vertex_bottom_up_incidences()) {
        assert((size_t)_fromVertex.idx() < outgoing_hes_per_vertex_.size());
        assert((size_t)_toVertex.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[_fromVertex].push_back(halfedge_handle(eh, 0));
        outgoing_hes_per_vertex_[_toVertex].push_back(halfedge_handle(eh, 1));
    }

    // Create item for edge bottom-up incidences
    if(has_edge_bottom_up_incidences()) {
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
    if(has_edge_bottom_up_incidences()) {

        for (const auto heh: face_halfedges(fh)) {
            auto opp = opposite_halfedge_handle(heh);

            assert((size_t)heh.idx() < incident_hfs_per_he_.size());
            assert((size_t)opp.idx() < incident_hfs_per_he_.size());

            incident_hfs_per_he_[heh].push_back(halfface_handle(fh, 0));
            incident_hfs_per_he_[opp].push_back(halfface_handle(fh, 1));
        }
    }

    // Create item for face bottom-up incidences
    if(has_face_bottom_up_incidences()) {
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

void TopologyKernel::reorder_incident_halffaces(EdgeHandle _eh) {

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
    auto &incident_hfs = incident_hfs_per_he_[heh];

    const size_t n_hfs = incident_hfs.size();

    if(n_hfs < 2)
        return;

    std::vector<HalfFaceHandle> new_halffaces;
    new_halffaces.reserve(n_hfs);

    // Start with one incident halfface and go into the first direction
    auto start_hf = incident_hfs.front();
    auto cur_hf = start_hf;
    auto cur_heh = heh;

    do {
        new_halffaces.push_back(cur_hf);
        if (new_halffaces.size() > incident_hfs.size()) {
            //std::cerr << "reorder_incident_halffaces(" << _eh.idx() << "): weird topology, aborting." << std::endl;
            return;
        };

        if (incident_cell(cur_hf) == InvalidCellHandle
                || is_deleted(incident_cell(cur_hf)))
            break;

        cur_hf = adjacent_halfface_in_cell(cur_hf, cur_heh);
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
        cur_heh = opposite_halfedge_handle(heh);

        while(true) {
            cur_hf = opposite_halfface_handle(cur_hf);

            if (incident_cell(cur_hf) == InvalidCellHandle
                    || is_deleted(incident_cell(cur_hf))) {
                // reached the other boundary
                break;
            }

            cur_hf = adjacent_halfface_in_cell(cur_hf, cur_heh);
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
                incident_hfs_per_he_[opposite_halfedge_handle(heh)].begin(),
                opposite_halfface_handle);
    }
#if 0
    else {
        std::cerr << "reorder_incident_halffaces: found " << new_halffaces.size() << " of " << incident_hfs.size()
            << " incident halffaces, likely the edge has more than one boundary! Currently not supported, not reordering." << std::endl;
        // TODO FIXME: we should support this case.
    }
#endif

}//========================================================================================

/// Add cell via incident halffaces
CellHandle TopologyKernel::add_cell(std::vector<HalfFaceHandle> _halffaces, bool _topologyCheck) {

#ifndef NDEBUG
    // Assert that halffaces have valid indices
    for(std::vector<HalfFaceHandle>::const_iterator it = _halffaces.begin(),
            end = _halffaces.end(); it != end; ++it)
        assert(it->is_valid() && ((size_t)it->idx() < faces_.size() * 2u) && !is_deleted(*it));
#endif


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
    }

    // Perform topology chec
    // Create new cell
    cells_.emplace_back(std::move(_halffaces));
    cell_deleted_.push_back(false);

    // Resize props
    resize_cprops(n_cells());

    CellHandle ch((int)cells_.size()-1);

    // Update face bottom-up incidences
    if(has_face_bottom_up_incidences()) {

        const auto &cell_halffaces = cells_[ch].halffaces();
        std::set<EdgeHandle> cell_edges;
        for(const auto &hfh: cell_halffaces) {
            assert((size_t)hfh.idx() < incident_cell_per_hf_.size());

#ifndef NDEBUG
            if(_topologyCheck) {
                if(incident_cell_per_hf_[hfh] != InvalidCellHandle) {
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
            incident_cell_per_hf_[hfh] = ch;

            // Collect all edges of cell
            for(const auto &eh: face_edges(face_handle(hfh))) {
                cell_edges.insert(eh);
            }
        }

        if(has_edge_bottom_up_incidences()) {

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
void TopologyKernel::set_edge(EdgeHandle _eh, VertexHandle _fromVertex, VertexHandle _toVertex) {

    assert(_fromVertex.is_valid() && (size_t)_fromVertex.idx() < n_vertices() && !is_deleted(_fromVertex));
    assert(_toVertex.is_valid() && (size_t)_toVertex.idx() < n_vertices() && !is_deleted(_toVertex));

    Edge& e = edge(_eh);

    // Update bottom-up entries
    if(has_vertex_bottom_up_incidences()) {

        VertexHandle fv = e.from_vertex();
        VertexHandle tv = e.to_vertex();

        const HalfEdgeHandle heh0 = halfedge_handle(_eh, 0);
        const HalfEdgeHandle heh1 = halfedge_handle(_eh, 1);

        std::vector<HalfEdgeHandle>::iterator h_end =
        		std::remove(outgoing_hes_per_vertex_[fv].begin(), outgoing_hes_per_vertex_[fv].end(), heh0);
        outgoing_hes_per_vertex_[fv].resize(h_end - outgoing_hes_per_vertex_[fv].begin());

        h_end = std::remove(outgoing_hes_per_vertex_[tv].begin(), outgoing_hes_per_vertex_[tv].end(), heh1);
        outgoing_hes_per_vertex_[tv].resize(h_end - outgoing_hes_per_vertex_[tv].begin());

        outgoing_hes_per_vertex_[_fromVertex].push_back(heh0);
        outgoing_hes_per_vertex_[_toVertex].push_back(heh1);
    }

    e.set_from_vertex(_fromVertex);
    e.set_to_vertex(_toVertex);
}

//========================================================================================

/// Set the half-edges of a face
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_face(FaceHandle _fh, const std::vector<HalfEdgeHandle>& _hes) {

    Face& f = face(_fh);

    if(has_edge_bottom_up_incidences()) {

        const HalfFaceHandle hf0 = halfface_handle(_fh, 0);
        const HalfFaceHandle hf1 = halfface_handle(_fh, 1);

        const std::vector<HalfEdgeHandle>& hes = f.halfedges();

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

        	std::vector<HalfFaceHandle>::iterator h_end =
        			std::remove(incident_hfs_per_he_[*he_it].begin(),
                        		incident_hfs_per_he_[*he_it].end(), hf0);
            incident_hfs_per_he_[*he_it].resize(h_end - incident_hfs_per_he_[*he_it].begin());

            h_end =  std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].begin(),
                        		 incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].end(), hf1);
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].resize(h_end - incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].begin());
        }

        for(std::vector<HalfEdgeHandle>::const_iterator he_it = _hes.begin(),
                he_end = _hes.end(); he_it != he_end; ++he_it) {

            incident_hfs_per_he_[*he_it].push_back(hf0);
            incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].push_back(hf1);
        }

        // TODO: Reorder incident half-faces
    }

    f.set_halfedges(_hes);
}

//========================================================================================

/// Set the half-faces of a cell
// cppcheck-suppress unusedFunction ; public interface
void TopologyKernel::set_cell(CellHandle _ch, const std::vector<HalfFaceHandle>& _hfs) {

    Cell& c = cell(_ch);

    if(has_face_bottom_up_incidences()) {

        const std::vector<HalfFaceHandle>& hfs = c.halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {

            incident_cell_per_hf_[*hf_it] = InvalidCellHandle;
        }

        for(std::vector<HalfFaceHandle>::const_iterator hf_it = _hfs.begin(),
                hf_end = _hfs.end(); hf_it != hf_end; ++hf_it) {

            incident_cell_per_hf_[*hf_it] = _ch;
        }
    }

    c.set_halffaces(_hfs);
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
VertexIter TopologyKernel::delete_vertex(VertexHandle _h) {

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
EdgeIter TopologyKernel::delete_edge(EdgeHandle _h) {

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
FaceIter TopologyKernel::delete_face(FaceHandle _h) {

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
CellIter TopologyKernel::delete_cell(CellHandle _h) {

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

    deferred_deletion_ = false;

    for (int i = (int)n_cells(); i > 0; --i) {
        auto ch = CellHandle(i - 1);
        if (is_deleted(ch)) {
            cell_deleted_[ch] = false;
            delete_cell_core(ch);
        }
    }
    n_deleted_cells_ = 0;

    for (int i = (int)n_faces(); i > 0; --i) {
        auto fh = FaceHandle(i - 1);
        if (is_deleted(fh)) {
            face_deleted_[fh] = false;
            delete_face_core(fh);
        }
    }
    n_deleted_faces_ = 0;

    for (int i = (int)n_edges(); i > 0; --i) {
        auto eh = EdgeHandle(i - 1);
        if (is_deleted(eh)) {
            edge_deleted_[eh] = false;
            delete_edge_core(eh);
        }
    }
    n_deleted_edges_ = 0;

    for (int i = (int)n_vertices(); i > 0; --i) {
        auto vh = VertexHandle(i - 1);
        if (is_deleted(vh)) {
            vertex_deleted_[vh] = false;
            delete_vertex_core(vh);
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

            const std::vector<HalfEdgeHandle>& inc_hes = outgoing_hes_per_vertex_[*v_it];

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
VertexIter TopologyKernel::delete_vertex_core(VertexHandle _h) {

    VertexHandle h = _h;
    assert(h.is_valid() && (size_t)h.idx() < n_vertices());

    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last not deleted vertex
    {
        VertexHandle last_undeleted_vertex = VertexHandle((int)n_vertices()-1);
        assert(!vertex_deleted_[last_undeleted_vertex]);
        swap_vertex_indices(h, last_undeleted_vertex);
        h = last_undeleted_vertex;
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_vertices_;
        vertex_deleted_[h] = true;
//        deleted_vertices_.push_back(h);

        // Iterator to next element in vertex list
//        return (vertices_begin() + h.idx()+1);
        return VertexIter(this, VertexHandle(h.idx()+1));
    }
    else
    {
        // 1)
        if(has_vertex_bottom_up_incidences()) {

            // Decrease all vertex handles >= _h in all edge definitions
            for(int i = h.idx(), end = (int)n_vertices(); i < end; ++i) {
                const std::vector<HalfEdgeHandle>& hes = outgoing_hes_per_vertex_[VertexHandle(i)];
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

        if(has_vertex_bottom_up_incidences()) {
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
EdgeIter TopologyKernel::delete_edge_core(EdgeHandle _h) {

    EdgeHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < edges_.size());

    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        EdgeHandle last_edge = EdgeHandle((int)edges_.size()-1);
        assert(!edge_deleted_[last_edge]);
        swap_edge_indices(h, last_edge);
        h = last_edge;
    }


    // 1)
    if(has_vertex_bottom_up_incidences()) {

        VertexHandle v0 = edge(h).from_vertex();
        VertexHandle v1 = edge(h).to_vertex();
        assert(v0.is_valid() && (size_t)v0.idx() < outgoing_hes_per_vertex_.size());
        assert(v1.is_valid() && (size_t)v1.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[v0].erase(
                std::remove(outgoing_hes_per_vertex_[v0].begin(),
                            outgoing_hes_per_vertex_[v0].end(),
                            halfedge_handle(h, 0)),
                            outgoing_hes_per_vertex_[v0].end());

        outgoing_hes_per_vertex_[v1].erase(
                std::remove(outgoing_hes_per_vertex_[v1].begin(),
                            outgoing_hes_per_vertex_[v1].end(),
                            halfedge_handle(h, 1)),
                            outgoing_hes_per_vertex_[v1].end());
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_edges_;
        edge_deleted_[h] = true;
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
            if(has_edge_bottom_up_incidences()) {

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

    #if defined(__clang_major__) && (__clang_major__ >= 5)
                    for(std::vector<HalfEdgeHandle>::iterator it = hes.begin(), end = hes.end();
                        it != end; ++it) {
                        cor.correctValue(*it);
                    }
    #else
                    std::for_each(hes.begin(), hes.end(),
                                  std::bind(&HEHandleCorrection::correctValue, &cor, std::placeholders::_1));
    #endif
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
    #if defined(__clang_major__) && (__clang_major__ >= 5)
                    for(std::vector<HalfEdgeHandle>::iterator it = hes.begin(), end = hes.end();
                        it != end; ++it) {
                        cor.correctValue(*it);
                    }
    #else
                    std::for_each(hes.begin(), hes.end(),
                                  std::bind(&HEHandleCorrection::correctValue, &cor, std::placeholders::_1));
    #endif
                    face(*f_it).set_halfedges(hes);
                }
            }
        }

        // 3)

        if(has_edge_bottom_up_incidences()) {
            assert((size_t)halfedge_handle(h, 1).idx() < incident_hfs_per_he_.size());

            incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 1).idx());
            incident_hfs_per_he_.erase(incident_hfs_per_he_.begin() + halfedge_handle(h, 0).idx());
        }

        if (!fast_deletion_enabled())
        {
            // 4)
            if(has_vertex_bottom_up_incidences()) {
                HEHandleCorrection cor(halfedge_handle(h, 1));
    #if defined(__clang_major__) && (__clang_major__ >= 5)
                for(std::vector<std::vector<HalfEdgeHandle> >::iterator it = outgoing_hes_per_vertex_.begin(),
                    end = outgoing_hes_per_vertex_.end(); it != end; ++it) {
                    cor.correctVecValue(*it);
                }
    #else
                std::for_each(outgoing_hes_per_vertex_.begin(),
                              outgoing_hes_per_vertex_.end(),
                              std::bind(&HEHandleCorrection::correctVecValue, &cor, std::placeholders::_1));
    #endif
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
FaceIter TopologyKernel::delete_face_core(FaceHandle _h) {

    FaceHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < faces_.size());


    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last one
    {
        FaceHandle last_face = FaceHandle((int)faces_.size()-1);
        assert(!face_deleted_[last_face]);
        swap_face_indices(h, last_face);
        h = last_face;
    }

    // 1)
    if(has_edge_bottom_up_incidences()) {

        const std::vector<HalfEdgeHandle>& hes = face(h).halfedges();
        for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin(),
                he_end = hes.end(); he_it != he_end; ++he_it) {

            assert((size_t)std::max(he_it->idx(), opposite_halfedge_handle(*he_it).idx()) < incident_hfs_per_he_.size());

            incident_hfs_per_he_[*he_it].erase(
                    std::remove(incident_hfs_per_he_[*he_it].begin(),
                                incident_hfs_per_he_[*he_it].end(),
                                halfface_handle(h, 0)), incident_hfs_per_he_[*he_it].end());


            incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].erase(
                    std::remove(incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].begin(),
                                incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].end(),
                                halfface_handle(h, 1)), incident_hfs_per_he_[opposite_halfedge_handle(*he_it)].end());

            reorder_incident_halffaces(edge_handle(*he_it));
        }
    }

    if (deferred_deletion_enabled())
    {
        ++n_deleted_faces_;
        face_deleted_[h] = true;
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
            if(has_face_bottom_up_incidences()) {

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
#if defined(__clang_major__) && (__clang_major__ >= 5)
                    for(std::vector<HalfFaceHandle>::iterator it = hfs.begin(),
                        end = hfs.end(); it != end; ++it) {
                        cor.correctValue(*it);
                    }
#else
                    std::for_each(hfs.begin(), hfs.end(),
                                  std::bind(&HFHandleCorrection::correctValue, &cor, std::placeholders::_1));
#endif
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
#if defined(__clang_major__) && (__clang_major__ >= 5)
                    for(std::vector<HalfFaceHandle>::iterator it = hfs.begin(),
                        end = hfs.end(); it != end; ++it) {
                        cor.correctValue(*it);
                    }
#else
                    std::for_each(hfs.begin(), hfs.end(),
                                  std::bind(&HFHandleCorrection::correctValue, &cor, std::placeholders::_1));
#endif
                    cell(*c_it).set_halffaces(hfs);
                }
            }
        }


        // 3)
        if(has_face_bottom_up_incidences()) {
            assert((size_t)halfface_handle(h, 1).idx() < incident_cell_per_hf_.size());

            incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 1).idx());
            incident_cell_per_hf_.erase(incident_cell_per_hf_.begin() + halfface_handle(h, 0).idx());
        }


        if (!fast_deletion_enabled())
        {
            // 4)
            if(has_edge_bottom_up_incidences()) {
                HFHandleCorrection cor(halfface_handle(h, 1));
#if defined(__clang_major__) && (__clang_major__ >= 5)
                for(std::vector<std::vector<HalfFaceHandle> >::iterator it = incident_hfs_per_he_.begin(), end = incident_hfs_per_he_.end(); it != end; ++it) {
                    cor.correctVecValue(*it);
                }
#else
                std::for_each(incident_hfs_per_he_.begin(),
                              incident_hfs_per_he_.end(),
                              std::bind(&HFHandleCorrection::correctVecValue, &cor, std::placeholders::_1));
#endif
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
CellIter TopologyKernel::delete_cell_core(CellHandle _h) {

    CellHandle h = _h;

    assert(h.is_valid() && (size_t)h.idx() < cells_.size());


    if (fast_deletion_enabled() && !deferred_deletion_enabled()) // for fast deletion swap handle with last not deleted cell
    {
        CellHandle last_undeleted_cell = CellHandle((int)cells_.size()-1);
        assert(!cell_deleted_[last_undeleted_cell]);
        swap_cell_indices(h, last_undeleted_cell);
        h = last_undeleted_cell;
    }


    // 1)
    if(has_face_bottom_up_incidences()) {
        const std::vector<HalfFaceHandle>& hfs = cell(h).halffaces();
        for(std::vector<HalfFaceHandle>::const_iterator hf_it = hfs.begin(),
                hf_end = hfs.end(); hf_it != hf_end; ++hf_it) {
            assert((size_t)hf_it->idx() < incident_cell_per_hf_.size());
            if (incident_cell_per_hf_[*hf_it] == h)
                incident_cell_per_hf_[*hf_it] = InvalidCellHandle;
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
        cell_deleted_[h] = true;
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
            if(has_face_bottom_up_incidences()) {
                CHandleCorrection cor(h);
#if defined(__clang_major__) && (__clang_major__ >= 5)
                for(std::vector<CellHandle>::iterator it = incident_cell_per_hf_.begin(),
                    end = incident_cell_per_hf_.end(); it != end; ++it) {
                    cor.correctValue(*it);
                }
#else
                std::for_each(incident_cell_per_hf_.begin(),
                              incident_cell_per_hf_.end(),
                              std::bind(&CHandleCorrection::correctValue, &cor, std::placeholders::_1));
#endif
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

    // correct pointers to those cells
    for (const auto hfh: cells_[_h1].halffaces()) {
        if (incident_cell_per_hf_[hfh] == _h1)
            incident_cell_per_hf_[hfh] = _h2;
    }
    for (const auto hfh: cells_[_h2].halffaces()) {
        if (incident_cell_per_hf_[hfh] == _h2)
            incident_cell_per_hf_[hfh] = _h1;
    }

    // swap vector entries
    std::swap(cells_[_h1], cells_[_h2]);
    detail::swap_bool(cell_deleted_[_h1], cell_deleted_[_h2]);
    swap_property_elements(_h1, _h2);
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
                CellHandle ch = incident_cell_per_hf_[hfh];
                if (!ch.is_valid())
                    continue;



                if (processed_cells.find(ch.idx()) == processed_cells.end())
                {

                    Cell& c = cells_[ch];

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
            Cell& c = cells_[CellHandle(i)];

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

                    std::vector<HalfFaceHandle>& incident_halffaces = incident_hfs_per_he_[heh];
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
    std::swap(faces_[_h1], faces_[_h2]);
    detail::swap_bool(face_deleted_[_h1], face_deleted_[_h2]);
    std::swap(incident_cell_per_hf_[_h1.halfface_handle(0)], incident_cell_per_hf_[_h2.halfface_handle(0)]);
    std::swap(incident_cell_per_hf_[_h1.halfface_handle(1)], incident_cell_per_hf_[_h2.halfface_handle(1)]);
    swap_property_elements(_h1, _h2);
    swap_property_elements(halfface_handle(_h1, 0), halfface_handle(_h2, 0));
    swap_property_elements(halfface_handle(_h1, 1), halfface_handle(_h2, 1));

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


            std::vector<HalfFaceHandle>& incident_halffaces = incident_hfs_per_he_[heh];
            for (unsigned int j = 0; j < incident_halffaces.size(); ++j) // for each incident halfface
            {
                HalfFaceHandle hfh = incident_halffaces[j];

                unsigned int f_id = hfh.idx() / 2;

                if (processed_faces.find(f_id) == processed_faces.end())
                {

                    Face& f = faces_[FaceHandle(f_id)];

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
            Face& f = faces_[FaceHandle(i)];

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

                std::vector<HalfEdgeHandle>& outgoing_hes = outgoing_hes_per_vertex_[vhs[j]];
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
    std::swap(edges_[_h1], edges_[_h2]);
    detail::swap_bool(edge_deleted_[_h1], edge_deleted_[_h2]);
    std::swap(incident_hfs_per_he_[_h1.halfedge_handle(0)], incident_hfs_per_he_[_h2.halfedge_handle(0)]);
    std::swap(incident_hfs_per_he_[_h1.halfedge_handle(1)], incident_hfs_per_he_[_h2.halfedge_handle(1)]);
    swap_property_elements(_h1, _h2);
    swap_property_elements(halfedge_handle(_h1, 0), halfedge_handle(_h2, 0));
    swap_property_elements(halfedge_handle(_h1, 1), halfedge_handle(_h2, 1));

}

void TopologyKernel::swap_vertex_indices(VertexHandle _h1, VertexHandle _h2)
{
    assert(is_valid(_h1));
    assert(is_valid(_h2));

    if (_h1 == _h2)
        return;

    std::array<VH, 2> ids {_h1, _h2};

    // correct pointers to those vertices

    if (has_vertex_bottom_up_incidences())
    {
        std::set<EH> processed_edges; // to ensure ids are only swapped once (in the case that the two swapped vertices are connected by an edge)
        for (unsigned int i = 0; i < 2; ++i) // For both swapped vertices
        {
            for (const auto heh: outgoing_hes_per_vertex_[ids[i]]) {
                const auto eh = heh.edge_handle();

                if (processed_edges.find(eh) == processed_edges.end())
                {
                    Edge& e = edges_[eh];
                    if (e.from_vertex() == ids[0])
                        e.set_from_vertex(ids[1]);
                    else if (e.from_vertex() == ids[1])
                        e.set_from_vertex(ids[0]);

                    if (e.to_vertex() == ids[0])
                        e.set_to_vertex(ids[1]);
                    else if (e.to_vertex() == ids[1])
                        e.set_to_vertex(ids[0]);

                    processed_edges.insert(eh);
                }
            }
        }

    }
    else
    {
        // search for all edges containing a swapped vertex

        for (auto &e: edges_)
        {
            if (e.from_vertex() == ids[0])
                e.set_from_vertex(ids[1]);
            else if (e.from_vertex() == ids[1])
                e.set_from_vertex(ids[0]);

            if (e.to_vertex() == ids[0])
                e.set_to_vertex(ids[1]);
            else if (e.to_vertex() == ids[1])
                e.set_to_vertex(ids[0]);
        }
    }

    // swap vector entries
    detail::swap_bool(vertex_deleted_[_h1], vertex_deleted_[_h2]);
    std::swap(outgoing_hes_per_vertex_[_h1], outgoing_hes_per_vertex_[_h2]);
    swap_property_elements(_h1, _h2);
}


void TopologyKernel::enable_deferred_deletion(bool _enable)
{
    if (deferred_deletion_ && !_enable)
        collect_garbage();

    deferred_deletion_ = _enable;
}

//========================================================================================

/// Get edge with handle _edgeHandle
const OpenVolumeMeshEdge& TopologyKernel::edge(EdgeHandle _edgeHandle) const
{
    assert(is_valid(_edgeHandle));
    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
const OpenVolumeMeshFace& TopologyKernel::face(FaceHandle _faceHandle) const
{
    assert(is_valid(_faceHandle));
    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
const OpenVolumeMeshCell& TopologyKernel::cell(CellHandle _cellHandle) const
{
    assert(is_valid(_cellHandle));
    return cells_[_cellHandle];
}

//========================================================================================

/// Get edge with handle _edgeHandle
OpenVolumeMeshEdge& TopologyKernel::edge(EdgeHandle _edgeHandle)
{
    assert(is_valid(_edgeHandle));
    return edges_[_edgeHandle];
}

//========================================================================================

/// Get face with handle _faceHandle
OpenVolumeMeshFace& TopologyKernel::face(FaceHandle _faceHandle)
{

    assert(is_valid(_faceHandle));
    return faces_[_faceHandle];
}

//========================================================================================

/// Get cell with handle _cellHandle
OpenVolumeMeshCell& TopologyKernel::cell(CellHandle _cellHandle)
{
    assert(is_valid(_cellHandle));
    return cells_[_cellHandle];
}

//========================================================================================

/// Get edge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::halfedge(HalfEdgeHandle _halfEdgeHandle) const
{
    assert(is_valid(_halfEdgeHandle));

    // In case the handle is even, just return the corresponding edge
    // Otherwise return the opposite halfedge via opposite()
    if(_halfEdgeHandle.subidx() == 0)
        return edges_[_halfEdgeHandle.edge_handle()];
    else
        return opposite_halfedge(edges_[_halfEdgeHandle.edge_handle()]);
}

//========================================================================================

/// Get face that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::halfface(HalfFaceHandle _halfFaceHandle) const {

    assert(is_valid(_halfFaceHandle));

    // In case the handle is even, just return the corresponding face
    // Otherwise return the opposite halfface via opposite()
    if(_halfFaceHandle.idx() % 2 == 0)
        return faces_[_halfFaceHandle.face_handle()];
    else
        return opposite_halfface(faces_[_halfFaceHandle.face_handle()]);
}

//========================================================================================

/// Get opposite halfedge that corresponds to halfedge with handle _halfEdgeHandle
OpenVolumeMeshEdge TopologyKernel::opposite_halfedge(HalfEdgeHandle _halfEdgeHandle) const
{
    return halfedge(_halfEdgeHandle.opposite_handle());
}

//========================================================================================

/// Get opposite halfface that corresponds to halfface with handle _halfFaceHandle
OpenVolumeMeshFace TopologyKernel::opposite_halfface(HalfFaceHandle _halfFaceHandle) const
{
    return halfface(_halfFaceHandle.opposite_handle());
}


//========================================================================================


HalfEdgeHandle TopologyKernel::halfedge(VertexHandle _vh1, VertexHandle _vh2) const
{
  return find_halfedge(_vh1, _vh2);
}


//========================================================================================


HalfEdgeHandle TopologyKernel::find_halfedge(VertexHandle _vh1, VertexHandle _vh2) const
{
    assert(is_valid(_vh1));
    assert(is_valid(_vh2));

    for(VertexOHalfEdgeIter voh_it = voh_iter(_vh1); voh_it.valid(); ++voh_it) {
        if(halfedge(*voh_it).to_vertex() == _vh2) {
            return *voh_it;
        }
    }

    return InvalidHalfEdgeHandle;
}

//========================================================================================

HalfEdgeHandle TopologyKernel::find_halfedge_in_cell(VertexHandle _vh1, VertexHandle _vh2, CellHandle _ch) const
{
  assert(is_valid(_vh1));
  assert(is_valid(_vh2));

  for( auto hfh : cell(_ch).halffaces())
  {
    for(auto heh : halfface(hfh).halfedges() )
    {
      if(from_vertex_handle(heh) == _vh1 && to_vertex_handle(heh) == _vh2)
        return heh;
      if(from_vertex_handle(heh) == _vh2 && to_vertex_handle(heh) == _vh1)
        return opposite_halfedge_handle(heh);
     }
  }

  return InvalidHalfEdgeHandle;
}

//========================================================================================

HalfFaceHandle TopologyKernel::halfface(const std::vector<VertexHandle>& _vs) const
{
  return find_halfface(_vs);
}

//========================================================================================

HalfFaceHandle TopologyKernel::find_halfface(const std::vector<VertexHandle>& _vs) const
{
    assert(_vs.size() > 2);

    VertexHandle v0 = _vs[0], v1 = _vs[1], v2 = _vs[2];

    assert(v0.is_valid() && v1.is_valid() && v2.is_valid());

    HalfEdgeHandle he0 = find_halfedge(v0, v1);
    if(!he0.is_valid()) return InvalidHalfFaceHandle;
    HalfEdgeHandle he1 = find_halfedge(v1, v2);
    if(!he1.is_valid()) return InvalidHalfFaceHandle;

    std::vector<HalfEdgeHandle> hes;
    hes.push_back(he0);
    hes.push_back(he1);

    return find_halfface(hes);
}

//========================================================================================

HalfFaceHandle TopologyKernel::find_halfface_in_cell(const std::vector<VertexHandle>& _vs, CellHandle _ch) const
{
  assert(_vs.size() > 2);

  VertexHandle v0 = _vs[0], v1 = _vs[1], v2 = _vs[2];

  assert(v0.is_valid() && v1.is_valid() && v2.is_valid());

  // check all halfedges of cell until (v0 -> v1) is found and then verify (v0 -> v1 -> v2)
  for( auto hfh : cell(_ch).halffaces())
  {
    for(auto heh : halfface(hfh).halfedges() )
    {
      if(from_vertex_handle(heh) == v0 && to_vertex_handle(heh) == v1)
        if(to_vertex_handle(next_halfedge_in_halfface(heh,hfh)) == v2)
          return hfh;

      // check if opposite halfedge gives desired halfedge
      if(from_vertex_handle(heh) == v1 && to_vertex_handle(heh) == v0)
      {
        HalfEdgeHandle heh_opp = opposite_halfedge_handle(heh);
        HalfFaceHandle hfh_opp = adjacent_halfface_in_cell(hfh,heh);
        if(to_vertex_handle(next_halfedge_in_halfface(heh_opp,hfh_opp)) == v2)
          return hfh_opp;
      }
    }
  }

  return InvalidHalfFaceHandle;
}


//========================================================================================


HalfFaceHandle TopologyKernel::halfface_extensive(const std::vector<VertexHandle>& _vs) const
{
  return find_halfface_extensive(_vs);
}

//========================================================================================


HalfFaceHandle TopologyKernel::find_halfface_extensive(const std::vector<VertexHandle>& _vs) const
{
  //TODO: schöner machen

  assert(_vs.size() > 2);

  VertexHandle v0 = _vs[0];
  VertexHandle v1 = _vs[1];

  assert(v0.is_valid() && v1.is_valid());

  HalfEdgeHandle he0 = find_halfedge(v0, v1);
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

HalfFaceHandle TopologyKernel::halfface(const std::vector<HalfEdgeHandle>& _hes) const
{
  return find_halfface(_hes);
}

//========================================================================================

HalfFaceHandle TopologyKernel::find_halfface(const std::vector<HalfEdgeHandle>& _hes) const
{

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

HalfEdgeHandle TopologyKernel::next_halfedge_in_halfface(HalfEdgeHandle _heh, HalfFaceHandle _hfh) const {

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

HalfEdgeHandle TopologyKernel::prev_halfedge_in_halfface(HalfEdgeHandle _heh, HalfFaceHandle _hfh) const {

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


std::vector<VertexHandle> TopologyKernel::get_halfface_vertices(HalfFaceHandle hfh) const
{
  assert(!halfface(hfh).halfedges().empty());

  return get_halfface_vertices(hfh, halfface(hfh).halfedges().front());
}


//========================================================================================


std::vector<VertexHandle> TopologyKernel::get_halfface_vertices(HalfFaceHandle hfh, VertexHandle vh) const
{
  Face hf = halfface(hfh);
  for (size_t i = 0; i < hf.halfedges().size(); ++i)
    if (halfedge(hf.halfedges()[i]).from_vertex() == vh)
      return get_halfface_vertices(hfh, hf.halfedges()[i]);

  return std::vector<VertexHandle>();
}


//========================================================================================


std::vector<VertexHandle> TopologyKernel::get_halfface_vertices(HalfFaceHandle hfh, HalfEdgeHandle heh) const
{
  std::vector<VertexHandle> vertices;

  // add vertices of halfface
  for (unsigned int i = 0; i < halfface(hfh).halfedges().size(); ++i)
  {
    vertices.push_back(halfedge(heh).from_vertex());
    heh = next_halfedge_in_halfface(heh, hfh);
  }

  return vertices;
}


//========================================================================================


bool
TopologyKernel::
is_incident( FaceHandle _fh, EdgeHandle _eh) const
{
  const auto f = face(_fh);
  for (HalfEdgeHandle heh: f.halfedges())
    if(edge_handle(heh) == _eh)
      return true;

  return false;
}


//========================================================================================


HalfFaceHandle
TopologyKernel::adjacent_halfface_in_cell(HalfFaceHandle _halfFaceHandle,
                                          HalfEdgeHandle _halfEdgeHandle) const
{
  assert(_halfFaceHandle.is_valid() && (size_t)_halfFaceHandle.idx() < faces_.size() * 2u);
  assert(_halfEdgeHandle.is_valid() && (size_t)_halfEdgeHandle.idx() < edges_.size() * 2u);
  assert(has_face_bottom_up_incidences());

  const auto ch = incident_cell(_halfFaceHandle);
  if(!ch.is_valid()) {
    // Specified halfface is on the outside of the complex
    return InvalidHalfFaceHandle;
  }

  // Make sure that _halfFaceHandle is incident to _halfEdgeHandle
  bool skipped = false;
  HalfFaceHandle idx = InvalidHalfFaceHandle;

  // For face-selfadjacent cells, we have to ensure the actual halfedge information
  // is used here, BUT...
  // To support legacy code, we have to flip the halfedge if it's the wrong one
  HalfEdgeHandle hehOpp = opposite_halfedge_handle(_halfEdgeHandle);
  bool hasHalfedge = false;
  bool hasOppHalfedge = false;
  const auto hf = halfface(_halfFaceHandle);
  for (HalfEdgeHandle heh: hf.halfedges()) {
    if (heh == hehOpp)
      hasOppHalfedge = true;
    else if (heh == _halfEdgeHandle)
      hasHalfedge = true;
  }
  if (!hasHalfedge) {
    if (hasOppHalfedge)
      _halfEdgeHandle = hehOpp;
    else
      return InvalidHalfFaceHandle;
  }

  for(const auto &hfh: cell(ch).halffaces()) {
    if(hfh == _halfFaceHandle) {
      assert(!skipped); // a halfface may only appear once in a cell!
      skipped = true;
      if (idx.is_valid()) {
        return idx;
      }
    } else {
      const auto hf_cur = halfface(hfh);
      for (const auto heh: hf_cur.halfedges()) {
        // For face-selfadjacent cells, we look for a halfface that
        // contains the opposite halfedge but isnt the opposite halfface
        if(opposite_halfedge_handle(heh) == _halfEdgeHandle && hfh != opposite_halfface_handle(_halfFaceHandle)) {
          if (idx.is_valid()) {
            // we found two(!) other halffaces that contain the given edge.
            // likely the given halfedge is not part of the given halfface
            return InvalidHalfFaceHandle;
          }
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

CellHandle TopologyKernel::incident_cell(HalfFaceHandle _halfFaceHandle) const
{
    assert((size_t)_halfFaceHandle.idx() < incident_cell_per_hf_.size() && _halfFaceHandle.is_valid());
    assert(has_face_bottom_up_incidences());

    return incident_cell_per_hf_[_halfFaceHandle];
}

//========================================================================================

void TopologyKernel::compute_vertex_bottom_up_incidences() {

    // Clear incidences
    outgoing_hes_per_vertex_.clear();
    outgoing_hes_per_vertex_.resize(n_vertices());
    VertexVector<int> n_edges_per_vertex(n_vertices(), 0);

    for (const auto &eh: edges()) {
        const auto &e = edge(eh);
        ++n_edges_per_vertex[e.from_vertex()];
        ++n_edges_per_vertex[e.to_vertex()];
    }
    for (const auto &vh: vertices()) {
        outgoing_hes_per_vertex_[vh].reserve(n_edges_per_vertex[vh]);
    }

    // Store outgoing halfedges per vertex
    for (const auto eh: edges()) {
        const auto &e = edge(eh);
        VertexHandle from = e.from_vertex();
        // If this condition is not fulfilled, it is out of caller's control and
        // definitely our bug, therefore an assert
        assert((size_t)from.idx() < outgoing_hes_per_vertex_.size());

        outgoing_hes_per_vertex_[from].push_back(halfedge_handle(eh, 0));

        VertexHandle to = e.to_vertex();
        assert((size_t)to.idx() < outgoing_hes_per_vertex_.size());

        // Store opposite halfedge handle
        outgoing_hes_per_vertex_[to].push_back(halfedge_handle(eh, 1));
    }
}

//========================================================================================

void TopologyKernel::compute_edge_bottom_up_incidences() {

    // Clear
    incident_hfs_per_he_.clear();
    incident_hfs_per_he_.resize(n_halfedges());

    EdgeVector<int> n_faces_per_edge(n_edges(), 0);
    for (const auto &fh: faces()) {
        for (const auto &heh: face_halfedges(fh)) {
            ++n_faces_per_edge[edge_handle(heh)];
        }
    }
    for (const auto &eh: edges()) {
        incident_hfs_per_he_[halfedge_handle(eh, 0)].reserve(n_faces_per_edge[eh]);
        incident_hfs_per_he_[halfedge_handle(eh, 1)].reserve(n_faces_per_edge[eh]);
    }
    // Store incident halffaces per halfedge
    for (const auto &fh: faces()) {
        for (const auto &heh: face_halfedges(fh)) {
            auto opp = opposite_halfedge_handle(heh);
            incident_hfs_per_he_[heh].push_back(halfface_handle(fh, 0));
            incident_hfs_per_he_[opp].push_back(halfface_handle(fh, 1));
        }
    }
}

//========================================================================================

void TopologyKernel::compute_face_bottom_up_incidences() {

    // Clear
    incident_cell_per_hf_.clear();
    incident_cell_per_hf_.resize(faces_.size() * 2u, InvalidCellHandle);

    for (const auto ch: cells()) {
        for (const auto hfh: cell_halffaces(ch)) {
            if(incident_cell_per_hf_[hfh] == InvalidCellHandle) {

                incident_cell_per_hf_[hfh] = ch;

            } else {
#ifndef NDEBUG
                std::cerr << "compute_face_bottom_up_incidences(): Detected non-three-manifold configuration!" << std::endl;
                std::cerr << "Connectivity probably won't work." << std::endl;
#endif
            }
        }
    }
}

} // Namespace OpenVolumeMesh
