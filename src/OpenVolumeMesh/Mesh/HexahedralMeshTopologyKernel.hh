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
#include <OpenVolumeMesh/Mesh/HexahedralMeshIterators.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh {

/**
 * \class HexahedralMeshTopologyKernel
 *
 * \brief A data structure basing on PolyhedralMesh with specializations for hexahedra.
 *
 * The hexahedron has an induced "virtual" coordinate system. This supposes
 * the incident half-faces to be given in a specific order.
 * See the following figure for an illustration of the induced
 * coordinate system.
 *
 * \image html induced_coordsys.png
 *
 * The abbreviations XF, XB, etc. are short for
 *
 * \li \c XF: X-axis front face
 * \li \c XB: X-axis back face
 * \li \c YF: Y-axis front face
 * \li \c ...
 *
 * The axes refer to the intrinsic "virtual" axes of the hexahedron.
 * The incident half-faces have to be defined in the following order:
 *
 * \li \c 1. XF
 * \li \c 2. XB
 * \li \c 3. YF
 * \li \c 4. YB
 * \li \c 5. ZF
 * \li \c 6. ZB
 */

class OVM_EXPORT HexahedralMeshTopologyKernel : public TopologyKernel {
public:

    // Orientation constants
    static const unsigned char XF = 0;
    static const unsigned char XB = 1;
    static const unsigned char YF = 2;
    static const unsigned char YB = 3;
    static const unsigned char ZF = 4;
    static const unsigned char ZB = 5;
    static const unsigned char INVALID = 6;

    static inline unsigned char opposite_orientation(const unsigned char _d) {
        return static_cast<unsigned char>(_d % 2 == 0 ? _d + 1 : _d - 1);
    }

    // Constructor
    HexahedralMeshTopologyKernel() = default;

    // Destructor
    ~HexahedralMeshTopologyKernel() override = default;

    // Overridden function
    FaceHandle add_face(std::vector<HalfEdgeHandle> _halfedges, bool _topologyCheck = false) override;

    // Overridden function
    FaceHandle add_face(const std::vector<VertexHandle>& _vertices) override;

    /// Overridden function
    CellHandle add_cell(std::vector<HalfFaceHandle> _halffaces, bool _topologyCheck = false) override;

private:

    bool check_halfface_ordering(const std::vector<HalfFaceHandle>& _hfs) const;

public:

    /** \brief Add cell via incident vertices
     *
     * Test whether all required faces are already defined
     * and, if not, create them.
     * Give vertices in the following order:
     *
     *      5-------6
     *     /|      /|
     *    / |     / |
     *   3-------2  |
     *   |  4----|--7
     *   | /     | /
     *   |/      |/
     *   0-------1
     *
     * @param _vertices A list of vertices in the correct order
     * @param _topologyCheck Flag if a topology check should be performed before adding the cell
     *
     * @return The new hexahedron's cell handle
     */
    CellHandle add_cell(const std::vector<VertexHandle>& _vertices, bool _topologyCheck = false);

    // ======================= Specialized Iterators =============================

    friend class CellSheetCellIter;
    friend class HalfFaceSheetHalfFaceIter;
    friend class HexVertexIter;

    typedef class CellSheetCellIter CellSheetCellIter;
    typedef class HalfFaceSheetHalfFaceIter HalfFaceSheetHalfFaceIter;
    typedef class HexVertexIter HexVertexIter;

    using TopologyKernel::is_valid;

    CellSheetCellIter csc_iter(CellHandle _ref_h, const unsigned char _orthDir, int _max_laps = 1) const {
        return CellSheetCellIter(_ref_h, _orthDir, this, _max_laps);
    }

    std::pair<CellSheetCellIter,CellSheetCellIter> cell_sheet_cells(CellHandle _ref_h, const unsigned char _orthDir, int _max_laps = 1) const {
        CellSheetCellIter begin = csc_iter(_ref_h, _orthDir, _max_laps);
        CellSheetCellIter end   = make_end_circulator(begin);
        return std::make_pair(begin, end);
    }

    HalfFaceSheetHalfFaceIter hfshf_iter(HalfFaceHandle _ref_h, int _max_laps = 1) const {
        return HalfFaceSheetHalfFaceIter(_ref_h, this, _max_laps);
    }

    std::pair<HalfFaceSheetHalfFaceIter,HalfFaceSheetHalfFaceIter> halfface_sheet_halffaces(HalfFaceHandle _ref_h, int _max_laps = 1) const {
        HalfFaceSheetHalfFaceIter begin = hfshf_iter(_ref_h, _max_laps);
        HalfFaceSheetHalfFaceIter end   = make_end_circulator(begin);
        return std::make_pair(begin, end);
    }

    HexVertexIter hv_iter(CellHandle _ref_h, int _max_laps = 1) const {
        return HexVertexIter(_ref_h, this, _max_laps);
    }

    std::pair<HexVertexIter,HexVertexIter> hex_vertices(CellHandle _ref_h, int _max_laps = 1) const {
        HexVertexIter begin = hv_iter(_ref_h, _max_laps);
        HexVertexIter end   = make_end_circulator(begin);
        return std::make_pair(begin, end);
    }

    // ======================= Connectivity functions =============================

    inline HalfFaceHandle opposite_halfface_handle_in_cell(HalfFaceHandle _hfh, CellHandle _ch) const {

        assert(is_valid(_hfh));
        assert(is_valid(_ch));

        if(orientation(_hfh, _ch) == XF) return xback_halfface(_ch);
        if(orientation(_hfh, _ch) == XB) return xfront_halfface(_ch);
        if(orientation(_hfh, _ch) == YF) return yback_halfface(_ch);
        if(orientation(_hfh, _ch) == YB) return yfront_halfface(_ch);
        if(orientation(_hfh, _ch) == ZF) return zback_halfface(_ch);
        if(orientation(_hfh, _ch) == ZB) return zfront_halfface(_ch);

        return TopologyKernel::InvalidHalfFaceHandle;
    }

    inline HalfFaceHandle xfront_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[XF];
    }

    inline HalfFaceHandle xback_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[XB];
    }

    inline HalfFaceHandle yfront_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[YF];
    }

    inline HalfFaceHandle yback_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[YB];
    }

    inline HalfFaceHandle zfront_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[ZF];
    }

    inline HalfFaceHandle zback_halfface(CellHandle _ch) const {
        assert(is_valid(_ch));

        return TopologyKernel::cell(_ch).halffaces()[ZB];
    }

    unsigned char orientation(HalfFaceHandle _hfh, CellHandle _ch) const {
        assert(is_valid(_ch));

        std::vector<HalfFaceHandle> halffaces = TopologyKernel::cell(_ch).halffaces();
        for(unsigned int i = 0; i < halffaces.size(); ++i) {
            if(halffaces[i] == _hfh) return (unsigned char)i;
        }

        return INVALID;
    }

    static inline unsigned char orthogonal_orientation(const unsigned char _o1, const unsigned char _o2) {

        if(_o1 == XF && _o2 == YF) return ZF;
        if(_o1 == XF && _o2 == YB) return ZB;
        if(_o1 == XF && _o2 == ZF) return YB;
        if(_o1 == XF && _o2 == ZB) return YF;
        if(_o1 == XB && _o2 == YF) return ZB;
        if(_o1 == XB && _o2 == YB) return ZF;
        if(_o1 == XB && _o2 == ZF) return YF;
        if(_o1 == XB && _o2 == ZB) return YB;

        if(_o1 == YF && _o2 == XF) return ZB;
        if(_o1 == YF && _o2 == XB) return ZF;
        if(_o1 == YF && _o2 == ZF) return XF;
        if(_o1 == YF && _o2 == ZB) return XB;
        if(_o1 == YB && _o2 == XF) return ZF;
        if(_o1 == YB && _o2 == XB) return ZB;
        if(_o1 == YB && _o2 == ZF) return XB;
        if(_o1 == YB && _o2 == ZB) return XF;

        if(_o1 == ZF && _o2 == YF) return XB;
        if(_o1 == ZF && _o2 == YB) return XF;
        if(_o1 == ZF && _o2 == XF) return YF;
        if(_o1 == ZF && _o2 == XB) return YB;
        if(_o1 == ZB && _o2 == YF) return XF;
        if(_o1 == ZB && _o2 == YB) return XB;
        if(_o1 == ZB && _o2 == XF) return YB;
        if(_o1 == ZB && _o2 == XB) return YF;

        return INVALID;

    }

    inline HalfFaceHandle get_oriented_halfface(const unsigned char _o, CellHandle _ch) const {

        if(_o == XF) return xfront_halfface(_ch);
        if(_o == XB) return xback_halfface(_ch);
        if(_o == YF) return yfront_halfface(_ch);
        if(_o == YB) return yback_halfface(_ch);
        if(_o == ZF) return zfront_halfface(_ch);
        if(_o == ZB) return zback_halfface(_ch);
        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_sheet(HalfFaceHandle _hfh, HalfEdgeHandle _heh) const {

        if(!TopologyKernel::has_face_bottom_up_incidences()) {
#ifndef NDEBUG
            std::cerr << "No bottom-up incidences computed so far, could not get adjacent halfface on sheet!" << std::endl;
#endif
            return TopologyKernel::InvalidHalfFaceHandle;
        }

        HalfFaceHandle n_hf = _hfh;
        HalfEdgeHandle n_he = _heh;

        // Try the 1st way
        while(true) {
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::opposite_halfface_handle(n_hf);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            else return n_hf;
        }

        n_hf = TopologyKernel::opposite_halfface_handle(_hfh);
        n_he = TopologyKernel::opposite_halfedge_handle(_heh);

        // Try the 2nd way
        while(true) {
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::opposite_halfface_handle(n_hf);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            n_hf = TopologyKernel::adjacent_halfface_in_cell(n_hf, n_he);
            if(n_hf == TopologyKernel::InvalidHalfFaceHandle) break;
            else return TopologyKernel::opposite_halfface_handle(n_hf);
        }

        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle adjacent_halfface_on_surface(HalfFaceHandle _hfh, HalfEdgeHandle _heh) const {

        for(OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = TopologyKernel::hehf_iter(_heh);
                hehf_it.valid(); ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(TopologyKernel::is_boundary(*hehf_it)) {
                return *hehf_it;
            }
            if(TopologyKernel::is_boundary(TopologyKernel::opposite_halfface_handle(*hehf_it))) {
                return TopologyKernel::opposite_halfface_handle(*hehf_it);
            }
        }
        return TopologyKernel::InvalidHalfFaceHandle;
    }

    HalfFaceHandle neighboring_outside_halfface(HalfFaceHandle _hfh, HalfEdgeHandle _heh) const {

        if(!TopologyKernel::has_face_bottom_up_incidences()) {
#ifndef NDEBUG
            std::cerr << "No bottom-up incidences computed so far, could not get neighboring outside halfface!" << std::endl;
#endif
            return TopologyKernel::InvalidHalfFaceHandle;
        }

        for(OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = TopologyKernel::hehf_iter(_heh);
                hehf_it; ++hehf_it) {
            if(*hehf_it == _hfh) continue;
            if(TopologyKernel::is_boundary(*hehf_it)) return *hehf_it;
            if(TopologyKernel::is_boundary(TopologyKernel::opposite_halfface_handle(*hehf_it)))
                return TopologyKernel::opposite_halfface_handle(*hehf_it);
        }

        return TopologyKernel::InvalidHalfFaceHandle;
    }

private:

    HalfFaceHandle get_adjacent_halfface(HalfFaceHandle _hfh, HalfEdgeHandle _heh,
            const std::vector<HalfFaceHandle>& _halffaces) const;

};

} // Namespace OpenVolumeMesh
