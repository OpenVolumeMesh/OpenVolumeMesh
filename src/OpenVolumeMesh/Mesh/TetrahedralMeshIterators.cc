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

#include <set>

#include <OpenVolumeMesh/Mesh/TetrahedralMeshIterators.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Core/Iterators.hh>

namespace OpenVolumeMesh {

//================================================================================================
// TetVertexIter
//================================================================================================


TetVertexIter::TetVertexIter(CellHandle _ref_h,
        const TetrahedralMeshTopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    assert(_ref_h.is_valid());

    assert(_mesh->valence(_ref_h) == 4);

    const auto& cell_vhs = _mesh->get_cell_vertices(_ref_h);
    vertices_[0] = cell_vhs[0];
    vertices_[1] = cell_vhs[1];
    vertices_[2] = cell_vhs[2];
    vertices_[3] = cell_vhs[3];

    cur_index_ = 0;
    BaseIter::cur_handle(vertices_[cur_index_]);
}


TetVertexIter& TetVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = vertices_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(vertices_[cur_index_]);

    return *this;
}


TetVertexIter& TetVertexIter::operator++() {

    ++cur_index_;
    if(cur_index_ == vertices_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(vertices_[cur_index_]);

    return *this;
}

} // Namespace OpenVolumeMesh
