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


#include <OpenVolumeMesh/Core/Iterators.hh>
#include <OpenVolumeMesh/Core/Iterators/BoundaryHalfFaceHalfFaceIter.cc>
#include <OpenVolumeMesh/Core/Iterators/BoundaryItemIter.cc>
#include <OpenVolumeMesh/Core/Iterators/CellCellIter.cc>
#include <OpenVolumeMesh/Core/Iterators/CellIter.cc>
#include <OpenVolumeMesh/Core/Iterators/CellVertexIter.cc>
#include <OpenVolumeMesh/Core/Iterators/EdgeIter.cc>
#include <OpenVolumeMesh/Core/Iterators/FaceIter.cc>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeCellIter.cc>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeHalfFaceIter.cc>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeIter.cc>
#include <OpenVolumeMesh/Core/Iterators/HalfFaceIter.cc>
#include <OpenVolumeMesh/Core/Iterators/HalfFaceVertexIter.cc>
#include <OpenVolumeMesh/Core/Iterators/VertexCellIter.cc>
#include <OpenVolumeMesh/Core/Iterators/VertexFaceIter.cc>
#include <OpenVolumeMesh/Core/Iterators/VertexIter.cc>
#include <OpenVolumeMesh/Core/Iterators/VertexOHalfEdgeIter.cc>
#include <OpenVolumeMesh/Core/Iterators/VertexVertexIter.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/CellEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/CellFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeCellIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeHalfFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceHalfEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceVertexIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfEdgeFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceHalfEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexEdgeIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexHalfFaceIterImpl.cc>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexIHalfEdgeIterImpl.cc>

namespace OpenVolumeMesh {
template class OVM_EXPORT GenericCirculator<detail::VertexIHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::VertexEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::VertexHalfFaceIterImpl>;

template class OVM_EXPORT GenericCirculator<detail::HalfEdgeFaceIterImpl>;

template class OVM_EXPORT GenericCirculator<detail::EdgeHalfFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::EdgeFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::EdgeCellIterImpl>;

template class OVM_EXPORT GenericCirculator<detail::HalfFaceHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::HalfFaceEdgeIterImpl>;

template class OVM_EXPORT GenericCirculator<detail::FaceVertexIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::FaceHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::FaceEdgeIterImpl>;

template class OVM_EXPORT GenericCirculator<detail::CellHalfEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::CellEdgeIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::CellHalfFaceIterImpl>;
template class OVM_EXPORT GenericCirculator<detail::CellFaceIterImpl>;


} // namespace OpenVolumeMesh
