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


#include <type_traits>
#include <iterator>
#include <set>
#include <vector>

#ifndef NDEBUG
#include <iostream>
#endif

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>
#include <OpenVolumeMesh/Core/Iterators/GenericCirculator.hh>

#include <OpenVolumeMesh/Core/Iterators/detail/CellEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/CellFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeCellIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/EdgeHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/FaceVertexIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfEdgeFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/Iterators/detail/VertexIHalfEdgeIterImpl.hh>


#include <OpenVolumeMesh/Core/Iterators/BoundaryHalfFaceHalfFaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/BoundaryItemIter.hh>
#include <OpenVolumeMesh/Core/Iterators/CellCellIter.hh>
#include <OpenVolumeMesh/Core/Iterators/CellIter.hh>
#include <OpenVolumeMesh/Core/Iterators/CellVertexIter.hh>
#include <OpenVolumeMesh/Core/Iterators/EdgeIter.hh>
#include <OpenVolumeMesh/Core/Iterators/FaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeCellIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeHalfFaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfFaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfFaceVertexIter.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexCellIter.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexFaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexIter.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexOHalfEdgeIter.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexVertexIter.hh>



namespace OpenVolumeMesh {


extern template class OVM_EXPORT GenericCirculator<detail::VertexIHalfEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::VertexEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::VertexHalfFaceIterImpl>;

extern template class OVM_EXPORT GenericCirculator<detail::HalfEdgeFaceIterImpl>;

extern template class OVM_EXPORT GenericCirculator<detail::EdgeHalfFaceIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::EdgeFaceIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::EdgeCellIterImpl>;

extern template class OVM_EXPORT GenericCirculator<detail::HalfFaceHalfEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::HalfFaceEdgeIterImpl>;

extern template class OVM_EXPORT GenericCirculator<detail::FaceVertexIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::FaceHalfEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::FaceEdgeIterImpl>;

extern template class OVM_EXPORT GenericCirculator<detail::CellHalfEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::CellEdgeIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::CellHalfFaceIterImpl>;
extern template class OVM_EXPORT GenericCirculator<detail::CellFaceIterImpl>;





typedef GenericCirculator<detail::VertexIHalfEdgeIterImpl> VertexIHalfEdgeIter;
typedef GenericCirculator<detail::VertexEdgeIterImpl> VertexEdgeIter;
typedef GenericCirculator<detail::VertexHalfFaceIterImpl> VertexHalfFaceIter;

typedef GenericCirculator<detail::HalfEdgeFaceIterImpl> HalfEdgeFaceIter;

typedef GenericCirculator<detail::EdgeHalfFaceIterImpl> EdgeHalfFaceIter;
typedef GenericCirculator<detail::EdgeFaceIterImpl> EdgeFaceIter;
typedef GenericCirculator<detail::EdgeCellIterImpl> EdgeCellIter;

typedef GenericCirculator<detail::HalfFaceHalfEdgeIterImpl> HalfFaceHalfEdgeIter;
typedef GenericCirculator<detail::HalfFaceEdgeIterImpl> HalfFaceEdgeIter;

typedef GenericCirculator<detail::FaceVertexIterImpl> FaceVertexIter;
typedef GenericCirculator<detail::FaceHalfEdgeIterImpl> FaceHalfEdgeIter;
typedef GenericCirculator<detail::FaceEdgeIterImpl> FaceEdgeIter;

typedef GenericCirculator<detail::CellHalfEdgeIterImpl> CellHalfEdgeIter;
typedef GenericCirculator<detail::CellEdgeIterImpl> CellEdgeIter;
typedef GenericCirculator<detail::CellHalfFaceIterImpl> CellHalfFaceIter;
typedef GenericCirculator<detail::CellFaceIterImpl> CellFaceIter;

//===========================================================================

//===========================================================================

} // Namespace OpenVolumeMesh
