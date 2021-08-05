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


#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Core/Entities.hh>

namespace OpenVolumeMesh {

#if 0
template <class T>
const std::string typeName();
#endif

template<typename T>
const std::string typeName() {
    throw std::runtime_error("Serialization is not supported for this data type!");
}


template <> OVM_EXPORT const std::string typeName<int>();
template <> OVM_EXPORT const std::string typeName<unsigned int>();
template <> OVM_EXPORT const std::string typeName<short>();
template <> OVM_EXPORT const std::string typeName<long>();
template <> OVM_EXPORT const std::string typeName<unsigned long>();
template <> OVM_EXPORT const std::string typeName<char>();
template <> OVM_EXPORT const std::string typeName<unsigned char>();
template <> OVM_EXPORT const std::string typeName<bool>();
template <> OVM_EXPORT const std::string typeName<float>();
template <> OVM_EXPORT const std::string typeName<double>();
template <> OVM_EXPORT const std::string typeName<std::string>();
template <> OVM_EXPORT const std::string typeName<std::map<HalfEdgeHandle, int> >();
template <> OVM_EXPORT const std::string typeName<std::vector<double> >();
template <> OVM_EXPORT const std::string typeName<std::vector<VertexHandle> >();
template <> OVM_EXPORT const std::string typeName<std::vector<HalfFaceHandle> >();
template <> OVM_EXPORT const std::string typeName<std::vector<std::vector<HalfFaceHandle> > >();

template<typename Entity>
const std::string entityTypeName();

template <> OVM_EXPORT const std::string entityTypeName<Entity::Vertex>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::HalfEdge>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::Edge>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::Face>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::HalfFace>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::Cell>();
template <> OVM_EXPORT const std::string entityTypeName<Entity::Mesh>();

std::string entityTypeName(EntityType _type);

} // Namespace OpenVolumeMesh
