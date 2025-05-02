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
#include <memory>

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/System/Deprecation.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>
#include <OpenVolumeMesh/Core/HandleIndexing.hh>

namespace OpenVolumeMesh {

class ResourceManager;

class OVM_EXPORT BasePropertyPtr {
public:
    virtual ~BasePropertyPtr() = default;
    virtual std::string const& name() const& = 0;
};

/**
 * \class PropertyPtr
 *
 * Provides handle-type-safe user access to property contents.
 */

template <typename T, typename EntityTag>
class PropertyPtr : public HandleIndexing<EntityTag, PropertyStoragePtr<T>>
                  , public BasePropertyPtr
{
    static_assert(is_entity<EntityTag>::value);
    using PropertyStoragePtr<T>::storage;
public:
    using EntityTagT = EntityTag;
    // defined in ResourceManagerT_impl to avoid circular references
    PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def=T());

    std::string const& name() const& override {return PropertyStoragePtr<T>::name();};

    friend class ResourceManager;
    template<typename _T>
    friend class PropertyStorageT;

protected:
    using HandleIndexing<EntityTag, PropertyStoragePtr<T>>::HandleIndexing;
};

template<typename T>
std::unique_ptr<BasePropertyPtr>
PropertyStorageT<T>::make_property_ptr()
{
    return entitytag_dispatch(entity_type(), [&](auto entity){
        using EntityTag = decltype(entity);
        auto sp = std::static_pointer_cast<PropertyStorageT<T>>(shared_from_this());
        auto pp = new PropertyPtr<T, EntityTag>(std::move(sp));
        return std::unique_ptr<BasePropertyPtr>(pp);
    });
}

// aliases for legacy code:
template<typename T> using VertexPropertyT   = PropertyPtr<T, Entity::Vertex>;
template<typename T> using EdgePropertyT     = PropertyPtr<T, Entity::Edge>;
template<typename T> using HalfEdgePropertyT = PropertyPtr<T, Entity::HalfEdge>;
template<typename T> using FacePropertyT     = PropertyPtr<T, Entity::Face>;
template<typename T> using HalfFacePropertyT = PropertyPtr<T, Entity::HalfFace>;
template<typename T> using CellPropertyT     = PropertyPtr<T, Entity::Cell>;
template<typename T> using MeshPropertyT     = PropertyPtr<T, Entity::Mesh>;

// aliases with more consistent naming:
template<typename T> using VertexPropertyPtr   = PropertyPtr<T, Entity::Vertex>;
template<typename T> using EdgePropertyPtr     = PropertyPtr<T, Entity::Edge>;
template<typename T> using HalfEdgePropertyPtr = PropertyPtr<T, Entity::HalfEdge>;
template<typename T> using FacePropertyPtr     = PropertyPtr<T, Entity::Face>;
template<typename T> using HalfFacePropertyPtr = PropertyPtr<T, Entity::HalfFace>;
template<typename T> using CellPropertyPtr     = PropertyPtr<T, Entity::Cell>;
template<typename T> using MeshPropertyPtr     = PropertyPtr<T, Entity::Mesh>;

} // Namespace OpenVolumeMesh
