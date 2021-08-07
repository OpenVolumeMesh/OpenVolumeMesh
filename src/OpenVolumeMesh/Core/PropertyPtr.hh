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

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/System/Deprecation.hh>
#include <OpenVolumeMesh/Core/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>

namespace OpenVolumeMesh {

class ResourceManager;

/**
 * \class PropertyPtr
 *
 * Provides handle-type-safe user access to property contents.
 */

template <typename T, typename EntityTag>
class PropertyPtr : public PropertyStoragePtr<T>
{
    static_assert(is_entity<EntityTag>::value);
    using PropertyStoragePtr<T>::storage;
public:
    using PropertyStoragePtr<T>::begin;
    using PropertyStoragePtr<T>::end;
    using PropertyStoragePtr<T>::size;
    using PropertyStoragePtr<T>::serialize;
    using PropertyStoragePtr<T>::deserialize;
    using PropertyStoragePtr<T>::persistent;
    using PropertyStoragePtr<T>::anonymous;
    using PropertyStoragePtr<T>::n_elements;
    using PropertyStoragePtr<T>::typeNameWrapper;
    using PropertyStoragePtr<T>::entity_type;
    using PropertyStoragePtr<T>::operator bool;
    using PropertyStoragePtr<T>::def;
    using PropertyStoragePtr<T>::fill;
    using PropertyStoragePtr<T>::name;
    using reference = typename PropertyStoragePtr<T>::reference;
    using const_reference = typename PropertyStoragePtr<T>::const_reference;

    // defined in ResourceManagerT_impl to avoid circular references
    PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def);

    ~PropertyPtr() {
        // TODO: move this to extra parent class, then implement default assignment ops and copy constructors
        auto resman = storage()->resMan();
        if (resman != nullptr) {
            if (storage().use_count() == 1) {
                resman->template property_deleted<EntityTag>(storage().get());
            }
        }
    }

    friend class ResourceManager;
    friend class PropertyStorageT<T>;

    using EntityHandleT = HandleT<EntityTag>;

    /// No range check performed!
    reference operator[](const EntityHandleT& _h) { return (*storage())[_h.uidx()]; }

    /// No range check performed!
    const_reference operator[](const EntityHandleT& _h) const { return (*storage())[_h.uidx()]; }

    reference       at(const EntityHandleT& _h) { return storage()->at(_h.uidx()); }
    const_reference at(const EntityHandleT& _h) const { return storage()->at(_h.uidx()); }


    [[deprecated]]
    void copy(const EntityHandleT &_src, const EntityHandleT &_dst) {
        (*this)[_dst] = (*this)[_src];
    }


protected:
     PropertyPtr(std::shared_ptr<PropertyStorageT<T>> &&_ptr)
         : PropertyStoragePtr<T>(std::move(_ptr))
     {}
};

template<typename T> using VertexPropertyT   = PropertyPtr<T, Entity::Vertex>;
template<typename T> using EdgePropertyT     = PropertyPtr<T, Entity::Edge>;
template<typename T> using HalfEdgePropertyT = PropertyPtr<T, Entity::HalfEdge>;
template<typename T> using FacePropertyT     = PropertyPtr<T, Entity::Face>;
template<typename T> using HalfFacePropertyT = PropertyPtr<T, Entity::HalfFace>;
template<typename T> using CellPropertyT     = PropertyPtr<T, Entity::Cell>;
template<typename T> using MeshPropertyT     = PropertyPtr<T, Entity::Mesh>;

} // Namespace OpenVolumeMesh
