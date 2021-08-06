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

#include <OpenVolumeMesh/Core/BaseProperty.hh>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/System/Deprecation.hh>
#include <OpenVolumeMesh/Core/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>

namespace OpenVolumeMesh {

class ResourceManager;

/**
 * \class PropertyPtr
 *
 * Provides user access to property contents.
 */

template <class T, typename EntityTag>
class PropertyPtr : public BaseProperty
{
    static_assert(is_entity<EntityTag>::value);
public:
    // defined in ResourceManagerT_impl to avoid circular references
    PropertyPtr(ResourceManager *mesh, std::string _name, T const &_def);

    ~PropertyPtr() {
        auto resman = ptr_->resMan();
        if (resman != nullptr) {
            if (ptr_.use_count() == 1) {
                resman->template property_deleted<EntityTag>(ptr_.get());
            }
        }
    }

    operator bool() const {return ptr_->resMan() != nullptr;}

    friend class ResourceManager;
    friend class PropertyStorageT<T>;

    using PropStorageT = PropertyStorageT<T>;

    typedef typename PropStorageT::value_type                  value_type;
    typedef typename PropStorageT::vector_type::const_iterator const_iterator;
    typedef typename PropStorageT::vector_type::iterator       iterator;
    typedef typename PropStorageT::reference                   reference;
    typedef typename PropStorageT::const_reference             const_reference;

    using EntityHandleT = HandleT<EntityTag>;

    const std::string& name() const & override {
        // the string we return a reference to lives long enough, no warnings please:
        // cppcheck-suppress returnTempReference
        return ptr_->name();
    }

    const_iterator begin() const { return ptr_->begin(); }
    const_iterator end() const   { return ptr_->end(); }
    iterator begin() { return ptr_->begin(); }
    iterator end()   { return ptr_->end(); }
    size_t size() const { return ptr_->size(); }

    /// No range check performed!
    reference operator[](const EntityHandleT& _h) { return (*ptr_)[_h.uidx()]; }

    /// No range check performed!
    const_reference operator[](const EntityHandleT& _h) const { return (*ptr_)[_h.uidx()]; }

    reference       at(const EntityHandleT& _h) { return ptr_->at(_h.uidx()); }
    const_reference at(const EntityHandleT& _h) const { return ptr_->at(_h.uidx()); }

    void serialize(std::ostream& _ostr) const override { ptr_->serialize(_ostr); }
    void deserialize(std::istream& _istr) override { ptr_->deserialize(_istr); }

    bool persistent() const override { return ptr_->persistent(); }
    bool anonymous() const override { return ptr_->anonymous(); }
    size_t n_elements() const { return ptr_->n_elements(); }

    std::string typeNameWrapper() const override {return ptr_->typeNameWrapper(); }
    EntityType entity_type() const {return ptr_->entity_type();}

    [[deprecated]]
    void copy(const EntityHandleT &_src, const EntityHandleT &_dst) {
        (*this)[_dst] = (*this)[_src];
    }

    /// get default value.
    T const &def() const {return ptr_->def();}

    /// set all values to `val`.
    void fill(T const&val) { ptr_->fill(val); }


protected:
     PropertyPtr(std::shared_ptr<PropStorageT> &&_ptr)
         : ptr_(std::move(_ptr))
     {}
     std::shared_ptr<PropStorageT> const &ptr() const {return ptr_;}
     std::shared_ptr<PropStorageT> ptr_;
};

template<typename T>
PropertyStorageT<T>::operator std::unique_ptr<BaseProperty>()
{
    auto sp = std::static_pointer_cast<PropertyStorageT<T>>(shared_from_this());
    return entitytag_dispatch(entity_type(), [](auto entity_tag, auto &&sp)
    {
    // unfortunately we cannot use make_unique due to our protected constructor
            auto ptr = new PropertyPtr<T, decltype(entity_tag)>(std::move(sp));
            return std::unique_ptr<BaseProperty>(ptr);
    }, std::move(sp));
}

template<typename T> using VertexPropertyT   = PropertyPtr<T, Entity::Vertex>;
template<typename T> using EdgePropertyT     = PropertyPtr<T, Entity::Edge>;
template<typename T> using HalfEdgePropertyT = PropertyPtr<T, Entity::HalfEdge>;
template<typename T> using FacePropertyT     = PropertyPtr<T, Entity::Face>;
template<typename T> using HalfFacePropertyT = PropertyPtr<T, Entity::HalfFace>;
template<typename T> using CellPropertyT     = PropertyPtr<T, Entity::Cell>;
template<typename T> using MeshPropertyT     = PropertyPtr<T, Entity::Mesh>;

} // Namespace OpenVolumeMesh
