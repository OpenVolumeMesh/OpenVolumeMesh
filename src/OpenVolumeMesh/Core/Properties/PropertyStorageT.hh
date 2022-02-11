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


#include <cassert>
#include <istream>
#include <ostream>
#include <numeric>
#include <string>
#include <vector>
#include <memory>

#include <OpenVolumeMesh/FileManager/TypeNames.hh>
#include <OpenVolumeMesh/FileManager/Serializers.hh>

#include <OpenVolumeMesh/Core/Properties/PropertyStorageBase.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyStoragePtr.hh>

#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>

namespace OpenVolumeMesh {

/** \class PropertyStorageT
 *
 *  \brief Data storage for properties of type T.
 *
 */

template<class T>
class PropertyStorageT: public PropertyStorageBase {
public:

    template <class PropT, class Entity> friend class PropertyPtr;
    template <class PropT> friend class PropertyStoragePtr;
    friend class ResourceManager;

    typedef T                                         Value;
    typedef typename std::vector<T>                   vector_type;
    typedef T                                         value_type;
    typedef typename vector_type::reference           reference;
    typedef typename vector_type::const_reference     const_reference;

public:

    explicit PropertyStorageT(
            detail::Tracker<PropertyStorageBase> *tracker,
            std::string _name,
            EntityType _entity_type,
            const T _def,
            bool _shared)
        : PropertyStorageBase(tracker, std::move(_name), detail::internal_type_name<T>(), _entity_type, _shared),
          def_(std::move(_def))
    {}

public:
    void reserve(size_t _n) final {
        data_.reserve(_n);
    }
    void resize(size_t _n) final  {
        data_.resize(_n, def_);
    }
    size_t size() const final  {
        return data_.size();
    }
    void clear() final  {
        data_.clear();
    }
    void push_back() final  {
        data_.push_back(def_);
    }
    void swap(size_t _i0, size_t _i1) final  {
        std::swap(data_[_i0], data_[_i1]);
    }
    void copy(size_t _src_idx, size_t _dst_idx) final{
        data_[_dst_idx] = data_[_src_idx];
    }
    void delete_element(size_t _idx) final {
        assert(_idx < data_.size());
        data_.erase(data_.begin() + static_cast<long>(_idx));
    }
    void swap(std::vector<T> &_other) {
        if (data_.size() != _other.size()) {
            throw std::runtime_error("PropertyStorageT::swap: vector sizes don't match");
        }
        std::swap(data_, _other);
    }

    /// This is made for use in OpenFlipper's Plugin-PropertyVis,
    /// we recommend avoiding it. No guarantees on preserving this API.
    std::unique_ptr<BasePropertyPtr> make_property_ptr();

public:

    std::string typeNameWrapper() const final {return OpenVolumeMesh::typeName<T>();}

    // Function to serialize a property
    void serialize(std::ostream& _ostr) const final {
        for(typename vector_type::const_iterator it = data_.begin();
                it != data_.end(); ++it) {
            OpenVolumeMesh::serialize(_ostr, *it) << std::endl;
        }
    }

    // Function to deserialize a property
    void deserialize(std::istream& _istr) final {
        for(unsigned int i = 0; i < size(); ++i) {
            OpenVolumeMesh::deserialize(_istr, data_[i]);
        }
    }

public:
    // data access interface

    /// Get pointer to array (does not work for T==bool)
    const T* data() const {

        if (data_.empty())
            return 0;

        return data_.data();
    }

    /// Get reference to property vector (be careful, improper usage, e.g. resizing, may crash)
    vector_type& data_vector() {
        return data_;
    }

    const vector_type& data_vector() const {
        return data_;
    }

    /// Access the i'th element. No range check is performed!
    reference operator[](size_t _idx) {
        assert(_idx < data_.size());
        return data_[_idx];
    }

    /// Const access to the i'th element. No range check is performed!
    const_reference operator[](size_t _idx) const {
        assert(_idx < data_.size());
        return data_[_idx];
    }

    const_reference at(size_t _idx) const { return data_.at(_idx); }
    reference       at(size_t _idx)       { return data_.at(_idx); }

    std::shared_ptr<PropertyStorageBase> clone() const final {
        auto res = std::make_shared<PropertyStorageT<T>>(*this);
        res->set_tracker(nullptr);
        return res;
    }

    typename vector_type::const_iterator begin() const { return data_.begin(); }

    typename vector_type::iterator begin() { return data_.begin(); }

    typename vector_type::const_iterator end() const { return data_.end(); }

    typename vector_type::iterator end() { return data_.end(); }


    T const& def() const {return def_;}

    void fill(T const&val) {
        std::fill(data_.begin(), data_.end(), val);
    }


    void assign_values_from(const PropertyStorageBase *_other) final
    {
        const auto *other = _other->cast_to_StorageT<T>();
        data_ = other->data_;
        def_ = other->def_;

    }

    void move_values_from(PropertyStorageBase *_other) final
    {
        const auto *other = _other->cast_to_StorageT<T>();
        data_ = std::move(other->data_);
        def_ = std::move(other->def_);
    }

private:
    vector_type data_;
    T def_;
};



//-----------------------------------------------------------------------------
// Property specialization for bool type.
//-----------------------------------------------------------------------------

template<>
inline void PropertyStorageT<bool>::swap(size_t _i0, size_t _i1)
{
    // std::vector<bool>::swap(reference x, reference y) exists, but
    // on libstdc++ with _GLIBCXX_DEBUG it doesn't compile
    // (2018-02-26, libstdc++ 8.2.0)

    bool tmp = data_[_i0];
    data_[_i0] = data_[_i1];
    data_[_i1] = tmp;
}


template<>
inline void PropertyStorageT<bool>::deserialize(std::istream& _istr)
{
    for(unsigned int i = 0; i < size(); ++i) {
        value_type val;
        OpenVolumeMesh::deserialize(_istr, val);
        data_[i] = val;
    }
}
} // Namespace OpenVolumeMesh

