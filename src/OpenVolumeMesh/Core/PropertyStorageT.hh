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

#include <OpenVolumeMesh/Core/PropertyStorageBase.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
//#include <OpenVolumeMesh/Core/detail/Tracking.hh>

#include <OpenVolumeMesh/Core/Serializers.hh>

namespace OpenVolumeMesh {


template <class T>
class PropertyStorageT;

template <typename T>
/// convenience access to properties through a shared_ptr
class PropertyStoragePtr
                         //, public detail::Tracked<PropertyStoragePtr<T>>
{
public:
    using PropStorageT = PropertyStorageT<T>;
    typedef typename PropStorageT::value_type                  value_type;
    typedef typename PropStorageT::vector_type::const_iterator const_iterator;
    typedef typename PropStorageT::vector_type::iterator       iterator;
    typedef typename PropStorageT::reference                   reference;
    typedef typename PropStorageT::const_reference             const_reference;


    virtual ~PropertyStoragePtr();

    const_iterator begin() const { return storage()->begin(); }
    const_iterator end() const   { return storage()->end(); }
    iterator begin() { return storage()->begin(); }
    iterator end()   { return storage()->end(); }
    size_t size() const { return storage()->size(); }
    operator bool() const {return storage()->resMan() != nullptr;}


    /// No range check performed!
    reference operator[](size_t idx) { return (*storage())[idx];};
    /// No range check performed!
    const_reference operator[](size_t idx) const { return (*storage())[idx];};
    reference at(size_t idx) { return storage()->at(idx);};
    const_reference at(size_t idx) const { return storage()->at(idx);};

    void serialize(std::ostream& _ostr) const { storage()->serialize(_ostr); }
    void deserialize(std::istream& _istr) { storage()->deserialize(_istr); }

    bool persistent() const { return storage()->persistent(); }
    bool anonymous() const { return storage()->anonymous(); }
    size_t n_elements() const { return storage()->n_elements(); }

    std::string typeNameWrapper() const {return storage()->typeNameWrapper(); }
    EntityType entity_type() const {return storage()->entity_type();}

    const std::string& name() const & {
        // the string we return a reference to lives long enough, no warnings please:
        // cppcheck-suppress returnTempReference
        return storage()->name();
    }

    /// get default value.
    T const &def() const {return storage()->def();}

    /// set all values to `val`.
    void fill(T const&val) { storage()->fill(val); }



protected:
    PropertyStoragePtr(std::shared_ptr<PropertyStorageT<T>> &&_ptr = nullptr);

    std::shared_ptr<PropertyStorageT<T>> &storage();
    std::shared_ptr<PropertyStorageT<T>> const &storage() const;

    friend PropertyStorageT<T>;
    void set_storage(std::shared_ptr<PropertyStorageT<T>> &&_ptr);


private:
    std::shared_ptr<PropertyStorageT<T>> storage_;
};
template <class T, typename EntityTag>
class PropertyPtr;

//== CLASS DEFINITION =========================================================

/** \class PropertyStorageT
 *
 *  \brief Default property class for any type T.
 *
 *  The default property class for any type T.
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
            const std::string& _name,
            EntityType _entity_type,
            const T &_def = T())
        : PropertyStorageBase(tracker, _name, get_type_name<T>(), _entity_type),
          def_(_def)
    {}

public:
    void reserve(size_t _n) override{
		data_.reserve(_n);
	}
	void resize(size_t _n) override {
                data_.resize(_n, def_);
	}
	size_t size() const override {
		return data_.size();
	}
	void clear() override {
		data_.clear();
	}
	void push_back() override {
		data_.push_back(def_);
	}
	void swap(size_t _i0, size_t _i1) override {
        std::swap(data_[_i0], data_[_i1]);
    }
    void copy(size_t _src_idx, size_t _dst_idx) override{
		data_[_dst_idx] = data_[_src_idx];
	}
	void delete_element(size_t _idx) override {
        assert(_idx < data_.size());
		data_.erase(data_.begin() + static_cast<long>(_idx));
	}

public:

	size_t n_elements() const override {
		return data_.size();
	}

	// Function to serialize a property
    void serialize(std::ostream& _ostr) const override {
        for(typename vector_type::const_iterator it = data_.begin();
                it != data_.end(); ++it) {
            OpenVolumeMesh::serialize(_ostr, *it) << std::endl;
        }
    }

    // Function to deserialize a property
    void deserialize(std::istream& _istr) override {
        for(unsigned int i = 0; i < n_elements(); ++i) {
            OpenVolumeMesh::deserialize(_istr, data_[i]);
        }
    }

public:
	// data access interface

	/// Get pointer to array (does not work for T==bool)
	const T* data() const {

		if (data_.empty())
			return 0;

		return &data_[0];
	}

	/// Get reference to property vector (be careful, improper usage, e.g. resizing, may crash)
	vector_type& data_vector() {

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

    std::shared_ptr<PropertyStorageBase> clone() const override {
        auto res = std::make_shared<PropertyStorageT<T>>(*this);
        res->set_tracker(nullptr);
        return res;
	}

	typename vector_type::const_iterator begin() const { return data_.begin(); }

	typename vector_type::iterator begin() { return data_.begin(); }

	typename vector_type::const_iterator end() const { return data_.end(); }

    typename vector_type::iterator end() { return data_.end(); }

    std::string typeNameWrapper() const override {return OpenVolumeMesh::typeName<T>();}

    T const& def() const {return def_;}

    void fill(T const&val) {
        std::fill(data_.begin(), data_.end(), val);
    }


    void assign_values_from(const PropertyStorageBase *_other) override
    {
        const auto *other = _other->cast_to_StorageT<T>();
        data_ = other->data_;
        def_ = other->def_;

    }

    void move_values_from(PropertyStorageBase *_other) override
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

    auto tmp = data_[_i0];
    data_[_i0] = data_[_i1];
    data_[_i1] = tmp;
}


template<>
inline void PropertyStorageT<bool>::deserialize(std::istream& _istr)
{
    for(unsigned int i = 0; i < n_elements(); ++i) {
        value_type val;
        OpenVolumeMesh::deserialize(_istr, val);
        data_[i] = val;
    }
}
template<typename T>
PropertyStoragePtr<T>::PropertyStoragePtr(std::shared_ptr<PropertyStorageT<T> > &&_ptr)
    :
      //detail::Tracked<PropertyStoragePtr<T>>(&_ptr->pointer_tracker) ,
      storage_(std::move(_ptr))
{}

template<typename T>
void PropertyStoragePtr<T>::set_storage(std::shared_ptr<PropertyStorageT<T> > &&_ptr)
{
    storage_ = std::move(_ptr);
}

template<typename T>
PropertyStoragePtr<T>::~PropertyStoragePtr() = default;

template<typename T>
std::shared_ptr<PropertyStorageT<T> > &PropertyStoragePtr<T>::storage() {return storage_;}

template<typename T>
const std::shared_ptr<PropertyStorageT<T> > &PropertyStoragePtr<T>::storage() const {return storage_;}


} // Namespace OpenVolumeMesh

