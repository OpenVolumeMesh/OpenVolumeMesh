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

#include <iosfwd>
#include <string>
#include <memory>
#include <vector>

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/TypeName.hh>
#include <OpenVolumeMesh/Core/detail/Tracking.hh>

namespace OpenVolumeMesh {


template<typename T>
class PropertyStorageT;

/** \class PropertyStorageBase

 Abstract class defining the type-independent interface of a dynamic property.

 **/

class OVM_EXPORT PropertyStorageBase
        : public std::enable_shared_from_this<PropertyStorageBase>
        , public detail::Tracked<PropertyStorageBase>
{
public:
    template <class PropT, class HandleT> friend class PropertyPtr;

public:

    explicit PropertyStorageBase(
            detail::Tracker<PropertyStorageBase> *tracker,
            const std::string& _name,
            const std::string& _internal_type_name,
            EntityType _entity_type)
        : detail::Tracked<PropertyStorageBase>(tracker)
        , name_(_name)
        , internal_type_name_(_internal_type_name)
        , entity_type_(_entity_type)
        , persistent_(false)
    {}

    virtual ~PropertyStorageBase() = default;


public:

	/// Reserve memory for n elements.
	virtual void reserve(size_t _n) = 0;

	/// Resize storage to hold n elements.
	virtual void resize(size_t _n) = 0;

	/// Return underlying container size
	virtual size_t size() const = 0;

	/// Clear all elements and free memory.
	virtual void clear() = 0;

	/// Extend the number of elements by one.
	virtual void push_back() = 0;

	/// Let two elements swap their storage place.
	virtual void swap(size_t _i0, size_t _i1) = 0;

    /// Let two elements swap their storage place.
    virtual void copy(size_t _i0, size_t _i1) = 0;

	/// Erase an element of the vector
	virtual void delete_element(size_t _idx) = 0;

	/// Return a deep copy of self.
    virtual std::shared_ptr<PropertyStorageBase> clone() const = 0;

	/// Return the name of the property
	const std::string& name() const && = delete;
	const std::string& name() const & {
		return name_;
	}
    bool anonymous() const {return name_.empty();}

	const std::string& internal_type_name() const && = delete;
	const std::string& internal_type_name() const & {
		return internal_type_name_;
	}

	// Function to serialize a property
	virtual void serialize(std::ostream& /*_ostr*/) const {}

	// Function to deserialize a property
    virtual void deserialize(std::istream& /*_istr*/) {}
	// I/O support

	void set_persistent(bool _persistent) { persistent_ = _persistent; }

	bool persistent() const { return persistent_; }

	/// Number of elements in property
	virtual size_t n_elements() const = 0;

    virtual std::string typeNameWrapper() const = 0;

    EntityType entity_type() const {return entity_type_;}

    template<typename T>
    PropertyStorageT<T>* cast_to_StorageT()
    {
        if (get_type_name<T>() != internal_type_name()) {
            throw std::bad_cast();
        }
        return static_cast<PropertyStorageT<T>*>(this);
    }

    template<typename T>
    const PropertyStorageT<T>* cast_to_StorageT() const
    {
        if (get_type_name<T>() != internal_type_name()) {
            throw std::bad_cast();
        }
        return static_cast<const PropertyStorageT<T>*>(this);
    }

    /// Copy data from other property. `other` must point to an object with the same derived type as `this`!
    virtual void assign_values_from(const PropertyStorageBase *other) = 0;

    /// Move data from other property. `other` must point to an object with the same derived type as `this`!
    virtual void move_values_from(PropertyStorageBase *other) = 0;

private:
	std::string name_;
	std::string internal_type_name_;
    EntityType entity_type_;
	bool persistent_;
};

} // Namespace OpenVolumeMesh

