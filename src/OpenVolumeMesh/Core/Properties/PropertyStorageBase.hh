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

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>
#include <OpenVolumeMesh/Core/detail/Tracking.hh>

namespace OpenVolumeMesh {

class ResourceManager;
class BasePropertyPtr;

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
    PropertyStorageBase(
            detail::Tracker<PropertyStorageBase> *tracker,
            std::string _name,
            std::string _internal_type_name,
            EntityType _entity_type,
            bool _shared)
        : detail::Tracked<PropertyStorageBase>(tracker)
        , name_(std::move(_name))
        , internal_type_name_(std::move(_internal_type_name))
        , entity_type_(_entity_type)
        , persistent_(false)
        , shared_(_shared)
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

    /// Copy property values of src_idx to dst_idx
    virtual void copy(size_t src_idx, size_t dst_idx) = 0;

	/// Erase an element of the vector
	virtual void delete_element(size_t _idx) = 0;

	/// Return a deep copy of self.
    virtual std::shared_ptr<PropertyStorageBase> clone() const = 0;

	const std::string& name() const && = delete;
    /// Return the name of the property
    const std::string& name() const & {
		return name_;
	}
    bool anonymous() const {return name_.empty();}

    void set_name(std::string _name) {
        name_ = std::move(_name);
    }

	const std::string& internal_type_name() const && = delete;

    /// A platform-specific unique name for the contained type.
	const std::string& internal_type_name() const & {
		return internal_type_name_;
	}

    // ovm-ascii file name
    virtual std::string typeNameWrapper() const = 0;

	// Function to serialize a property
	virtual void serialize(std::ostream& /*_ostr*/) const {}

	// Function to deserialize a property
    virtual void deserialize(std::istream& /*_istr*/) {}
	// I/O support

	void set_persistent(bool _persistent) { persistent_ = _persistent; }
    void set_shared(bool _shared) { shared_ = _shared; }

	bool persistent() const { return persistent_; }
    bool shared() const { return shared_; }

    EntityType entity_type() const {return entity_type_;}

    template<typename T>
    PropertyStorageT<T>* cast_to_StorageT()
    {
        if (detail::internal_type_name<T>() != internal_type_name()) {
            throw std::bad_cast();
        }
        return static_cast<PropertyStorageT<T>*>(this);
    }

    template<typename T>
    const PropertyStorageT<T>* cast_to_StorageT() const
    {
        if (detail::internal_type_name<T>() != internal_type_name()) {
            throw std::bad_cast();
        }
        return static_cast<const PropertyStorageT<T>*>(this);
    }

    virtual std::unique_ptr<BasePropertyPtr> make_property_ptr() = 0;

    /// Copy data from other property. `other` must point to an object with the same derived type as `this`!
    virtual void assign_values_from(const PropertyStorageBase *other) = 0;

    /// Move data from other property. `other` must point to an object with the same derived type as `this`!
    virtual void move_values_from(PropertyStorageBase *other) = 0;

    /// Attach to a mesh. Does not resize property!
    void attach_to(const ResourceManager *resman);

    /// Is this property storage attached to a mesh?
    operator bool() const {return Tracked::has_tracker();}

private:
	std::string name_;
    // while we could compute this on the fly,
    // it does not take much storage and is useful in debugging
    std::string internal_type_name_;
    EntityType entity_type_;


    // TODO: an enum might be nicer for this:
    bool persistent_;
    bool shared_;
};

} // Namespace OpenVolumeMesh

