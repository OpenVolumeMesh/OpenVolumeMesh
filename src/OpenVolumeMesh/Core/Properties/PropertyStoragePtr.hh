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


#include <istream>
#include <ostream>
#include <numeric>
#include <string>
#include <memory>

#include <OpenVolumeMesh/Core/Properties/PropertyStorageBase.hh>

namespace OpenVolumeMesh {


template <typename T>
class PropertyStorageT;
namespace IO {
template <typename T, typename Codec>
class PropertyDecoderT;
} // namespace IO

template <typename T>
class PropertyStoragePtr
{
public:
    using PropStorageT = PropertyStorageT<T>;
    typedef typename PropStorageT::value_type                  value_type;
    typedef typename PropStorageT::vector_type::const_iterator const_iterator;
    typedef typename PropStorageT::vector_type::iterator       iterator;
    typedef typename PropStorageT::reference                   reference;
    typedef typename PropStorageT::const_reference             const_reference;


    virtual ~PropertyStoragePtr() = default;

    std::vector<T> const& data_vector() const { return storage()->data_vector(); }
    const_iterator begin() const { return storage()->begin(); }
    const_iterator end() const   { return storage()->end(); }
    iterator begin() { return storage()->begin(); }
    iterator end()   { return storage()->end(); }
    size_t size() const { return storage()->size(); }

    operator bool() const {return *storage();}

    reference operator[](size_t idx) { return (*storage())[idx];};
    const_reference operator[](size_t idx) const { return (*storage())[idx];};

    reference at(size_t idx) { return storage()->at(idx);};
    const_reference at(size_t idx) const { return storage()->at(idx);};

    // OVM-ASCII file IO:
    std::string typeNameWrapper() const {return storage()->typeNameWrapper(); }
    void serialize(std::ostream& _ostr) const { storage()->serialize(_ostr); }
    void deserialize(std::istream& _istr) { storage()->deserialize(_istr); }

    bool persistent() const { return storage()->persistent(); }
    bool shared() const { return storage()->shared(); }
    bool anonymous() const { return storage()->anonymous(); }

    EntityType entity_type() const {return storage()->entity_type();}

    const std::string& name() const & {
        // the string we return a reference to lives long enough, no warnings please:
        // cppcheck-suppress returnTempReference
        return storage()->name();
    }
    void set_name(std::string _name) {
        storage()->set_name(std::move(_name));
    }

    /// get default value.
    T const &def() const {return storage()->def();}

    /// set all values to `val`.
    void fill(T const&val) { storage()->fill(val); }

    void swap(std::vector<T> &_other) { storage()->swap(_other); }

protected:
    friend PropertyStorageT<T>;

    template<typename _T, typename Codec>
    friend class IO::PropertyDecoderT; // for storage() access

    PropertyStoragePtr(std::shared_ptr<PropertyStorageT<T>> _ptr = nullptr)
        : storage_(std::move(_ptr))
    {}

    std::shared_ptr<PropertyStorageT<T>> &storage() {return storage_;}
    std::shared_ptr<PropertyStorageT<T>> const &storage() const {return storage_;}

    void set_storage(std::shared_ptr<PropertyStorageT<T>> _ptr)
    {
        storage_ = std::move(_ptr);
    }

private:
    std::shared_ptr<PropertyStorageT<T>> storage_;
};












} // Namespace OpenVolumeMesh
