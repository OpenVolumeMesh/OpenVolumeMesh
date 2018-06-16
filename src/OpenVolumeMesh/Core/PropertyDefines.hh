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

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef PROPERTYDEFINES_HH_
#define PROPERTYDEFINES_HH_

#include <iosfwd>
#include <stdexcept>
#include <string>
#include <typeinfo>

#include "PropertyHandles.hh"
#include "PropertyPtr.hh"

namespace OpenVolumeMesh {

template <class T>
class OpenVolumeMeshPropertyT;

class ResourceManager;

template <class T>
const std::string typeName();

template <> const std::string typeName<int>();
template <> const std::string typeName<unsigned int>();
template <> const std::string typeName<short>();
template <> const std::string typeName<long>();
template <> const std::string typeName<unsigned long>();
template <> const std::string typeName<char>();
template <> const std::string typeName<unsigned char>();
template <> const std::string typeName<bool>();
template <> const std::string typeName<float>();
template <> const std::string typeName<double>();
template <> const std::string typeName<std::string>();
template <>  const std::string typeName<std::map<HalfEdgeHandle, int> >();
template <>  const std::string typeName<std::vector<double> >();
template <>  const std::string typeName<std::vector<VertexHandle> >();
template <>  const std::string typeName<std::vector<HalfFaceHandle> >();
template <>  const std::string typeName<std::vector<std::vector<HalfFaceHandle> > >();

/// Property classes for the different entity types
template<class T>
class VertexPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, VertexPropHandle> {
public:
    VertexPropertyT(const std::string& _name, ResourceManager& _resMan, VertexPropHandle _handle, const T _def = T());
    virtual ~VertexPropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "VProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    VertexPropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, VertexPropHandle _handle);
};
template<class T>
class EdgePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, EdgePropHandle> {
public:
    EdgePropertyT(const std::string& _name, ResourceManager& _resMan, EdgePropHandle _handle, const T _def = T());
    virtual ~EdgePropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "EProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    EdgePropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, EdgePropHandle _handle);
};
template<class T>
class HalfEdgePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfEdgePropHandle> {
public:
    HalfEdgePropertyT(const std::string& _name, ResourceManager& _resMan, HalfEdgePropHandle _handle, const T _def = T());
    virtual ~HalfEdgePropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "HEProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    HalfEdgePropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, HalfEdgePropHandle _handle);
};
template<class T>
class FacePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, FacePropHandle> {
public:
    FacePropertyT(const std::string& _name, ResourceManager& _resMan, FacePropHandle _handle, const T _def = T());
    virtual ~FacePropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "FProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    FacePropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, FacePropHandle _handle);
};
template<class T>
class HalfFacePropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, HalfFacePropHandle> {
public:
    HalfFacePropertyT(const std::string& _name, ResourceManager& _resMan, HalfFacePropHandle _handle, const T _def = T());
    virtual ~HalfFacePropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "HFProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    HalfFacePropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, HalfFacePropHandle _handle);
};
template<class T>
class CellPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, CellPropHandle> {
public:
    CellPropertyT(const std::string& _name, ResourceManager& _resMan, CellPropHandle _handle, const T _def = T());
    virtual ~CellPropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "CProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    CellPropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, CellPropHandle _handle);
};
template<class T>
class MeshPropertyT : public PropertyPtr<OpenVolumeMeshPropertyT<T>, MeshPropHandle> {
public:
    MeshPropertyT(const std::string& _name, ResourceManager& _resMan, MeshPropHandle _handle, const T _def = T());
    virtual ~MeshPropertyT() {}
    virtual BaseProperty* clone(ResourceManager &_resMan, const OpenVolumeMeshHandle _handle) const;
    virtual void serialize(std::ostream& _ostr) const;
    virtual void deserialize(std::istream& _istr);
    virtual const std::string entityType() const { return "MProp"; }
    virtual const std::string typeNameWrapper() const { return typeName<T>(); }
private:
    MeshPropertyT(OpenVolumeMeshPropertyT<T> *_prop, ResourceManager& _resMan, MeshPropHandle _handle);
};

} // Namespace OpenVolumeMesh

#if defined(INCLUDE_TEMPLATES) && !defined(PROPERTYDEFINEST_CC)
#include "PropertyDefinesT.cc"
#endif

#endif /* PROPERTYDEFINES_HH_ */
