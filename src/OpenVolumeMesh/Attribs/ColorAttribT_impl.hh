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


#define COLORATTRIBT_CC

#include <OpenVolumeMesh/Attribs/ColorAttrib.hh>

namespace OpenVolumeMesh {

template <class ColT>
ColorAttrib<ColT>::ColorAttrib(TopologyKernel& _kernel, const ColT _def) :
        vcolor_prop_(_kernel.request_vertex_property<ColT>("vertex_color", _def)),
        ecolor_prop_(_kernel.request_edge_property<ColT>("edge_color", _def)),
        hecolor_prop_(_kernel.request_halfedge_property<ColT>("halfedge_color", _def)),
        fcolor_prop_(_kernel.request_face_property<ColT>("face_color", _def)),
        hfcolor_prop_(_kernel.request_halfface_property<ColT>("halfface_color", _def)),
        ccolor_prop_(_kernel.request_cell_property<ColT>("cell_color", _def)),
        vertex_colors_available_(false),
        halfedge_colors_available_(false),
        edge_colors_available_(false),
        halfface_colors_available_(false),
        face_colors_available_(false),
        cell_colors_available_(false),
        default_color_ (_def)
{

}

template <class ColT>
ColorAttrib<ColT>::~ColorAttrib() {

}

template <class ColT>
void ColorAttrib<ColT>::clear_vertex_colors()
{
    vcolor_prop_.fill(default_color_);
    vertex_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_halfedge_colors()
{
    hecolor_prop_.fill(default_color_);
    halfedge_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_edge_colors()
{
    ecolor_prop_.fill(default_color_);
    edge_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_halfface_colors()
{
    hfcolor_prop_.fill(default_color_);
    halfface_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_face_colors()
{
    fcolor_prop_.fill(default_color_);
    face_colors_available_   = false;
}

template <class ColT>
void ColorAttrib<ColT>::clear_cell_colors()
{
    ccolor_prop_.fill(default_color_);
    cell_colors_available_   = false;
}

} // Namespace OpenVolumeMesh
