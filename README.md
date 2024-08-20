# OpenVolumeMesh

![OpenVolumeMesh Logo](/documentation/images/OpenVolumeMesh_text_128.png)

A generic C++ data structure for arbitrary polyhedral meshes,
with specializations for tetrahedral and hexahedral meshes.

[[_TOC_]]

## Introduction

The concepts of OVM are closely related to [OpenMesh](http://www.openmesh.org).
In particular, OpenVolumeMesh carries the general
idea of storing edges as so-called (directed) half-edges over to the face
definitions. So, faces are split up into so-called half-faces having opposing
orientations. Furthermore, in OpenVolumeMesh the data is arranged in a top-down
hierarchy, meaning that each entity of dimension n is defined through a
(ordered) tuple of entities of dimension (n-1). These are the intrinsic
adjacency relations of the volumentric meshes. One can additionally compute
bottom-up adjacencies which means that for each entity of dimension n, we also
store its local adjacencies to entities of dimension (n+1). These adjacency
relations have to be computed explicitly which can be performed in linear time
complexity. Both adjacency relations, the top-down and the bottom-up
adjacencies, are used to provide a set of iterators and circulators that are
comfortable in use. As in OpenMesh, OpenVolumeMesh provides an entirely generic
underlying property system that allows attaching properties of any kind to the
entities. In order to learn more about the implementational details of
OpenVolumeMesh, please refer to the online documentation available
at <http://openvolumemesh.org/Documentation/OpenVolumeMesh-Doc-Latest>.

OpenVolumeMesh is entirely written in C++, making heavy use of the
standard template library as well as template programming paradigms.
Although OpenVolumeMesh has been developed to the best of my knowledge,
it does not claim to be free from defects nor does it raises the claim to
have inveterate underlying implemented concepts. So, any ambitious developer
is invited to participate in the development process to make OpenVolumeMesh
a well-working, reliable, and useful library.

The initial version of OpenVolumeMesh is described in detail in [
Kremer, M., Bommes, D., Kobbelt, L. (2013). *OpenVolumeMesh – A Versatile Index-Based Data Structure for 3D Polytopal Complexes*. IMR 2012](https://graphics.rwth-aachen.de/media/papers/MKremer_OVM.pdf)


## Building and linking OpenVolumeMesh

You can build and install OpenVolumeMesh as static or shared library which can be linked
into C++ projects using any build system.

If your project uses CMake, you can integrate OpenVolumeMesh into
its build either using `add_subdirectory` to build OVM together with your project,
or using `find_package(OpenVolumeMesh)` to find a previously-installed standalone
library build.

Either option produces an `OpenVolumeMesh::OpenVolumeMesh` target.

### Build requirements

* [CMake](http://www.cmake.org) 3.12 or newer
* A C++17-compatible compiler (tested with GCC, Clang and MSVC)
* (Optional) [Doxyen](http://www.doxygen.org) for building the reference documentation
* (Optional, downloaded by CMake) [Google Test](https://github.com/google/googletest.git) for building the unit test suite.

### Command-line build

- Create a new build directory, e.g. "OVM-build-release", or whatever name might be suitable.

- Change into the recently created directory and type
  "cmake /path/to/OpenVolumeMesh/sources".

- If you want to change the build configuration, say from debug to release,
  type `ccmake .` or `cmake-gui .`.

- Once everything is configured to your satisfaction, type `cmake --build .` followed
  by `sudo cmake --build . -t install ` in order to build and install the library.

### MacOS/XCode build

- You can use cmake to generate XCode project files

### Windows / Visual Studio build

- You can either directly open the project using Visual Studio's CMake integration, or generate project files by using cmake directly.

#### Using Visual Studio's CMake integration
- TODO

#### Generating project files

On Windows, start the CMake gui tool and set the path to OpenVolumeMesh's
sources. Then select the target project type (Visual Studio 20xx)
and click on "Configure". Once everything is configured to your satisfaction,
click on "Generate". You will now find a Visual Studio project file
in the specified build folder (which is "Build" per default). Open this
file in Visual Studio and select "Build all".

## Concepts

### Supported mesh types
- Manifold volumetric meshes, limited support for non-manifold configurations
- self-adjacent entities

TODO: define exactly what configurations are allowed

### Entities

OpenVolumeMesh Entities are comprised of
- Vertices
- Edges and half-edges (directed edges), made of two vertices
- Faces and half-faces (directed faces), made up of a sequence of half-edges
- Cells, consisting of a set of half-faces

![Half-entities diagram](/documentation/images/halfedge_halfface.png)
![Entity Hierarchy diagram](/documentation/images/volume_mesh_hierarchy.png)

Entities are always stored in a top-down manner, i.e., entities of dimension $n$
are defined by their constituient entities of dimension $n-1$.

For every entity, there is a `Handle` type, which is a type-safe way of indexing meshes.

### Bottom-Up Incidence relations

For efficient mesh nagivation, OpenVolumMesh also provides bottom-up incidence
relations, i.e., iterators for higher-dimensional entities.
This support is optional: The different bottom-up incidences can be disabled
for lower memory usage and higher performance.

### Properties

OpenVolumeMesh supports storing data with any entity type using properties.

TODO: show examples

### Attributes

TODO

### File I/O

- `.vtk`
- `.ovm`: ASCII based file format
- `.ovmb`: efficient binary file format with support for custom property (de)serialization

## Contributing

The main repository is located at
[VCI Aachen gitlab](https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh),
which is read-only for the general public.
We run an official [github mirror](https://github.com/OpenVolumeMesh/OpenVolumeMesh) where you
can open issues, start discussions and submit pull requests.

Alternatively you can e-mail us at <martin.heistermann@unibe.ch> and <moebius@cs.rwth-aachen.de>.

## Related projects

Notable users of OVM include:
- [OpenFlipper](https://openflipper.org) - Interactive framework and GUI for geometry processing, supports visualisation, generation and processing of volumetric meshes using OVM
- [libHexEx](https://www.graphics.rwth-aachen.de/software/libHexEx/) - Robust extraction of hexahedral meshes from integer-grid maps
- [AlgoHex](https://www.algohex.eu/) (Automated hexahedral meshing): [libAlgoHex](https://github.com/cgg-bern/AlgoHex)
- Volume Parametrization Quantization for Hexahedral Meshing:
    - [TrulySeamless3D](https://github.com/HendrikBrueckler/TrulySeamless3D): Numerical sanitization of volumetric seamless maps
    - [MC3D](https://github.com/HendrikBrueckler/MC3D): Structured Volume Decomposition
    - [QGP3D](https://github.com/HendrikBrueckler/QGP3D): Quantization for volumetric motorcycle complexes

Other C++ libraries that aid in handling unstructured volumetric meshes include:

- [CGoGN](https://cgogn.github.io/) - combinatorial maps
- [CGAL Combinatorial Maps](https://doc.cgal.org/latest/Combinatorial_map/index.html)
- [CinoLib](https://github.com/mlivesu/cinolib)
- [MTet: A Mini Tetrahedral Mesh Data Structure](https://github.com/qnzhou/mtet)

## License Information

OpenVolumeMesh is free software licensed under the terms of the
GNU Lesser General Public License Version 3 as published by the Free Software
Foundation. You can redistribute and/or modify it as stated in the
above mentioned license terms. A copy of the license can be found
in the license sub-folder of this source-tree or under
<http://www.gnu.org/licenses>. By downloading and using the OpenVolumeMesh
library you automatically agree to these terms.

