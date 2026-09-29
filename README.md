# split_test

C++17 prototype for cutting one convex Q1 quadrilateral with nodal level-set values and remeshing the two phases.

It now supports two output states:

1. **Quad-dominant** mixed triangle/quad mesh.
2. **Triangle-only** mesh.

It also supports a user-controlled subdivision count.

## Division control

`divisions = N` means:

- every coarse quad becomes `N x N` quads;
- every coarse triangle becomes `N^2` triangles.

So increasing `N` gives a denser mesh without moving the straight level-set interface.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Run

Default example uses `divisions = 4` and writes both VTK states:

```bash
./build/split_demo
```

Outputs:

```text
split_mesh_quad_dominant.vtk
split_mesh_triangles.vtk
```

Custom level-set values, division count, and output prefix:

```bash
./build/split_demo -0.35 0.70 0.55 0.80 8 my_cut
```

Outputs:

```text
my_cut_quad_dominant.vtk
my_cut_triangles.vtk
```

## Library use

```cpp
#include "quad_levelset_mesher.hpp"

split_test::QuadInput q{
    {{{0,0}, {1,0}, {1,1}, {0,1}}},
    {{-1.0, 0.5, 1.0, 0.2}}
};

split_test::QuadLevelSetMesher mesher;

const std::size_t divisions = 6;

auto quad_mesh = mesher.remesh(
    q, divisions, split_test::RemeshMode::QuadDominant);

auto tri_mesh = mesher.remesh(
    q, divisions, split_test::RemeshMode::TriangleOnly);

split_test::write_legacy_vtk("quad_dominant.vtk", quad_mesh);
split_test::write_legacy_vtk("triangles.vtk", tri_mesh);
```

The VTK files contain nodal `phi`, cell `phase` (-1/+1), and `element_type` (3/4).
