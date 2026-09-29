# split_test

C++17 prototype for cutting one convex Q1 quadrilateral with nodal level-set values and remeshing the two phases with a **quad-dominant mixed triangle/quad mesh**.

## Features

- Takes 4 quadrilateral corner coordinates and 4 nodal `phi` values.
- Uses linear edge interpolation to locate `phi = 0`.
- Treats the ordinary two-edge intersection as a straight interface chord.
- Clips the element into positive/negative phase polygons.
- Keeps 4-sided regions as quads.
- Splits a 5-sided region into one quad + one triangle using a simple quality criterion.
- Handles the alternating-sign Q1 saddle case with center-sign disambiguation.
- Exports legacy ASCII VTK unstructured grids containing:
  - triangle and quad cells,
  - nodal `phi`,
  - cell `phase` (-1/+1),
  - cell `element_type` (3/4).

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Run

Default example:

```bash
./build/split_demo
```

Custom nodal values:

```bash
./build/split_demo -0.35 0.70 0.55 0.80 my_cut.vtk
```

Open the generated VTK file in ParaView.

## Library use

```cpp
#include "quad_levelset_mesher.hpp"

split_test::QuadInput q{
    {{{0,0}, {1,0}, {1,1}, {0,1}}},
    {{-1.0, 0.5, 1.0, 0.2}}
};

split_test::QuadLevelSetMesher mesher;
auto mesh = mesher.remesh(q);
split_test::write_legacy_vtk("cut.vtk", mesh);
```

The prototype intentionally targets a **single convex Q1 quadrilateral**. For a global mesh, the next step is to apply this kernel element-by-element and deduplicate interface nodes using the parent edge identity.
