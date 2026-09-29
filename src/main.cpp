#include "quad_levelset_mesher.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    using namespace split_test;

    QuadInput input{
        {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}}},
        {{-0.35, 0.70, 0.55, 0.80}}
    };

    std::string vtk_file = "split_mesh.vtk";
    if (argc == 5 || argc == 6) {
        for (int i = 0; i < 4; ++i) input.phi[static_cast<std::size_t>(i)] = std::atof(argv[i + 1]);
        if (argc == 6) vtk_file = argv[5];
    } else if (argc != 1) {
        std::cerr << "Usage: " << argv[0] << " [phi0 phi1 phi2 phi3 [output.vtk]]\n";
        return 2;
    }

    try {
        const QuadLevelSetMesher mesher;
        const MixedMesh mesh = mesher.remesh(input);
        write_legacy_vtk(vtk_file, mesh);

        std::cout << "nodes      : " << mesh.nodes.size() << '\n';
        std::cout << "cells      : " << mesh.cells.size() << '\n';
        std::cout << "quads      : " << mesh.quad_count() << '\n';
        std::cout << "triangles  : " << mesh.triangle_count() << '\n';
        std::cout << "total area : " << mesh.total_area() << '\n';
        std::cout << "VTK        : " << vtk_file << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
