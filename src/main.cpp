#include "quad_levelset_mesher.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    using namespace split_test;

    QuadInput input{
        {{{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}}},
        {{-0.35, 0.70, 0.55, 0.80}}
    };

    std::size_t divisions = 4;
    std::string prefix = "split_mesh";

    if (argc == 5 || argc == 6 || argc == 7) {
        for (int i = 0; i < 4; ++i) input.phi[static_cast<std::size_t>(i)] = std::atof(argv[i + 1]);
        if (argc >= 6) {
            const long parsed = std::strtol(argv[5], nullptr, 10);
            if (parsed < 1) throw std::invalid_argument("divisions must be >= 1");
            divisions = static_cast<std::size_t>(parsed);
        }
        if (argc == 7) prefix = argv[6];
    } else if (argc != 1) {
        std::cerr << "Usage: " << argv[0]
                  << " [phi0 phi1 phi2 phi3 [divisions [output_prefix]]]\n";
        return 2;
    }

    try {
        const QuadLevelSetMesher mesher;
        const MixedMesh quad_mesh = mesher.remesh(input, divisions, RemeshMode::QuadDominant);
        const MixedMesh tri_mesh  = mesher.remesh(input, divisions, RemeshMode::TriangleOnly);

        const std::string quad_file = prefix + "_quad_dominant.vtk";
        const std::string tri_file  = prefix + "_triangles.vtk";
        write_legacy_vtk(quad_file, quad_mesh);
        write_legacy_vtk(tri_file, tri_mesh);

        std::cout << "divisions: " << divisions << "\n\n";
        std::cout << "quad-dominant mesh\n";
        std::cout << "  nodes     : " << quad_mesh.nodes.size() << '\n';
        std::cout << "  cells     : " << quad_mesh.cells.size() << '\n';
        std::cout << "  quads     : " << quad_mesh.quad_count() << '\n';
        std::cout << "  triangles : " << quad_mesh.triangle_count() << '\n';
        std::cout << "  area      : " << quad_mesh.total_area() << '\n';
        std::cout << "  VTK       : " << quad_file << "\n\n";

        std::cout << "triangle-only mesh\n";
        std::cout << "  nodes     : " << tri_mesh.nodes.size() << '\n';
        std::cout << "  cells     : " << tri_mesh.cells.size() << '\n';
        std::cout << "  quads     : " << tri_mesh.quad_count() << '\n';
        std::cout << "  triangles : " << tri_mesh.triangle_count() << '\n';
        std::cout << "  area      : " << tri_mesh.total_area() << '\n';
        std::cout << "  VTK       : " << tri_file << '\n';
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
