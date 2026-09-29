#include "quad_levelset_mesher.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace split_test;

namespace {
QuadInput unit_quad(std::array<double, 4> phi) {
    return {{{{0.0,0.0}, {1.0,0.0}, {1.0,1.0}, {0.0,1.0}}}, phi};
}

void check_area(const MixedMesh& mesh, double expected = 1.0) {
    assert(std::abs(mesh.total_area() - expected) < 1.0e-10);
}
}

int main() {
    const QuadLevelSetMesher mesher;

    {
        const auto m = mesher.remesh(unit_quad({1.0, 2.0, 3.0, 4.0}));
        assert(m.cells.size() == 1);
        assert(m.quad_count() == 1);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({1.0, 2.0, 3.0, 4.0}), 3, RemeshMode::QuadDominant);
        assert(m.cells.size() == 9);
        assert(m.quad_count() == 9);
        assert(m.triangle_count() == 0);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({1.0, 2.0, 3.0, 4.0}), 3, RemeshMode::TriangleOnly);
        assert(m.cells.size() == 18);
        assert(m.quad_count() == 0);
        assert(m.triangle_count() == 18);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({-1.0, 1.0, 1.0, -1.0}), 2, RemeshMode::QuadDominant);
        assert(m.cells.size() == 8);
        assert(m.quad_count() == 8);
        assert(m.triangle_count() == 0);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({-1.0, 1.0, 1.0, 1.0}), 3, RemeshMode::QuadDominant);
        assert(m.cells.size() == 27);
        assert(m.quad_count() == 9);
        assert(m.triangle_count() == 18);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({-1.0, 1.0, 1.0, 1.0}), 3, RemeshMode::TriangleOnly);
        assert(m.quad_count() == 0);
        assert(m.triangle_count() == 36);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({1.0, -0.8, 1.2, -1.4}), 4, RemeshMode::TriangleOnly);
        assert(!m.cells.empty());
        assert(m.quad_count() == 0);
        check_area(m);
    }

    {
        bool threw = false;
        try {
            (void)mesher.remesh(unit_quad({-1.0, 1.0, 1.0, -1.0}), 0);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    std::cout << "all tests passed\n";
    return 0;
}
