#include "quad_levelset_mesher.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace split_test;

namespace {
QuadInput unit_quad(std::array<double, 4> phi) {
    return {{{{0.0,0.0}, {1.0,0.0}, {1.0,1.0}, {0.0,1.0}}}, phi};
}

void check_area(const MixedMesh& mesh, double expected = 1.0) {
    assert(std::abs(mesh.total_area() - expected) < 1.0e-11);
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
        const auto m = mesher.remesh(unit_quad({-1.0, 1.0, 1.0, -1.0}));
        assert(m.cells.size() == 2);
        assert(m.quad_count() == 2);
        assert(m.triangle_count() == 0);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({-1.0, 1.0, 1.0, 1.0}));
        assert(m.cells.size() == 3);
        assert(m.quad_count() == 1);
        assert(m.triangle_count() == 2);
        check_area(m);
    }

    {
        const auto m = mesher.remesh(unit_quad({1.0, -0.8, 1.2, -1.4}));
        assert(!m.cells.empty());
        check_area(m);
    }

    std::cout << "all tests passed\n";
    return 0;
}
