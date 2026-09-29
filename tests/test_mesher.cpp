#include "quad_levelset_mesher.hpp"

#include <algorithm>
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

double cell_area(const MixedMesh& mesh, const Cell& cell) {
    double a = 0.0;
    for (std::size_t i = 0; i < cell.nodes.size(); ++i) {
        const auto& p = mesh.nodes[cell.nodes[i]].p;
        const auto& q = mesh.nodes[cell.nodes[(i + 1) % cell.nodes.size()]].p;
        a += p.x*q.y - q.x*p.y;
    }
    return 0.5 * std::abs(a);
}

double max_cell_area(const MixedMesh& mesh) {
    double a = 0.0;
    for (const auto& cell : mesh.cells) {
        a = std::max(a, cell_area(mesh, cell));
    }
    return a;
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
        // A triangular phase region in QuadDominant mode must not remain
        // entirely triangular: retain one inner triangle and surround it
        // with three quads.
        const auto m = mesher.remesh(
            unit_quad({-1.0, 1.0, 1.0, 1.0}),
            1, RemeshMode::QuadDominant);

        std::size_t negative_quads = 0;
        std::size_t negative_triangles = 0;
        for (const auto& cell : m.cells) {
            if (cell.phase != -1) continue;
            if (cell.type == CellType::Quad) ++negative_quads;
            if (cell.type == CellType::Triangle) ++negative_triangles;
        }

        assert(negative_quads == 3);
        assert(negative_triangles == 1);
        assert(m.quad_count() > m.triangle_count());
        check_area(m);
    }

    {
        // The interface x=0.5 lies exactly on the 2x2 subdivision line.
        const auto m = mesher.remesh(
            unit_quad({-1.0, 1.0, 1.0, -1.0}),
            2, RemeshMode::QuadDominant);
        assert(m.cells.size() == 4);
        assert(m.quad_count() == 4);
        assert(m.triangle_count() == 0);
        check_area(m);
        assert(max_cell_area(m) <= 0.25 + 1.0e-12);
    }

    {
        // Refinement happens before cutting: no final element may be larger
        // than one parent 3x3 micro-cell.
        const auto m = mesher.remesh(
            unit_quad({-1.0, 1.0, 1.0, 1.0}),
            3, RemeshMode::QuadDominant);
        assert(!m.cells.empty());
        assert(m.quad_count() > 0);
        check_area(m);
        assert(max_cell_area(m) <= (1.0 / 9.0) + 1.0e-12);
    }

    {
        const auto m = mesher.remesh(
            unit_quad({-1.0, 1.0, 1.0, 1.0}),
            3, RemeshMode::TriangleOnly);
        assert(m.quad_count() == 0);
        assert(m.triangle_count() > 0);
        check_area(m);
        assert(max_cell_area(m) <= (1.0 / 9.0) + 1.0e-12);
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
