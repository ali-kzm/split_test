#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace split_test {

struct Point2 {
    double x{};
    double y{};
};

struct Vertex {
    Point2 p{};
    double phi{};
};

enum class CellType {
    Triangle = 3,
    Quad = 4
};

enum class RemeshMode {
    QuadDominant,
    TriangleOnly
};

struct Cell {
    CellType type{CellType::Triangle};
    std::vector<std::size_t> nodes;
    int phase{}; // -1 for phi < 0, +1 for phi > 0
};

struct MixedMesh {
    std::vector<Vertex> nodes;
    std::vector<Cell> cells;

    [[nodiscard]] std::size_t triangle_count() const;
    [[nodiscard]] std::size_t quad_count() const;
    [[nodiscard]] double total_area() const;
};

struct QuadInput {
    // Nodes must be ordered around the quadrilateral boundary (CW or CCW).
    std::array<Point2, 4> nodes{};
    std::array<double, 4> phi{};
};

class QuadLevelSetMesher {
public:
    explicit QuadLevelSetMesher(double epsilon = 1.0e-12);

    // divisions = N subdivides each coarse quad into N x N quads and each
    // coarse triangle into N^2 triangles. N=1 preserves the coarse topology.
    [[nodiscard]] MixedMesh remesh(
        const QuadInput& input,
        std::size_t divisions = 1,
        RemeshMode mode = RemeshMode::QuadDominant) const;

private:
    double eps_;
};

void write_legacy_vtk(const std::string& filename, const MixedMesh& mesh);

} // namespace split_test
