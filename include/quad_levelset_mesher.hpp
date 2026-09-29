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

    // Split one convex Q1 quadrilateral using linear interpolation of phi on edges.
    // Typical 2-edge cuts produce a quad-dominant triangle/quad mesh.
    // The alternating-sign Q1 saddle case is resolved by a center-sign diagonal.
    [[nodiscard]] MixedMesh remesh(const QuadInput& input) const;

private:
    double eps_;
};

void write_legacy_vtk(const std::string& filename, const MixedMesh& mesh);

} // namespace split_test
