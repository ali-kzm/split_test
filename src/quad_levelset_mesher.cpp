#include "quad_levelset_mesher.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace split_test {
namespace {

struct PolyVertex {
    Point2 p{};
    double phi{};
};

using Polygon = std::vector<PolyVertex>;

[[nodiscard]] double cross(const Point2& a, const Point2& b, const Point2& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

[[nodiscard]] double signed_area(const Polygon& p) {
    double a = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const auto& u = p[i].p;
        const auto& v = p[(i + 1) % p.size()].p;
        a += u.x * v.y - v.x * u.y;
    }
    return 0.5 * a;
}

[[nodiscard]] double polygon_area(const Polygon& p) {
    return std::abs(signed_area(p));
}

[[nodiscard]] bool nearly_equal(double a, double b, double eps) {
    return std::abs(a - b) <= eps;
}

[[nodiscard]] bool same_point(const Point2& a, const Point2& b, double eps) {
    return nearly_equal(a.x, b.x, eps) && nearly_equal(a.y, b.y, eps);
}

void remove_duplicate_vertices(Polygon& p, double eps) {
    if (p.empty()) return;

    Polygon cleaned;
    cleaned.reserve(p.size());
    for (const auto& v : p) {
        if (cleaned.empty() || !same_point(cleaned.back().p, v.p, eps)) {
            cleaned.push_back(v);
        }
    }
    if (cleaned.size() > 1 && same_point(cleaned.front().p, cleaned.back().p, eps)) {
        cleaned.pop_back();
    }
    p = std::move(cleaned);
}

[[nodiscard]] PolyVertex interpolate_zero(const PolyVertex& a, const PolyVertex& b) {
    const double denom = a.phi - b.phi;
    if (std::abs(denom) <= std::numeric_limits<double>::epsilon()) {
        return {{0.5 * (a.p.x + b.p.x), 0.5 * (a.p.y + b.p.y)}, 0.0};
    }
    const double t = a.phi / denom;
    return {{a.p.x + t * (b.p.x - a.p.x), a.p.y + t * (b.p.y - a.p.y)}, 0.0};
}

[[nodiscard]] Polygon clip_by_sign(const Polygon& in, int phase, double eps) {
    Polygon out;
    if (in.empty()) return out;

    auto inside = [phase, eps](double phi) {
        return phase > 0 ? phi >= -eps : phi <= eps;
    };

    PolyVertex a = in.back();
    bool a_inside = inside(a.phi);
    for (const auto& b : in) {
        const bool b_inside = inside(b.phi);
        if (a_inside && b_inside) {
            out.push_back(b);
        } else if (a_inside && !b_inside) {
            out.push_back(interpolate_zero(a, b));
        } else if (!a_inside && b_inside) {
            out.push_back(interpolate_zero(a, b));
            out.push_back(b);
        }
        a = b;
        a_inside = b_inside;
    }

    remove_duplicate_vertices(out, eps);
    if (out.size() >= 3 && signed_area(out) < 0.0) {
        std::reverse(out.begin(), out.end());
    }
    return out;
}

[[nodiscard]] double triangle_quality(const Point2& a, const Point2& b, const Point2& c) {
    const double area2 = std::abs(cross(a, b, c));
    const double ab2 = (a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y);
    const double bc2 = (b.x-c.x)*(b.x-c.x) + (b.y-c.y)*(b.y-c.y);
    const double ca2 = (c.x-a.x)*(c.x-a.x) + (c.y-a.y)*(c.y-a.y);
    const double denom = ab2 + bc2 + ca2;
    return denom > 0.0 ? (2.0 * std::sqrt(3.0) * area2 / denom) : 0.0;
}

[[nodiscard]] double quad_quality(const std::array<Point2, 4>& q) {
    double score = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 4; ++i) {
        const auto& prev = q[(i + 3) % 4];
        const auto& cur  = q[i];
        const auto& next = q[(i + 1) % 4];
        score = std::min(score, triangle_quality(prev, cur, next));
    }
    return std::isfinite(score) ? score : 0.0;
}

[[nodiscard]] std::size_t add_node(MixedMesh& mesh, const PolyVertex& v, double eps) {
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (same_point(mesh.nodes[i].p, v.p, eps)) {
            if (std::abs(v.phi) < std::abs(mesh.nodes[i].phi)) {
                mesh.nodes[i].phi = v.phi;
            }
            return i;
        }
    }
    mesh.nodes.push_back({v.p, v.phi});
    return mesh.nodes.size() - 1;
}

void add_cell(MixedMesh& mesh, const Polygon& poly, int phase, double eps) {
    if (poly.size() != 3 && poly.size() != 4) {
        throw std::runtime_error("Internal error: add_cell expects a triangle or quadrilateral");
    }
    if (polygon_area(poly) <= eps * eps) return;

    Cell cell;
    cell.type = poly.size() == 3 ? CellType::Triangle : CellType::Quad;
    cell.phase = phase;
    cell.nodes.reserve(poly.size());
    for (const auto& v : poly) {
        cell.nodes.push_back(add_node(mesh, v, eps));
    }
    mesh.cells.push_back(std::move(cell));
}

void add_quad_dominant_triangle(
    MixedMesh& mesh, Polygon tri, int phase, double eps) {
    if (tri.size() != 3) {
        throw std::runtime_error(
            "Internal error: quad-dominant triangle template expects 3 vertices");
    }
    if (polygon_area(tri) <= eps * eps) return;
    if (signed_area(tri) < 0.0) {
        std::reverse(tri.begin(), tri.end());
    }

    // Keep exactly one triangle in the interior and surround it with 3 quads.
    // No extra points are inserted on the outer boundary, so this template is
    // conforming with neighboring coarse/micro cells.
    constexpr double inner_scale = 0.35;

    PolyVertex center;
    center.p.x = (tri[0].p.x + tri[1].p.x + tri[2].p.x) / 3.0;
    center.p.y = (tri[0].p.y + tri[1].p.y + tri[2].p.y) / 3.0;
    center.phi = (tri[0].phi + tri[1].phi + tri[2].phi) / 3.0;

    std::array<PolyVertex, 3> inner{};
    for (std::size_t i = 0; i < 3; ++i) {
        inner[i].p.x =
            center.p.x + inner_scale * (tri[i].p.x - center.p.x);
        inner[i].p.y =
            center.p.y + inner_scale * (tri[i].p.y - center.p.y);
        inner[i].phi =
            center.phi + inner_scale * (tri[i].phi - center.phi);
    }

    add_cell(mesh, Polygon{tri[0], tri[1], inner[1], inner[0]}, phase, eps);
    add_cell(mesh, Polygon{tri[1], tri[2], inner[2], inner[1]}, phase, eps);
    add_cell(mesh, Polygon{tri[2], tri[0], inner[0], inner[2]}, phase, eps);
    add_cell(mesh, Polygon{inner[0], inner[1], inner[2]}, phase, eps);
}

void decompose_polygon(
    MixedMesh& mesh, Polygon poly, int phase, double eps,
    bool quad_dominant) {
    remove_duplicate_vertices(poly, eps);
    if (poly.size() < 3 || polygon_area(poly) <= eps * eps) return;
    if (signed_area(poly) < 0.0) std::reverse(poly.begin(), poly.end());

    if (poly.size() == 3) {
        if (quad_dominant) {
            add_quad_dominant_triangle(mesh, poly, phase, eps);
        } else {
            add_cell(mesh, poly, phase, eps);
        }
        return;
    }

    if (poly.size() == 4) {
        add_cell(mesh, poly, phase, eps);
        return;
    }

    if (poly.size() == 5) {
        double best_score = -1.0;
        std::size_t best_i = 0;
        for (std::size_t i = 0; i < 5; ++i) {
            const std::size_t i0 = i;
            const std::size_t i1 = (i + 1) % 5;
            const std::size_t i2 = (i + 2) % 5;
            const std::size_t i3 = (i + 3) % 5;
            const std::size_t i4 = (i + 4) % 5;

            const double tq = triangle_quality(poly[i0].p, poly[i1].p, poly[i2].p);
            const std::array<Point2, 4> q{poly[i0].p, poly[i2].p, poly[i3].p, poly[i4].p};
            const double qq = quad_quality(q);
            const double score = std::min(tq, qq);
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }

        const std::size_t i0 = best_i;
        const std::size_t i1 = (best_i + 1) % 5;
        const std::size_t i2 = (best_i + 2) % 5;
        const std::size_t i3 = (best_i + 3) % 5;
        const std::size_t i4 = (best_i + 4) % 5;
        const Polygon triangle{poly[i0], poly[i1], poly[i2]};
        if (quad_dominant) {
            add_quad_dominant_triangle(mesh, triangle, phase, eps);
        } else {
            add_cell(mesh, triangle, phase, eps);
        }
        add_cell(mesh, Polygon{poly[i0], poly[i2], poly[i3], poly[i4]}, phase, eps);
        return;
    }

    for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
        const Polygon triangle{poly[0], poly[i], poly[i + 1]};
        if (quad_dominant) {
            add_quad_dominant_triangle(mesh, triangle, phase, eps);
        } else {
            add_cell(mesh, triangle, phase, eps);
        }
    }
}

[[nodiscard]] bool is_convex_quad(const Polygon& q, double eps) {
    if (q.size() != 4) return false;
    int sign = 0;
    for (int i = 0; i < 4; ++i) {
        const double c = cross(q[i].p, q[(i + 1) % 4].p, q[(i + 2) % 4].p);
        if (std::abs(c) <= eps) return false;
        const int s = c > 0.0 ? 1 : -1;
        if (sign == 0) sign = s;
        else if (s != sign) return false;
    }
    return true;
}

[[nodiscard]] int strict_sign(double v, double eps) {
    if (v > eps) return 1;
    if (v < -eps) return -1;
    return 0;
}

void process_triangle(
    MixedMesh& mesh, const Polygon& tri, double eps,
    bool quad_dominant) {
    const Polygon pos = clip_by_sign(tri, +1, eps);
    const Polygon neg = clip_by_sign(tri, -1, eps);
    decompose_polygon(mesh, pos, +1, eps, quad_dominant);
    decompose_polygon(mesh, neg, -1, eps, quad_dominant);
}

[[nodiscard]] PolyVertex bilinear_vertex(
    const QuadInput& input, double u, double v) {
    const double w0 = (1.0-u)*(1.0-v);
    const double w1 = u*(1.0-v);
    const double w2 = u*v;
    const double w3 = (1.0-u)*v;

    PolyVertex out;
    out.p.x = w0*input.nodes[0].x + w1*input.nodes[1].x
            + w2*input.nodes[2].x + w3*input.nodes[3].x;
    out.p.y = w0*input.nodes[0].y + w1*input.nodes[1].y
            + w2*input.nodes[2].y + w3*input.nodes[3].y;
    out.phi = w0*input.phi[0] + w1*input.phi[1]
            + w2*input.phi[2] + w3*input.phi[3];
    return out;
}

[[nodiscard]] MixedMesh triangulate_mesh(const MixedMesh& input, double eps) {
    MixedMesh out;
    for (const auto& cell : input.cells) {
        if (cell.type == CellType::Triangle) {
            Polygon tri;
            for (const auto id : cell.nodes) {
                tri.push_back({input.nodes[id].p, input.nodes[id].phi});
            }
            add_cell(out, tri, cell.phase, eps);
            continue;
        }

        const PolyVertex q0{input.nodes[cell.nodes[0]].p, input.nodes[cell.nodes[0]].phi};
        const PolyVertex q1{input.nodes[cell.nodes[1]].p, input.nodes[cell.nodes[1]].phi};
        const PolyVertex q2{input.nodes[cell.nodes[2]].p, input.nodes[cell.nodes[2]].phi};
        const PolyVertex q3{input.nodes[cell.nodes[3]].p, input.nodes[cell.nodes[3]].phi};

        const double score02 = std::min(
            triangle_quality(q0.p, q1.p, q2.p),
            triangle_quality(q0.p, q2.p, q3.p));
        const double score13 = std::min(
            triangle_quality(q1.p, q2.p, q3.p),
            triangle_quality(q1.p, q3.p, q0.p));

        if (score02 >= score13) {
            add_cell(out, Polygon{q0, q1, q2}, cell.phase, eps);
            add_cell(out, Polygon{q0, q2, q3}, cell.phase, eps);
        } else {
            add_cell(out, Polygon{q1, q2, q3}, cell.phase, eps);
            add_cell(out, Polygon{q1, q3, q0}, cell.phase, eps);
        }
    }
    return out;
}

void append_mesh(MixedMesh& destination, const MixedMesh& source, double eps) {
    for (const auto& cell : source.cells) {
        Polygon poly;
        poly.reserve(cell.nodes.size());
        for (const auto id : cell.nodes) {
            poly.push_back({source.nodes[id].p, source.nodes[id].phi});
        }
        add_cell(destination, poly, cell.phase, eps);
    }
}

[[nodiscard]] MixedMesh remesh_coarse_quad(
    const QuadInput& input, double eps, bool quad_dominant) {
    Polygon q;
    q.reserve(4);
    for (std::size_t i = 0; i < 4; ++i) {
        q.push_back({input.nodes[i], input.phi[i]});
    }

    if (!is_convex_quad(q, eps)) {
        throw std::invalid_argument(
            "Input nodes must form a non-degenerate convex quadrilateral in boundary order");
    }
    if (signed_area(q) < 0.0) {
        std::reverse(q.begin(), q.end());
    }

    int npos = 0;
    int nneg = 0;
    for (const auto& v : q) {
        const int s = strict_sign(v.phi, eps);
        if (s > 0) ++npos;
        if (s < 0) ++nneg;
    }

    MixedMesh mesh;
    if (nneg == 0) {
        decompose_polygon(mesh, q, +1, eps, quad_dominant);
        return mesh;
    }
    if (npos == 0) {
        decompose_polygon(mesh, q, -1, eps, quad_dominant);
        return mesh;
    }

    const std::array<int, 4> s{
        strict_sign(q[0].phi, eps), strict_sign(q[1].phi, eps),
        strict_sign(q[2].phi, eps), strict_sign(q[3].phi, eps)};
    const bool alternating =
        s[0] != 0 && s[1] != 0 && s[2] != 0 && s[3] != 0 &&
        s[0] == s[2] && s[1] == s[3] && s[0] != s[1];

    if (alternating) {
        const double center_phi =
            0.25 * (q[0].phi + q[1].phi + q[2].phi + q[3].phi);
        const int raw_center_sign = strict_sign(center_phi, eps);
        const int center_sign = raw_center_sign != 0 ? raw_center_sign : s[0];

        if (s[0] == center_sign) {
            process_triangle(mesh, Polygon{q[0], q[1], q[2]}, eps, quad_dominant);
            process_triangle(mesh, Polygon{q[0], q[2], q[3]}, eps, quad_dominant);
        } else {
            process_triangle(mesh, Polygon{q[1], q[2], q[3]}, eps, quad_dominant);
            process_triangle(mesh, Polygon{q[1], q[3], q[0]}, eps, quad_dominant);
        }
        return mesh;
    }

    const Polygon pos = clip_by_sign(q, +1, eps);
    const Polygon neg = clip_by_sign(q, -1, eps);
    decompose_polygon(mesh, pos, +1, eps, quad_dominant);
    decompose_polygon(mesh, neg, -1, eps, quad_dominant);
    return mesh;
}

} // namespace

std::size_t MixedMesh::triangle_count() const {
    return static_cast<std::size_t>(std::count_if(cells.begin(), cells.end(), [](const Cell& c) {
        return c.type == CellType::Triangle;
    }));
}

std::size_t MixedMesh::quad_count() const {
    return static_cast<std::size_t>(std::count_if(cells.begin(), cells.end(), [](const Cell& c) {
        return c.type == CellType::Quad;
    }));
}

double MixedMesh::total_area() const {
    double area = 0.0;
    for (const auto& c : cells) {
        Polygon p;
        p.reserve(c.nodes.size());
        for (const auto id : c.nodes) {
            p.push_back({nodes.at(id).p, nodes.at(id).phi});
        }
        area += polygon_area(p);
    }
    return area;
}

QuadLevelSetMesher::QuadLevelSetMesher(double epsilon) : eps_(epsilon) {
    if (!(eps_ > 0.0)) {
        throw std::invalid_argument("epsilon must be positive");
    }
}

MixedMesh QuadLevelSetMesher::remesh(
    const QuadInput& input, std::size_t divisions, RemeshMode mode) const {
    if (divisions == 0) {
        throw std::invalid_argument("divisions must be >= 1");
    }

    // Validate the parent Q1 quadrilateral once.
    Polygon parent;
    parent.reserve(4);
    for (std::size_t i = 0; i < 4; ++i) {
        parent.push_back({input.nodes[i], input.phi[i]});
    }
    if (!is_convex_quad(parent, eps_)) {
        throw std::invalid_argument(
            "Input nodes must form a non-degenerate convex quadrilateral in boundary order");
    }

    // Important: refine the parent element FIRST.  Each micro-Q1 cell is then
    // cut by its locally linearized phi=0 interface.  This avoids the old
    // behavior where tiny elements clustered in one phase polygon while long
    // coarse cells survived elsewhere.
    MixedMesh mesh;
    for (std::size_t j = 0; j < divisions; ++j) {
        const double v0 = static_cast<double>(j) /
                          static_cast<double>(divisions);
        const double v1 = static_cast<double>(j + 1) /
                          static_cast<double>(divisions);

        for (std::size_t i = 0; i < divisions; ++i) {
            const double u0 = static_cast<double>(i) /
                              static_cast<double>(divisions);
            const double u1 = static_cast<double>(i + 1) /
                              static_cast<double>(divisions);

            const PolyVertex p00 = bilinear_vertex(input, u0, v0);
            const PolyVertex p10 = bilinear_vertex(input, u1, v0);
            const PolyVertex p11 = bilinear_vertex(input, u1, v1);
            const PolyVertex p01 = bilinear_vertex(input, u0, v1);

            QuadInput sub;
            sub.nodes = {p00.p, p10.p, p11.p, p01.p};
            sub.phi   = {p00.phi, p10.phi, p11.phi, p01.phi};

            const MixedMesh local = remesh_coarse_quad(
                sub, eps_, mode == RemeshMode::QuadDominant);
            append_mesh(mesh, local, eps_);
        }
    }

    if (mode == RemeshMode::TriangleOnly) {
        mesh = triangulate_mesh(mesh, eps_);
    }
    return mesh;
}

void write_legacy_vtk(const std::string& filename, const MixedMesh& mesh) {
    std::ofstream out(filename);
    if (!out) throw std::runtime_error("Cannot open VTK file: " + filename);

    out << "# vtk DataFile Version 3.0\n";
    out << "Q1 level-set split mesh\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << std::setprecision(17);

    out << "POINTS " << mesh.nodes.size() << " double\n";
    for (const auto& n : mesh.nodes) {
        out << n.p.x << ' ' << n.p.y << " 0\n";
    }

    std::size_t connectivity_size = 0;
    for (const auto& c : mesh.cells) connectivity_size += 1 + c.nodes.size();
    out << "CELLS " << mesh.cells.size() << ' ' << connectivity_size << "\n";
    for (const auto& c : mesh.cells) {
        out << c.nodes.size();
        for (const auto id : c.nodes) out << ' ' << id;
        out << '\n';
    }

    out << "CELL_TYPES " << mesh.cells.size() << "\n";
    for (const auto& c : mesh.cells) {
        out << (c.type == CellType::Triangle ? 5 : 9) << '\n';
    }

    out << "POINT_DATA " << mesh.nodes.size() << "\n";
    out << "SCALARS phi double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& n : mesh.nodes) out << n.phi << '\n';

    out << "CELL_DATA " << mesh.cells.size() << "\n";
    out << "SCALARS phase int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) out << c.phase << '\n';

    out << "SCALARS element_type int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& c : mesh.cells) {
        out << (c.type == CellType::Triangle ? 3 : 4) << '\n';
    }
}

} // namespace split_test
