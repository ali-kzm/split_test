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

void decompose_polygon(MixedMesh& mesh, Polygon poly, int phase, double eps) {
    remove_duplicate_vertices(poly, eps);
    if (poly.size() < 3 || polygon_area(poly) <= eps * eps) return;
    if (signed_area(poly) < 0.0) std::reverse(poly.begin(), poly.end());

    if (poly.size() == 3 || poly.size() == 4) {
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
        add_cell(mesh, Polygon{poly[i0], poly[i1], poly[i2]}, phase, eps);
        add_cell(mesh, Polygon{poly[i0], poly[i2], poly[i3], poly[i4]}, phase, eps);
        return;
    }

    for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
        add_cell(mesh, Polygon{poly[0], poly[i], poly[i + 1]}, phase, eps);
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

void process_triangle(MixedMesh& mesh, const Polygon& tri, double eps) {
    const Polygon pos = clip_by_sign(tri, +1, eps);
    const Polygon neg = clip_by_sign(tri, -1, eps);
    decompose_polygon(mesh, pos, +1, eps);
    decompose_polygon(mesh, neg, -1, eps);
}

[[nodiscard]] PolyVertex lerp_triangle(
    const PolyVertex& a, const PolyVertex& b, const PolyVertex& c,
    std::size_t i, std::size_t j, std::size_t n) {
    const double wb = static_cast<double>(i) / static_cast<double>(n);
    const double wc = static_cast<double>(j) / static_cast<double>(n);
    const double wa = 1.0 - wb - wc;
    return {{wa*a.p.x + wb*b.p.x + wc*c.p.x,
             wa*a.p.y + wb*b.p.y + wc*c.p.y},
            wa*a.phi + wb*b.phi + wc*c.phi};
}

[[nodiscard]] PolyVertex lerp_quad(
    const PolyVertex& q0, const PolyVertex& q1,
    const PolyVertex& q2, const PolyVertex& q3,
    std::size_t i, std::size_t j, std::size_t n) {
    const double u = static_cast<double>(i) / static_cast<double>(n);
    const double v = static_cast<double>(j) / static_cast<double>(n);
    const double w0 = (1.0-u)*(1.0-v);
    const double w1 = u*(1.0-v);
    const double w2 = u*v;
    const double w3 = (1.0-u)*v;
    return {{w0*q0.p.x + w1*q1.p.x + w2*q2.p.x + w3*q3.p.x,
             w0*q0.p.y + w1*q1.p.y + w2*q2.p.y + w3*q3.p.y},
            w0*q0.phi + w1*q1.phi + w2*q2.phi + w3*q3.phi};
}

[[nodiscard]] MixedMesh triangulate_mesh(const MixedMesh& input, double eps) {
    MixedMesh out;
    for (const auto& cell : input.cells) {
        if (cell.type == CellType::Triangle) {
            Polygon tri;
            for (const auto id : cell.nodes) tri.push_back({input.nodes[id].p, input.nodes[id].phi});
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

[[nodiscard]] MixedMesh subdivide_mesh(const MixedMesh& input, std::size_t divisions, double eps) {
    if (divisions == 1) return input;

    MixedMesh out;
    for (const auto& cell : input.cells) {
        if (cell.type == CellType::Triangle) {
            const PolyVertex a{input.nodes[cell.nodes[0]].p, input.nodes[cell.nodes[0]].phi};
            const PolyVertex b{input.nodes[cell.nodes[1]].p, input.nodes[cell.nodes[1]].phi};
            const PolyVertex c{input.nodes[cell.nodes[2]].p, input.nodes[cell.nodes[2]].phi};

            for (std::size_t j = 0; j < divisions; ++j) {
                for (std::size_t i = 0; i + j < divisions; ++i) {
                    const auto p00 = lerp_triangle(a, b, c, i, j, divisions);
                    const auto p10 = lerp_triangle(a, b, c, i + 1, j, divisions);
                    const auto p01 = lerp_triangle(a, b, c, i, j + 1, divisions);
                    add_cell(out, Polygon{p00, p10, p01}, cell.phase, eps);

                    if (i + j + 1 < divisions) {
                        const auto p11 = lerp_triangle(a, b, c, i + 1, j + 1, divisions);
                        add_cell(out, Polygon{p10, p11, p01}, cell.phase, eps);
                    }
                }
            }
        } else {
            const PolyVertex q0{input.nodes[cell.nodes[0]].p, input.nodes[cell.nodes[0]].phi};
            const PolyVertex q1{input.nodes[cell.nodes[1]].p, input.nodes[cell.nodes[1]].phi};
            const PolyVertex q2{input.nodes[cell.nodes[2]].p, input.nodes[cell.nodes[2]].phi};
            const PolyVertex q3{input.nodes[cell.nodes[3]].p, input.nodes[cell.nodes[3]].phi};

            for (std::size_t j = 0; j < divisions; ++j) {
                for (std::size_t i = 0; i < divisions; ++i) {
                    const auto p00 = lerp_quad(q0, q1, q2, q3, i, j, divisions);
                    const auto p10 = lerp_quad(q0, q1, q2, q3, i + 1, j, divisions);
                    const auto p11 = lerp_quad(q0, q1, q2, q3, i + 1, j + 1, divisions);
                    const auto p01 = lerp_quad(q0, q1, q2, q3, i, j + 1, divisions);
                    add_cell(out, Polygon{p00, p10, p11, p01}, cell.phase, eps);
                }
            }
        }
    }
    return out;
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

    Polygon q;
    q.reserve(4);
    for (std::size_t i = 0; i < 4; ++i) q.push_back({input.nodes[i], input.phi[i]});

    if (!is_convex_quad(q, eps_)) {
        throw std::invalid_argument("Input nodes must form a non-degenerate convex quadrilateral in boundary order");
    }
    if (signed_area(q) < 0.0) std::reverse(q.begin(), q.end());

    auto finish = [&](MixedMesh mesh) {
        if (mode == RemeshMode::TriangleOnly) {
            mesh = triangulate_mesh(mesh, eps_);
        }
        return subdivide_mesh(mesh, divisions, eps_);
    };

    int npos = 0;
    int nneg = 0;
    for (const auto& v : q) {
        const int s = strict_sign(v.phi, eps_);
        if (s > 0) ++npos;
        if (s < 0) ++nneg;
    }

    MixedMesh mesh;
    if (nneg == 0) {
        decompose_polygon(mesh, q, +1, eps_);
        return finish(std::move(mesh));
    }
    if (npos == 0) {
        decompose_polygon(mesh, q, -1, eps_);
        return finish(std::move(mesh));
    }

    const std::array<int, 4> s{
        strict_sign(q[0].phi, eps_), strict_sign(q[1].phi, eps_),
        strict_sign(q[2].phi, eps_), strict_sign(q[3].phi, eps_)};
    const bool alternating = s[0] != 0 && s[1] != 0 && s[2] != 0 && s[3] != 0 &&
                             s[0] == s[2] && s[1] == s[3] && s[0] != s[1];

    if (alternating) {
        const double center_phi = 0.25 * (q[0].phi + q[1].phi + q[2].phi + q[3].phi);
        const int raw_center_sign = strict_sign(center_phi, eps_);
        const int center_sign = raw_center_sign != 0 ? raw_center_sign : s[0];

        if (s[0] == center_sign) {
            process_triangle(mesh, Polygon{q[0], q[1], q[2]}, eps_);
            process_triangle(mesh, Polygon{q[0], q[2], q[3]}, eps_);
        } else {
            process_triangle(mesh, Polygon{q[1], q[2], q[3]}, eps_);
            process_triangle(mesh, Polygon{q[1], q[3], q[0]}, eps_);
        }
        return finish(std::move(mesh));
    }

    const Polygon pos = clip_by_sign(q, +1, eps_);
    const Polygon neg = clip_by_sign(q, -1, eps_);
    decompose_polygon(mesh, pos, +1, eps_);
    decompose_polygon(mesh, neg, -1, eps_);
    return finish(std::move(mesh));
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
