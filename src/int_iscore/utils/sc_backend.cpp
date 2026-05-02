#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

struct Atom {
    double x;
    double y;
    double z;
    double radius;
    int molecule;
};

struct SurfacePoint {
    std::array<double, 3> pos{};
    std::array<double, 3> normal{};
    std::size_t atom_index{0};
    int molecule{0};
};

struct SurfaceStats {
    std::size_t n_surface_m1{0};
    std::size_t n_surface_m2{0};
    std::size_t n_buried_m1{0};
    std::size_t n_buried_m2{0};
    std::size_t n_trimmed_m1{0};
    std::size_t n_trimmed_m2{0};
};

struct CellKey {
    int x;
    int y;
    int z;

    bool operator==(const CellKey& other) const noexcept {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& key) const noexcept {
        std::size_t h1 = std::hash<int>{}(key.x * 73856093);
        std::size_t h2 = std::hash<int>{}(key.y * 19349663);
        std::size_t h3 = std::hash<int>{}(key.z * 83492791);
        return h1 ^ h2 ^ h3;
    }
};

struct GridSpec {
    std::array<double, 3> origin{};
    std::array<int, 3> dims{};
    double spacing{0.25};
};

double squared_distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

double squared_distance_xyz(const Atom& atom, const std::array<double, 3>& p) {
    const double dx = atom.x - p[0];
    const double dy = atom.y - p[1];
    const double dz = atom.z - p[2];
    return dx * dx + dy * dy + dz * dz;
}

std::array<double, 3> normalize(std::array<double, 3> v) {
    const double norm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm < 1e-12) {
        return {0.0, 0.0, 1.0};
    }
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
    return v;
}

class SurfaceGridIndex {
public:
    SurfaceGridIndex(const std::vector<SurfacePoint>& points, double cell_size)
        : points_(points), cell_size_(cell_size) {
        buckets_.reserve(points.size());
        for (std::size_t i = 0; i < points.size(); ++i) {
            buckets_[cell_key(points[i].pos)].push_back(i);
        }
    }

    std::pair<std::size_t, double> nearest(const SurfacePoint& query, double max_distance) const {
        if (points_.empty()) {
            return {0, std::numeric_limits<double>::infinity()};
        }

        const CellKey qcell = cell_key(query.pos);
        const int max_shell = static_cast<int>(std::ceil(max_distance / cell_size_));
        double best_d2 = max_distance * max_distance;
        std::size_t best_idx = 0;
        bool found = false;

        for (int shell = 0; shell <= max_shell; ++shell) {
            const double shell_lb = shell == 0 ? 0.0 : (shell - 1) * cell_size_;
            if (found && shell_lb * shell_lb >= best_d2) {
                break;
            }

            for (int dx = -shell; dx <= shell; ++dx) {
                for (int dy = -shell; dy <= shell; ++dy) {
                    for (int dz = -shell; dz <= shell; ++dz) {
                        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != shell) {
                            continue;
                        }

                        CellKey key{qcell.x + dx, qcell.y + dy, qcell.z + dz};
                        const double cell_lb2 = min_distance_sq_to_cell(query.pos, key);
                        if (cell_lb2 >= best_d2) {
                            continue;
                        }

                        auto it = buckets_.find(key);
                        if (it == buckets_.end()) {
                            continue;
                        }

                        for (std::size_t idx : it->second) {
                            const double d2 = squared_distance(query.pos, points_[idx].pos);
                            if (d2 < best_d2) {
                                best_d2 = d2;
                                best_idx = idx;
                                found = true;
                            }
                        }
                    }
                }
            }
        }

        if (!found) {
            for (std::size_t i = 0; i < points_.size(); ++i) {
                const double d2 = squared_distance(query.pos, points_[i].pos);
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_idx = i;
                    found = true;
                }
            }
        }

        return {best_idx, std::sqrt(best_d2)};
    }

private:
    CellKey cell_key(const std::array<double, 3>& pos) const {
        return CellKey{
            static_cast<int>(std::floor(pos[0] / cell_size_)),
            static_cast<int>(std::floor(pos[1] / cell_size_)),
            static_cast<int>(std::floor(pos[2] / cell_size_)),
        };
    }

    double min_distance_sq_to_cell(const std::array<double, 3>& pos, const CellKey& key) const {
        const double min_x = key.x * cell_size_;
        const double min_y = key.y * cell_size_;
        const double min_z = key.z * cell_size_;
        const double max_x = min_x + cell_size_;
        const double max_y = min_y + cell_size_;
        const double max_z = min_z + cell_size_;

        const double dx = pos[0] < min_x ? (min_x - pos[0]) : (pos[0] > max_x ? pos[0] - max_x : 0.0);
        const double dy = pos[1] < min_y ? (min_y - pos[1]) : (pos[1] > max_y ? pos[1] - max_y : 0.0);
        const double dz = pos[2] < min_z ? (min_z - pos[2]) : (pos[2] > max_z ? pos[2] - max_z : 0.0);
        return dx * dx + dy * dy + dz * dz;
    }

    const std::vector<SurfacePoint>& points_;
    double cell_size_;
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> buckets_;
};

std::vector<Atom> load_atoms(
    py::array_t<double, py::array::c_style | py::array::forcecast> coords,
    py::array_t<double, py::array::c_style | py::array::forcecast> radii,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> molecule_ids) {

    if (coords.ndim() != 2 || coords.shape(1) != 3) {
        throw std::runtime_error("coords must be an Nx3 float64 array");
    }
    if (radii.ndim() != 1 || molecule_ids.ndim() != 1) {
        throw std::runtime_error("radii and molecule_ids must be 1D arrays");
    }
    if (coords.shape(0) != radii.shape(0) || coords.shape(0) != molecule_ids.shape(0)) {
        throw std::runtime_error("coords, radii, and molecule_ids must have the same length");
    }

    const auto c = coords.unchecked<2>();
    const auto r = radii.unchecked<1>();
    const auto m = molecule_ids.unchecked<1>();

    std::vector<Atom> atoms;
    atoms.reserve(static_cast<std::size_t>(coords.shape(0)));
    for (ssize_t i = 0; i < coords.shape(0); ++i) {
        atoms.push_back(Atom{c(i, 0), c(i, 1), c(i, 2), r(i), m(i)});
    }
    return atoms;
}

GridSpec make_grid(const std::vector<Atom>& atoms, double probe_radius, double dot_density) {
    if (atoms.empty()) {
        throw std::runtime_error("No atoms supplied");
    }

    double min_x = std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    for (const auto& atom : atoms) {
        const double extent = atom.radius + probe_radius + 2.0;
        min_x = std::min(min_x, atom.x - extent);
        min_y = std::min(min_y, atom.y - extent);
        min_z = std::min(min_z, atom.z - extent);
        max_x = std::max(max_x, atom.x + extent);
        max_y = std::max(max_y, atom.y + extent);
        max_z = std::max(max_z, atom.z + extent);
    }

    // Tie grid resolution to requested dot density. This is still an approximation,
    // but it keeps the backend parameterization compatible with the SC API.
    double spacing = 1.0 / std::sqrt(std::max(1.0, dot_density));
    spacing = std::clamp(spacing, 0.20, 0.50);

    GridSpec grid;
    grid.origin = {min_x, min_y, min_z};
    grid.spacing = spacing;
    grid.dims = {
        static_cast<int>(std::ceil((max_x - min_x) / spacing)) + 3,
        static_cast<int>(std::ceil((max_y - min_y) / spacing)) + 3,
        static_cast<int>(std::ceil((max_z - min_z) / spacing)) + 3,
    };
    return grid;
}

inline std::size_t flatten(const GridSpec& grid, int ix, int iy, int iz) {
    return static_cast<std::size_t>(ix) * static_cast<std::size_t>(grid.dims[1]) * static_cast<std::size_t>(grid.dims[2]) +
           static_cast<std::size_t>(iy) * static_cast<std::size_t>(grid.dims[2]) +
           static_cast<std::size_t>(iz);
}

std::vector<uint8_t> build_vdw_volume(const std::vector<Atom>& atoms, const GridSpec& grid) {
    const std::size_t total = static_cast<std::size_t>(grid.dims[0]) * grid.dims[1] * grid.dims[2];
    std::vector<uint8_t> volume(total, 0);

    for (const auto& atom : atoms) {
        const int rv = static_cast<int>(std::ceil(atom.radius / grid.spacing));
        const int cx = static_cast<int>(std::floor((atom.x - grid.origin[0]) / grid.spacing));
        const int cy = static_cast<int>(std::floor((atom.y - grid.origin[1]) / grid.spacing));
        const int cz = static_cast<int>(std::floor((atom.z - grid.origin[2]) / grid.spacing));

        for (int ix = std::max(0, cx - rv - 1); ix <= std::min(grid.dims[0] - 1, cx + rv + 1); ++ix) {
            for (int iy = std::max(0, cy - rv - 1); iy <= std::min(grid.dims[1] - 1, cy + rv + 1); ++iy) {
                for (int iz = std::max(0, cz - rv - 1); iz <= std::min(grid.dims[2] - 1, cz + rv + 1); ++iz) {
                    const std::array<double, 3> p{
                        grid.origin[0] + ix * grid.spacing,
                        grid.origin[1] + iy * grid.spacing,
                        grid.origin[2] + iz * grid.spacing,
                    };
                    if (squared_distance_xyz(atom, p) <= atom.radius * atom.radius) {
                        volume[flatten(grid, ix, iy, iz)] = 1;
                    }
                }
            }
        }
    }
    return volume;
}

std::vector<uint8_t> dilate_spherical(const std::vector<uint8_t>& input, const GridSpec& grid, double radius) {
    const int rv = std::max(1, static_cast<int>(std::ceil(radius / grid.spacing)));
    std::vector<std::array<int, 3>> offsets;
    offsets.reserve(static_cast<std::size_t>((2 * rv + 1) * (2 * rv + 1) * (2 * rv + 1)));
    for (int dx = -rv; dx <= rv; ++dx) {
        for (int dy = -rv; dy <= rv; ++dy) {
            for (int dz = -rv; dz <= rv; ++dz) {
                const double dist2 = static_cast<double>(dx * dx + dy * dy + dz * dz) * grid.spacing * grid.spacing;
                if (dist2 <= radius * radius) {
                    offsets.push_back({dx, dy, dz});
                }
            }
        }
    }

    std::vector<uint8_t> out(input.size(), 0);
    for (int ix = 0; ix < grid.dims[0]; ++ix) {
        for (int iy = 0; iy < grid.dims[1]; ++iy) {
            for (int iz = 0; iz < grid.dims[2]; ++iz) {
                if (!input[flatten(grid, ix, iy, iz)]) {
                    continue;
                }
                for (const auto& off : offsets) {
                    const int nx = ix + off[0];
                    const int ny = iy + off[1];
                    const int nz = iz + off[2];
                    if (nx >= 0 && nx < grid.dims[0] && ny >= 0 && ny < grid.dims[1] && nz >= 0 && nz < grid.dims[2]) {
                        out[flatten(grid, nx, ny, nz)] = 1;
                    }
                }
            }
        }
    }
    return out;
}

std::vector<SurfacePoint> extract_boundary_surface(
    const std::vector<uint8_t>& solid,
    const GridSpec& grid,
    const std::vector<Atom>& atoms,
    double probe_radius) {

    std::vector<SurfacePoint> points;
    points.reserve(solid.size() / 10);

    const std::array<std::array<int, 3>, 6> neigh{{
        {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}
    }};

    for (int ix = 1; ix < grid.dims[0] - 1; ++ix) {
        for (int iy = 1; iy < grid.dims[1] - 1; ++iy) {
            for (int iz = 1; iz < grid.dims[2] - 1; ++iz) {
                const auto idx = flatten(grid, ix, iy, iz);
                if (!solid[idx]) {
                    continue;
                }
                bool boundary = false;
                std::array<double, 3> grad{0.0, 0.0, 0.0};
                for (const auto& n : neigh) {
                    const int nx = ix + n[0];
                    const int ny = iy + n[1];
                    const int nz = iz + n[2];
                    const auto nidx = flatten(grid, nx, ny, nz);
                    if (!solid[nidx]) {
                        boundary = true;
                    }
                    grad[0] += static_cast<double>(solid[flatten(grid, ix + 1, iy, iz)] - solid[flatten(grid, ix - 1, iy, iz)]);
                    grad[1] += static_cast<double>(solid[flatten(grid, ix, iy + 1, iz)] - solid[flatten(grid, ix, iy - 1, iz)]);
                    grad[2] += static_cast<double>(solid[flatten(grid, ix, iy, iz + 1)] - solid[flatten(grid, ix, iy, iz - 1)]);
                }
                if (!boundary) {
                    continue;
                }

                std::array<double, 3> pos{
                    grid.origin[0] + ix * grid.spacing,
                    grid.origin[1] + iy * grid.spacing,
                    grid.origin[2] + iz * grid.spacing,
                };
                std::size_t nearest_atom = 0;
                double best = std::numeric_limits<double>::infinity();
                for (std::size_t a = 0; a < atoms.size(); ++a) {
                    const double d2 = squared_distance_xyz(atoms[a], pos);
                    if (d2 < best) {
                        best = d2;
                        nearest_atom = a;
                    }
                }
                points.push_back(SurfacePoint{pos, normalize({-grad[0], -grad[1], -grad[2]}), nearest_atom, atoms[nearest_atom].molecule});
            }
        }
    }
    return points;
}

double compute_sc_score(
    const std::vector<SurfacePoint>& m1,
    const std::vector<SurfacePoint>& m2,
    const std::vector<SurfacePoint>& combined,
    double weight,
    double trim,
    SurfaceStats& stats) {

    stats.n_surface_m1 = m1.size();
    stats.n_surface_m2 = m2.size();

    auto nearest_surface = [](const std::vector<SurfacePoint>& src, const std::vector<SurfacePoint>& other) {
        std::vector<double> distances(src.size(), std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i < src.size(); ++i) {
            for (std::size_t j = 0; j < other.size(); ++j) {
                const double d2 = squared_distance(src[i].pos, other[j].pos);
                if (d2 < distances[i]) {
                    distances[i] = d2;
                }
            }
            distances[i] = std::sqrt(distances[i]);
        }
        return distances;
    };

    const auto nearest1c = nearest_surface(m1, combined);
    const auto nearest2c = nearest_surface(m2, combined);

    auto classify = [&](const std::vector<SurfacePoint>& src,
                        const std::vector<double>& nearest_distances_to_complex,
                        std::size_t& buried_count,
                        std::size_t& trimmed_count) {
        std::vector<const SurfacePoint*> buried;
        std::vector<const SurfacePoint*> accessible;
        buried.reserve(src.size());
        accessible.reserve(src.size());

        for (std::size_t i = 0; i < src.size(); ++i) {
            // Points that disappear from the isolated molecular surface when the
            // complex surface is built are treated as buried. The threshold is tied
            // to the implicit dot spacing used by the backend.
            if (nearest_distances_to_complex[i] > 0.35) {
                buried.push_back(&src[i]);
            } else {
                accessible.push_back(&src[i]);
            }
        }

        buried_count = buried.size();

        std::vector<const SurfacePoint*> trimmed;
        trimmed.reserve(buried.size());
        for (const auto* point : buried) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto* acc : accessible) {
                best = std::min(best, std::sqrt(squared_distance(point->pos, acc->pos)));
            }
            if (best >= trim) {
                trimmed.push_back(point);
            }
        }
        trimmed_count = trimmed.size();
        return trimmed;
    };

    auto trimmed1 = classify(m1, nearest1c, stats.n_buried_m1, stats.n_trimmed_m1);
    auto trimmed2 = classify(m2, nearest2c, stats.n_buried_m2, stats.n_trimmed_m2);

    if (trimmed1.empty() || trimmed2.empty()) {
        return 0.0;
    }

    auto nearest_score = [&](const std::vector<const SurfacePoint*>& src, const std::vector<SurfacePoint>& target) {
        std::vector<double> values;
        values.reserve(src.size());
        SurfaceGridIndex index(target, 2.5);
        for (const auto* point : src) {
            auto [best_idx, best_dist] = index.nearest(*point, 10.0);
            const SurfacePoint* best_point = &target[best_idx];
            const double best = best_dist * best_dist;
            const double dot = -(point->normal[0] * best_point->normal[0] + point->normal[1] * best_point->normal[1] + point->normal[2] * best_point->normal[2]);
            values.push_back(dot * std::exp(-weight * best));
        }
        std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
        return values[values.size() / 2];
    };

    const double s12 = nearest_score(trimmed1, m2);
    const double s21 = nearest_score(trimmed2, m1);
    return 0.5 * (s12 + s21);
}

}  // namespace

PYBIND11_MODULE(sc_backend, m) {
    m.doc() = "Connolly-style SC backend prototype";

    m.def(
        "calculate_sc_backend",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> coords,
           py::array_t<double, py::array::c_style | py::array::forcecast> radii,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> molecule_ids,
           double probe_radius = 1.7,
           double dot_density = 15.0,
           double weight = 0.5,
           double trim = 1.5,
           double interface_distance = 8.0) {
            (void)interface_distance;
            auto atoms = load_atoms(coords, radii, molecule_ids);
            std::vector<Atom> atoms_m1;
            std::vector<Atom> atoms_m2;
            atoms_m1.reserve(atoms.size());
            atoms_m2.reserve(atoms.size());
            for (const auto& atom : atoms) {
                if (atom.molecule == 1) {
                    atoms_m1.push_back(atom);
                } else if (atom.molecule == 2) {
                    atoms_m2.push_back(atom);
                }
            }
            if (atoms_m1.empty() || atoms_m2.empty()) {
                throw std::runtime_error("Both molecule groups must be non-empty");
            }

            const auto grid1 = make_grid(atoms_m1, probe_radius, dot_density);
            const auto grid2 = make_grid(atoms_m2, probe_radius, dot_density);
            auto solid1 = build_vdw_volume(atoms_m1, grid1);
            auto solid2 = build_vdw_volume(atoms_m2, grid2);
            auto expanded1 = dilate_spherical(solid1, grid1, probe_radius);
            auto expanded2 = dilate_spherical(solid2, grid2, probe_radius);
            auto surface1 = extract_boundary_surface(expanded1, grid1, atoms_m1, probe_radius);
            auto surface2 = extract_boundary_surface(expanded2, grid2, atoms_m2, probe_radius);
            SurfaceStats stats;
            std::vector<SurfacePoint> combined_surface = surface1;
            combined_surface.insert(combined_surface.end(), surface2.begin(), surface2.end());
            const double sc = compute_sc_score(surface1, surface2, combined_surface, weight, trim, stats);

            py::dict out;
            out["n_surface_m1"] = py::int_(stats.n_surface_m1);
            out["n_surface_m2"] = py::int_(stats.n_surface_m2);
            out["n_buried_m1"] = py::int_(stats.n_buried_m1);
            out["n_buried_m2"] = py::int_(stats.n_buried_m2);
            out["n_trimmed_m1"] = py::int_(stats.n_trimmed_m1);
            out["n_trimmed_m2"] = py::int_(stats.n_trimmed_m2);
            out["grid_spacing_m1"] = py::float_(grid1.spacing);
            out["grid_spacing_m2"] = py::float_(grid2.spacing);
            out["n_total_surface"] = py::int_(surface1.size() + surface2.size());
            return py::make_tuple(sc, out);
        },
        py::arg("coords"),
        py::arg("radii"),
        py::arg("molecule_ids"),
        py::arg("probe_radius") = 1.7,
        py::arg("dot_density") = 15.0,
        py::arg("weight") = 0.5,
        py::arg("trim") = 1.5,
        py::arg("interface_distance") = 8.0);
}
