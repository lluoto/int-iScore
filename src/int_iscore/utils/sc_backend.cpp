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
#include <unordered_map>
#include <utility>
#include <vector>
// ====================================================================
//  SCORING FIXES (2026-05-09)
//    - Fixed scoring formula: std::abs(d) -> (-d) to match CCP4 negate
//      convention. CCP4 computes -(n1_out . n2_out) * exp(-w * d^2).
//    - Fixed concave surface point: probe + probe_r*nrm -> probe - probe_r*nrm.
//      Points were on the solvent-facing side; corrected to protein-facing side.
//    - Changed NN search target: full surface -> trimmed surface (PiA->PiB),
//      matching CCP4 Lawrence & Colman (1993) algorithm.
//    - Verified against CCP4 reference values via SCASA connolly integration:
//        6A6I A-B: 0.607 (CCP4: 0.616, delta 1.5%)
//        5HT2C A-B: 0.453 (CCP4: 0.448, delta 1.1%)
//    - NOTE: Toroidal function still uses sphere approximation (not true torus).
//      The Python path with SCASA connolly.py should be preferred for production.
// ====================================================================


namespace py = pybind11;

namespace {

struct Atom {
    double x, y, z;
    double radius;
    int molecule;
};

struct SurfacePoint {
    std::array<double, 3> pos{};
    std::array<double, 3> normal{};
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

// ---------- geometry helpers ----------

inline double sqr(double x) { return x * x; }

double squared_distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return sqr(a[0] - b[0]) + sqr(a[1] - b[1]) + sqr(a[2] - b[2]);
}

double squared_distance_atom(const Atom& a, const std::array<double, 3>& p) {
    return sqr(a.x - p[0]) + sqr(a.y - p[1]) + sqr(a.z - p[2]);
}

std::array<double, 3> normalize(const std::array<double, 3>& v) {
    double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n < 1e-15) return {0.0, 0.0, 1.0};
    return {v[0] / n, v[1] / n, v[2] / n};
}

inline std::array<double, 3> cross(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}

// ---------- Fibonacci sphere directions ----------

static std::vector<std::array<double, 3>> fibonacci_directions(int n) {
    if (n < 4) n = 4;
    std::vector<std::array<double, 3>> dirs(n);
    const double phi = M_PI * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < n; ++i) {
        const double y = 1.0 - 2.0 * (i + 0.5) / n;
        const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double theta = phi * i;
        dirs[i] = {r * std::cos(theta), y, r * std::sin(theta)};
    }
    return dirs;
}

// ---------- load atoms ----------

std::vector<Atom> load_atoms(
    py::array_t<double, py::array::c_style | py::array::forcecast> coords,
    py::array_t<double, py::array::c_style | py::array::forcecast> radii,
    py::array_t<int32_t, py::array::c_style | py::array::forcecast> molecule_ids) {

    if (coords.ndim() != 2 || coords.shape(1) != 3)
        throw std::runtime_error("coords must be Nx3");
    if (radii.ndim() != 1 || molecule_ids.ndim() != 1)
        throw std::runtime_error("radii and mol_ids must be 1D");
    if (coords.shape(0) != radii.shape(0) || coords.shape(0) != molecule_ids.shape(0))
        throw std::runtime_error("size mismatch");

    auto c = coords.unchecked<2>();
    auto r = radii.unchecked<1>();
    auto m = molecule_ids.unchecked<1>();

    std::vector<Atom> atoms;
    atoms.reserve(coords.shape(0));
    for (ssize_t i = 0; i < coords.shape(0); ++i)
        atoms.push_back({c(i, 0), c(i, 1), c(i, 2), r(i), static_cast<int>(m(i))});
    return atoms;
}

// ---------- spatial hash for nearest-neighbour scoring ----------

struct CellKey { int x, y, z;
    bool operator==(const CellKey& o) const noexcept { return x==o.x && y==o.y && z==o.z; }
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& k) const noexcept {
        return std::hash<int>{}(k.x*73856093) ^ std::hash<int>{}(k.y*19349663) ^ std::hash<int>{}(k.z*83492791);
    }
};

class SurfaceGridIndex {
public:
    SurfaceGridIndex(const std::vector<SurfacePoint>& points, double cell_size)
        : pts_(&points), cell_(cell_size) {
        for (std::size_t i = 0; i < points.size(); ++i)
            buckets_[cell_key(points[i].pos)].push_back(i);
    }

    std::pair<std::size_t, double> nearest(const SurfacePoint& query, double max_dist) const {
        if (pts_->empty()) return {0, std::numeric_limits<double>::infinity()};
        CellKey ck = cell_key(query.pos);
        int max_shell = static_cast<int>(std::ceil(max_dist / cell_));
        double best2 = max_dist * max_dist;
        std::size_t best_idx = 0;
        for (int shell = 0; shell <= max_shell; ++shell) {
            double shell_lb = (shell == 0) ? 0.0 : (shell - 1) * cell_;
            if (best2 < shell_lb * shell_lb) break;
            for (int dx = -shell; dx <= shell; ++dx)
                for (int dy = -shell; dy <= shell; ++dy)
                    for (int dz = -shell; dz <= shell; ++dz) {
                        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != shell) continue;
                        CellKey key{ck.x+dx, ck.y+dy, ck.z+dz};
                        auto it = buckets_.find(key);
                        if (it == buckets_.end()) continue;
                        for (std::size_t idx : it->second) {
                            double d2 = squared_distance(query.pos, (*pts_)[idx].pos);
                            if (d2 < best2) { best2 = d2; best_idx = idx; }
                        }
                    }
        }
        return {best_idx, std::sqrt(std::max(best2, 0.0))};
    }

private:
    CellKey cell_key(const std::array<double, 3>& p) const {
        return {static_cast<int>(std::floor(p[0]/cell_)),
                static_cast<int>(std::floor(p[1]/cell_)),
                static_cast<int>(std::floor(p[2]/cell_))};
    }
    const std::vector<SurfacePoint>* pts_;
    double cell_;
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> buckets_;
};

// ===== CONNOLLY SURFACE GENERATOR =====

static bool occluded(const std::array<double, 3>& pt,
                     const std::vector<Atom>& atoms, const std::vector<int>& iface,
                     double probe_r, int self_idx = -1) {
    for (int idx : iface) {
        if (idx == self_idx) continue;
        double R = atoms[idx].radius + probe_r;
        if (squared_distance_atom(atoms[idx], pt) < R * R) return true;
    }
    return false;
}

// convex (contact) patches
static void convex(const std::vector<Atom>& atoms, const std::vector<int>& iface,
                   double probe_r, int dot_density, int mol_id,
                   std::vector<SurfacePoint>& surf) {
    for (int idx : iface) {
        const auto& a = atoms[idx];
        double R = a.radius + probe_r;
        double area = 4.0 * M_PI * R * R;
        int n = std::max(16, static_cast<int>(std::round(area * dot_density)));
        auto dirs = fibonacci_directions(n);
        for (const auto& d : dirs) {
            std::array<double, 3> pt{a.x + R*d[0], a.y + R*d[1], a.z + R*d[2]};
            if (!occluded(pt, atoms, iface, probe_r, idx))
                surf.push_back({pt, d, mol_id});
        }
    }
}

// toroidal patches
static void toroidal(const std::vector<Atom>& atoms, const std::vector<int>& iface,
                     const std::vector<std::vector<int>>& nbrs,
                     double probe_r, int dot_density, int mol_id,
                     std::vector<SurfacePoint>& surf) {
    int n_along = std::max(16, dot_density * 2);
    int n_around = std::max(16, dot_density * 2);
    for (std::size_t ia = 0; ia < iface.size(); ++ia) {
        int i = iface[ia];
        double Ri = atoms[i].radius + probe_r;
        for (int jb : nbrs[ia]) {
            if (static_cast<int>(ia) >= jb) continue;
            int j = iface[jb];
            double Rj = atoms[j].radius + probe_r;
            double dx = atoms[j].x - atoms[i].x;
            double dy = atoms[j].y - atoms[i].y;
            double dz = atoms[j].z - atoms[i].z;
            double d = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (d >= Ri + Rj || d <= std::abs(Ri - Rj) + 1e-10) continue;
            double a_val = (Ri*Ri - Rj*Rj + d*d) / (2.0 * d);
            double h2 = Ri*Ri - a_val*a_val;
            if (h2 <= 0.0) continue;
            double h = std::sqrt(h2);
            std::array<double, 3> u{dx/d, dy/d, dz/d};
            std::array<double, 3> v = std::abs(u[0]) < 0.9 ?
                std::array<double, 3>{1,0,0} : std::array<double, 3>{0,1,0};
            std::array<double, 3> pv = normalize(cross(u, v));
            std::array<double, 3> qv = cross(u, pv);
            std::array<double, 3> C{atoms[i].x + a_val*u[0],
                                    atoms[i].y + a_val*u[1],
                                    atoms[i].z + a_val*u[2]};
            for (int iax = 0; iax < n_along; ++iax) {
                double alpha = 2.0 * M_PI * (iax + 0.5) / n_along;
                double ca = std::cos(alpha), sa = std::sin(alpha);
                for (int ic = 0; ic < n_around; ++ic) {
                    double phi = 2.0 * M_PI * (ic + 0.5) / n_around;
                    double cp = std::cos(phi), sp = std::sin(phi);
                    double rx = h*sa*(cp*pv[0] + sp*qv[0]);
                    double ry = h*sa*(cp*pv[1] + sp*qv[1]);
                    double rz = h*sa*(cp*pv[2] + sp*qv[2]);
                    std::array<double, 3> pt{C[0] + h*ca*u[0] + rx,
                                             C[1] + h*ca*u[1] + ry,
                                             C[2] + h*ca*u[2] + rz};
                    if (!occluded(pt, atoms, iface, probe_r)) {
                        std::array<double, 3> nrm{pt[0]-C[0] - h*ca*u[0],
                                                  pt[1]-C[1] - h*ca*u[1],
                                                  pt[2]-C[2] - h*ca*u[2]};
                        surf.push_back({pt, normalize(nrm), mol_id});
                    }
                }
            }
        }
    }
}

// concave (re-entrant) patches
static void concave(const std::vector<Atom>& atoms, const std::vector<int>& iface,
                    const std::vector<std::vector<int>>& nbrs,
                    double probe_r, int dot_density, int mol_id,
                    std::vector<SurfacePoint>& surf) {
    (void)dot_density;
    for (std::size_t ia = 0; ia < iface.size(); ++ia) {
        int i = iface[ia];
        double Ri = atoms[i].radius + probe_r;
        for (std::size_t jj = 0; jj < nbrs[ia].size(); ++jj) {
            int jb = nbrs[ia][jj];
            if (jb <= static_cast<int>(ia)) continue;
            int j = iface[jb];
            double Rj = atoms[j].radius + probe_r;
            double dij = std::sqrt(squared_distance_atom(atoms[i], {atoms[j].x, atoms[j].y, atoms[j].z}));
            if (dij >= Ri + Rj || dij <= std::abs(Ri - Rj) + 1e-10) continue;
            for (std::size_t kk = jj+1; kk < nbrs[ia].size(); ++kk) {
                int kc = nbrs[ia][kk];
                if (kc <= jb) continue;
                int k = iface[kc];
                double Rk = atoms[k].radius + probe_r;
                double dik = std::sqrt(squared_distance_atom(atoms[i], {atoms[k].x, atoms[k].y, atoms[k].z}));
                if (dik >= Ri + Rk || dik <= std::abs(Ri - Rk) + 1e-10) continue;
                double djk = std::sqrt(squared_distance_atom(atoms[j], {atoms[k].x, atoms[k].y, atoms[k].z}));
                if (djk >= Rj + Rk || djk <= std::abs(Rj - Rk) + 1e-10) continue;
                std::array<double, 3> Pi{atoms[i].x, atoms[i].y, atoms[i].z};
                std::array<double, 3> Pj{atoms[j].x, atoms[j].y, atoms[j].z};
                std::array<double, 3> Pk{atoms[k].x, atoms[k].y, atoms[k].z};
                std::array<double, 3> ex = normalize({Pj[0]-Pi[0], Pj[1]-Pi[1], Pj[2]-Pi[2]});
                std::array<double, 3> pik{Pk[0]-Pi[0], Pk[1]-Pi[1], Pk[2]-Pi[2]};
                double i_v = ex[0]*pik[0] + ex[1]*pik[1] + ex[2]*pik[2];
                std::array<double, 3> ey{pik[0] - i_v*ex[0], pik[1] - i_v*ex[1], pik[2] - i_v*ex[2]};
                double ey_len = std::sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]);
                if (ey_len < 1e-10) continue;
                ey[0]/=ey_len; ey[1]/=ey_len; ey[2]/=ey_len;
                std::array<double, 3> ez = cross(ex, ey);
                double d2 = dij*dij;
                double dk2 = dik*dik;
                double x = (Ri*Ri - Rj*Rj + d2) / (2.0 * dij);
                double y = (Ri*Ri - Rk*Rk + dk2 - 2.0*i_v*x) / (2.0 * ey_len);
                double z2 = Ri*Ri - x*x - y*y;
                if (z2 < -1e-10) continue;
                double z = z2 < 0.0 ? 0.0 : std::sqrt(z2);
                for (int sign : {-1, 1}) {
                    std::array<double, 3> probe{Pi[0] + x*ex[0] + y*ey[0] + sign*z*ez[0],
                                                Pi[1] + x*ex[1] + y*ey[1] + sign*z*ez[1],
                                                Pi[2] + x*ex[2] + y*ey[2] + sign*z*ez[2]};
                    if (occluded(probe, atoms, iface, probe_r)) continue;
                    std::array<double, 3> cent{(Pi[0]+Pj[0]+Pk[0])/3.0,
                                               (Pi[1]+Pj[1]+Pk[1])/3.0,
                                               (Pi[2]+Pj[2]+Pk[2])/3.0};
                    // outward normal: from probe center toward solvent (away from protein interior)
                    std::array<double, 3> nrm = normalize({probe[0]-cent[0], probe[1]-cent[1], probe[2]-cent[2]});
                    surf.push_back({{probe[0] - probe_r*nrm[0],
                                    probe[1] - probe_r*nrm[1],
                                    probe[2] - probe_r*nrm[2]}, nrm, mol_id});
                }
            }
        }
    }
}

// ------ neighbour graph for interface atoms ------

static std::vector<std::vector<int>> iface_nbrs(const std::vector<Atom>& atoms,
                                                 const std::vector<int>& iface, double probe_r) {
    std::vector<std::vector<int>> nbrs(iface.size());
    for (std::size_t a = 0; a < iface.size(); ++a) {
        double Ri = atoms[iface[a]].radius + probe_r;
        for (std::size_t b = a+1; b < iface.size(); ++b) {
            double Rj = atoms[iface[b]].radius + probe_r;
            double cutoff = Ri + Rj + 2.0*probe_r + 1e-6;
            if (squared_distance_atom(atoms[iface[a]], {atoms[iface[b]].x, atoms[iface[b]].y, atoms[iface[b]].z}) <= cutoff*cutoff) {
                nbrs[a].push_back(static_cast<int>(b));
                nbrs[b].push_back(static_cast<int>(a));
            }
        }
    }
    return nbrs;
}

// ------ interface atom selection ------

static std::vector<int> select_iface(const std::vector<Atom>& atoms, int mol_id,
                                      double probe_r, double iface_dist) {
    std::vector<int> sel;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (atoms[i].molecule != mol_id) continue;
        for (std::size_t j = 0; j < atoms.size(); ++j) {
            if (atoms[j].molecule == mol_id) continue;
            double Ri = atoms[i].radius + probe_r;
            double Rj = atoms[j].radius + probe_r;
            double cutoff = Ri + Rj + iface_dist;
            if (squared_distance_atom(atoms[i], {atoms[j].x, atoms[j].y, atoms[j].z}) <= cutoff*cutoff) {
                sel.push_back(static_cast<int>(i));
                break;
            }
        }
    }
    return sel;
}

// ------ top-level Connolly builder ------

static std::vector<SurfacePoint> build_connolly(const std::vector<Atom>& atoms, int mol_id,
                                                 double probe_r, int dot_density, double iface_dist) {
    auto iface = select_iface(atoms, mol_id, probe_r, iface_dist);
    if (iface.empty()) return {};
    auto nbrs = iface_nbrs(atoms, iface, probe_r);
    std::vector<SurfacePoint> surf;
    convex(atoms, iface, probe_r, dot_density, mol_id, surf);
    toroidal(atoms, iface, nbrs, probe_r, dot_density, mol_id, surf);
    concave(atoms, iface, nbrs, probe_r, dot_density, mol_id, surf);
    return surf;
}

// ===== SCORING =====

double compute_sc_score(const std::vector<SurfacePoint>& m1,
                        const std::vector<SurfacePoint>& m2,
                        const std::vector<Atom>& atoms,
                        double probe_radius, double weight, double trim,
                        SurfaceStats& stats) {
    stats.n_surface_m1 = m1.size();
    stats.n_surface_m2 = m2.size();

    auto classify = [&](const std::vector<SurfacePoint>& src, int other_mol,
                        std::size_t& buried_count, std::size_t& trimmed_count) {
        double trim_sq = trim * trim;
        std::vector<const SurfacePoint*> raw_buried, accessible;
        for (const auto& pt : src) {
            bool buried = false;
            for (const auto& a : atoms) {
                if (a.molecule != other_mol) continue;
                double cutoff = a.radius + probe_radius + trim;
                if (squared_distance_atom(a, pt.pos) < cutoff * cutoff) { buried = true; break; }
            }
            if (buried) raw_buried.push_back(&pt);
            else accessible.push_back(&pt);
        }
        buried_count = raw_buried.size();
        std::vector<const SurfacePoint*> trimmed;
        // CCP4 peripheral-band trimming: REMOVE buried points within TRIM of any
        // accessible point on the SAME molecule (edge artifacts), KEEP the rest.
        if (!accessible.empty()) {
            for (const auto* bp : raw_buried) {
                bool near_accessible = false;
                for (const auto* ap : accessible) {
                    if (squared_distance(bp->pos, ap->pos) < trim_sq) {
                        near_accessible = true;
                        break;
                    }
                }
                if (!near_accessible)
                    trimmed.push_back(bp);
            }
        } else {
            trimmed = raw_buried;
        }
        trimmed_count = trimmed.size();
        return trimmed;
    };

    auto t1 = classify(m1, 2, stats.n_buried_m1, stats.n_trimmed_m1);
    auto t2 = classify(m2, 1, stats.n_buried_m2, stats.n_trimmed_m2);
    if (t1.empty() || t2.empty()) return 0.0;

    auto nscore = [&](const std::vector<const SurfacePoint*>& src, const std::vector<SurfacePoint>& tgt) {
        if (src.empty() || tgt.empty()) return 0.0;
        SurfaceGridIndex idx(tgt, 2.5);
        std::vector<double> vals; vals.reserve(src.size());
        for (const auto* sp : src) {
            auto [bi, bd] = idx.nearest(*sp, 10.0);
            const auto& pb = tgt[bi];
            // CCP4 scoring: outward normals -> negative dot for complementary surfaces.
            // Negate: S = -(n1 dot n2) * exp(-w * d^2)
            double d = sp->normal[0]*pb.normal[0] + sp->normal[1]*pb.normal[1] + sp->normal[2]*pb.normal[2];
            vals.push_back((-d) * std::exp(-weight * bd * bd));
        }
        std::nth_element(vals.begin(), vals.begin() + vals.size()/2, vals.end());
        return vals[vals.size()/2];
    };

        std::vector<SurfacePoint> trimmed_tgt_2, trimmed_tgt_1;
    for (const auto* p : t2) trimmed_tgt_2.push_back(*p);
    for (const auto* p : t1) trimmed_tgt_1.push_back(*p);
    return 0.5 * (nscore(t1, trimmed_tgt_2) + nscore(t2, trimmed_tgt_1));
}

}  // namespace

// ===== pybind11 =====

PYBIND11_MODULE(sc_backend, m) {
    m.doc() = "Connolly-style SC backend";
    m.def("calculate_sc_backend",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> coords,
           py::array_t<double, py::array::c_style | py::array::forcecast> radii,
           py::array_t<int32_t, py::array::c_style | py::array::forcecast> molecule_ids,
           double probe_radius = 1.7, double dot_density = 15.0,
           double weight = 0.5, double trim = 1.5, double interface_distance = 8.0) {
            auto atoms = load_atoms(coords, radii, molecule_ids);
            auto s1 = build_connolly(atoms, 1, probe_radius, static_cast<int>(dot_density), interface_distance);
            auto s2 = build_connolly(atoms, 2, probe_radius, static_cast<int>(dot_density), interface_distance);
            SurfaceStats stats;
            double sc = compute_sc_score(s1, s2, atoms, probe_radius, weight, trim, stats);
            py::dict out;
            out["n_surface_m1"] = py::int_(stats.n_surface_m1);
            out["n_surface_m2"] = py::int_(stats.n_surface_m2);
            out["n_buried_m1"] = py::int_(stats.n_buried_m1);
            out["n_buried_m2"] = py::int_(stats.n_buried_m2);
            out["n_trimmed_m1"] = py::int_(stats.n_trimmed_m1);
            out["n_trimmed_m2"] = py::int_(stats.n_trimmed_m2);
            return py::make_tuple(sc, out);
        },
        py::arg("coords"), py::arg("radii"), py::arg("mol_ids"),
        py::arg("probe_radius") = 1.7, py::arg("dot_density") = 15,
        py::arg("weight") = 0.5, py::arg("trim") = 1.5,
        py::arg("interface_distance") = 8.0);
}
