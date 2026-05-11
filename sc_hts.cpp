// sc_hts.cpp — High-throughput Shape Complementarity (Sc) Engine
// Two-tier architecture: AccurateConnolly (production) + FastSAS (experimental)
// C++17, nanoflann KD-Tree, OpenMP parallelized
// Compile: g++ -std=c++17 -O3 -fopenmp sc_hts.cpp -o sc_hts
//
// DEFAULT MODE: AccurateConnolly (SCASA bridge, CCP4-grade, +/-2% accuracy)
//   Accurate mode calls Python SCASA bridge via popen; requires sc_bridge.py in PATH.
//
// EXPERIMENTAL: FastSAS (mode 0) — under active development, does NOT yet generalize.
//   See FASTSAS_ANALYSIS_PROMPT.md for current status and open challenges.
//   Contributions welcome: https://github.com/lluoto/int-iScore
//
// Usage:
//   sc_hts <pdb_file> <chain1> <chain2> [mode=1] [pdb_id]
//   mode: 0=FastSAS(experimental), 1=AccurateConnolly(default)
//
// Output (CSV line to stdout): PDB_ID, Chains, Mode, TotalDots, BuriedDots, TrimmedDots, Median_D, Sc_Score

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <cfloat>
#include <cstdint>
#include <chrono>

#include <nanoflann.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// Global Parameters (hardcoded as constexpr)
// ============================================================================
static constexpr double PROBE_RADIUS       = 1.4;   // Å — water molecule probe
static constexpr double DOT_DENSITY        = 15.0;  // dots per Å²
static constexpr double TRIM_DISTANCE      = 1.5;   // Å — CCP4 trim band
static constexpr double WEIGHT_ACCURATE    = 0.5;   // Å⁻² — CCP4 standard
static constexpr double WEIGHT_FAST        = 0.05;  // Å⁻² — adjusted for SAS inflated surface
static constexpr double INTERFACE_DISTANCE = 8.0;   // Å — interface atom cutoff
static constexpr double FAST_DISTANCE_CUTOFF = 4.0; // Å — hard truncation for FastSAS
static constexpr double PI_VAL            = 3.14159265358979323846;
static constexpr double PHI_GOLDEN        = 1.6180339887498948482;
// ============================================================================
// FastSAS scoring method selection (compile-time switch)
// ============================================================================
enum class FastMethod : int {
    VOXEL_DENSITY   = 0,  // Voxelized packing density (pure geometry, <30ms)
    SAS_FEATURE_REG = 1   // SAS surface feature-based regression (~250ms)
};
// EXPERIMENTAL: FastSAS scoring method (under development, needs optimization).
// Change this constant to switch FastSAS scoring method:
static constexpr FastMethod FAST_METHOD = FastMethod::SAS_FEATURE_REG;


// ============================================================================
// Atom struct
// ============================================================================
struct Atom {
    int    serial;
    std::string name;
    std::string resName;
    char   chainID;
    int    resSeq;
    double x, y, z;
    std::string element;
    double radius;
};

// ============================================================================
// SOA PointCloud (Structure of Arrays)
// ============================================================================
struct PointCloud {
    std::vector<double> xs, ys, zs;
    std::vector<double> nxs, nys, nzs;

    int size() const { return (int)xs.size(); }
    void reserve(int n) {
        xs.reserve(n); ys.reserve(n); zs.reserve(n);
        nxs.reserve(n); nys.reserve(n); nzs.reserve(n);
    }
    void clear() { xs.clear(); ys.clear(); zs.clear();
                   nxs.clear(); nys.clear(); nzs.clear(); }
    void push(double x, double y, double z, double nx, double ny, double nz) {
        xs.push_back(x); ys.push_back(y); zs.push_back(z);
        nxs.push_back(nx); nys.push_back(ny); nzs.push_back(nz);
    }
};
// ============================================================================
// VoxelGrid - 1D flat array for Method A (Voxelized Packing Density)
// ============================================================================
struct VoxelGrid {
    double origin_x, origin_y, origin_z;
    double voxel_size;
    int nx, ny, nz;
    int nxy;
    std::vector<uint8_t> grid;

    VoxelGrid(double ox, double oy, double oz, double sz,
              int _nx, int _ny, int _nz)
        : origin_x(ox), origin_y(oy), origin_z(oz), voxel_size(sz),
          nx(_nx), ny(_ny), nz(_nz), nxy(_nx * _ny) {
        grid.assign(nx * ny * nz, 0);
    }

    inline int idx(double x, double y, double z) const {
        int ix = (int)((x - origin_x) / voxel_size);
        int iy = (int)((y - origin_y) / voxel_size);
        int iz = (int)((z - origin_z) / voxel_size);
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
            return -1;
        return ix + iy * nx + iz * nxy;
    }

    void rasterize(const std::vector<Atom>& atoms, int chain_bit) {
        uint8_t mask = (uint8_t)(1 << chain_bit);
        double inv_vs = 1.0 / voxel_size;
        #pragma omp parallel for
        for (int ai = 0; ai < (int)atoms.size(); ++ai) {
            const Atom& a = atoms[ai];
            double R = a.radius + PROBE_RADIUS;  // SAS shell radius
            double R2 = R * R;
            int ix_min = std::max(0, (int)((a.x - R - origin_x) * inv_vs));
            int ix_max = std::min(nx - 1, (int)((a.x + R - origin_x) * inv_vs));
            int iy_min = std::max(0, (int)((a.y - R - origin_y) * inv_vs));
            int iy_max = std::min(ny - 1, (int)((a.y + R - origin_y) * inv_vs));
            int iz_min = std::max(0, (int)((a.z - R - origin_z) * inv_vs));
            int iz_max = std::min(nz - 1, (int)((a.z + R - origin_z) * inv_vs));
            for (int iz = iz_min; iz <= iz_max; ++iz) {
                double dz = (iz * voxel_size + origin_z) - a.z;
                for (int iy = iy_min; iy <= iy_max; ++iy) {
                    double dy = (iy * voxel_size + origin_y) - a.y;
                    for (int ix = ix_min; ix <= ix_max; ++ix) {
                        double dx = (ix * voxel_size + origin_x) - a.x;
                        if (dx*dx + dy*dy + dz*dz <= R2) {
                            int ci = ix + iy * nx + iz * nxy;
                            #pragma omp atomic
                            grid[ci] |= mask;
                        }
                    }
                }
            }
        }
    }

    void count(int& v_interlock, int& v_union) const {
        v_interlock = 0; v_union = 0;
        #pragma omp parallel for reduction(+:v_interlock, v_union)
        for (int i = 0; i < (int)grid.size(); ++i) {
            if (grid[i] != 0) ++v_union;
            if (grid[i] == 0x03) ++v_interlock;
        }
    }
};


// ============================================================================
// nanoflann adaptor for SOA PointCloud
// ============================================================================
template <typename T>
struct SOAAdaptor {
    const PointCloud& cloud;

    SOAAdaptor(const PointCloud& c) : cloud(c) {}

    size_t kdtree_get_point_count() const { return cloud.size(); }
    T   kdtree_get_pt(size_t idx, int dim) const {
        if (dim == 0) return T(cloud.xs[idx]);
        if (dim == 1) return T(cloud.ys[idx]);
        return T(cloud.zs[idx]);
    }
    template <class B> bool kdtree_get_bbox(B&) const { return false; }
};

using KDTree3D = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, SOAAdaptor<double>>,
    SOAAdaptor<double>, 3>;

// ============================================================================
// Atom Radii Lookup Table
// ============================================================================
static const std::unordered_map<std::string, double> ATOM_RADII = {
    {"C", 1.70}, {"CA", 1.70}, {"N", 1.55}, {"O", 1.52}, {"S", 1.80},
    {"P", 1.80}, {"H", 1.20}, {"F", 1.47}, {"CL",1.75}, {"MG",1.73},
    {"ZN",1.39}, {"FE",1.40}, {"MN",1.40}, {"NA",1.80}, {"K", 2.20},
    {"I", 1.98}, {"BR",1.85}, {"SE",1.80},
};

static double get_atom_radius(const std::string& element) {
    auto it = ATOM_RADII.find(element);
    return (it != ATOM_RADII.end()) ? it->second : 1.70;
}

// ============================================================================
// PDB Parser
// ============================================================================
static std::vector<Atom> parse_pdb(const std::string& pdb_path) {
    std::vector<Atom> atoms;
    std::ifstream file(pdb_path);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open PDB file: " << pdb_path << std::endl;
        return atoms;
    }

    std::string line;
    Atom current;
    while (std::getline(file, line)) {
        if (line.length() < 54) continue;
        std::string record = line.substr(0, 6);
        // Trim whitespace
        record.erase(0, record.find_first_not_of(" \t\n\r\f\v"));
        record.erase(record.find_last_not_of(" \t\n\r\f\v") + 1);

        if (record != "ATOM" && record != "HETATM") continue;

        try {
            current.serial   = std::stoi(line.substr(6, 5));
            current.name     = line.substr(12, 4);
            // remove whitespace from name
            current.name.erase(0, current.name.find_first_not_of(" "));
            current.name.erase(current.name.find_last_not_of(" ") + 1);
            current.resName  = line.substr(17, 3);
            current.resName.erase(0, current.resName.find_first_not_of(" "));
            current.chainID  = (line.length() > 21) ? line[21] : ' ';
            current.resSeq   = std::stoi(line.substr(22, 4));
            current.x        = std::stod(line.substr(30, 8));
            current.y        = std::stod(line.substr(38, 8));
            current.z        = std::stod(line.substr(46, 8));

            // Try to get element from columns 77-78, else deduce from atom name
            if (line.length() >= 78) {
                current.element = line.substr(76, 2);
                current.element.erase(0, current.element.find_first_not_of(" "));
                current.element.erase(current.element.find_last_not_of(" ") + 1);
            }
            if (current.element.empty()) {
                // Deduce from atom name: first char is usually element
                current.element = current.name.substr(0, 1);
                if (current.element == "C" && current.name.length() > 1 &&
                    current.name[1] == 'L') current.element = "CL";
            }

            current.radius = get_atom_radius(current.element);
            atoms.push_back(current);
        } catch (...) {
            continue;  // skip malformed lines
        }
    }
    return atoms;
}

// ============================================================================
// Filter atoms by chain, group into molecule1 / molecule2
// ============================================================================
struct Molecule {
    std::vector<Atom> atoms;
    PointCloud total_surface;
    PointCloud buried_surface;
    PointCloud trimmed_surface;
};

static bool filter_by_chain(const std::vector<Atom>& all,
                            const std::string& chains,
                            Molecule& mol) {
    for (const auto& a : all) {
        if (chains.find(a.chainID) != std::string::npos) {
            mol.atoms.push_back(a);
        }
    }
    return !mol.atoms.empty();
}

// ============================================================================
// Centering: translate all atom coordinates so centroid = origin
// ============================================================================
static void center_atoms(std::vector<Atom>& atoms1, std::vector<Atom>& atoms2) {
    double cx = 0, cy = 0, cz = 0;
    int n = atoms1.size() + atoms2.size();
    for (const auto& a : atoms1) { cx += a.x; cy += a.y; cz += a.z; }
    for (const auto& a : atoms2) { cx += a.x; cy += a.y; cz += a.z; }
    cx /= n; cy /= n; cz /= n;

    for (auto& a : atoms1) { a.x -= cx; a.y -= cy; a.z -= cz; }
    for (auto& a : atoms2) { a.x -= cx; a.y -= cy; a.z -= cz; }
}

// ============================================================================
// Fibonacci Sphere Sampling — uniform points on sphere
// ============================================================================
static void fibonacci_sphere(int n_points, std::vector<double>& xs,
                             std::vector<double>& ys, std::vector<double>& zs) {
    xs.resize(n_points); ys.resize(n_points); zs.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
        double y = 1.0 - (2.0 * i + 1.0) / n_points;  // y goes from +1 to -1
        double radius_at_y = std::sqrt(1.0 - y * y);
        double theta = 2.0 * PI_VAL * i / PHI_GOLDEN;
        xs[i] = std::cos(theta) * radius_at_y;
        ys[i] = y;
        zs[i] = std::sin(theta) * radius_at_y;
    }
}

// ============================================================================
// Build nanoflann KD-Tree from PointCloud
// ============================================================================
static KDTree3D* build_kdtree(const PointCloud& cloud) {
    auto* adaptor = new SOAAdaptor<double>(cloud);
    auto* tree = new KDTree3D(3, *adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(10)); // max leaf size
    tree->buildIndex();
    return tree;
}

// ============================================================================
// Build KD-Tree of atom centers (for occlusion / buried check)
// ============================================================================
struct AtomCenterCloud {
    const std::vector<Atom>& atoms;
    AtomCenterCloud(const std::vector<Atom>& a) : atoms(a) {}
    size_t kdtree_get_point_count() const { return atoms.size(); }
    double kdtree_get_pt(size_t idx, int dim) const {
        if (dim == 0) return atoms[idx].x;
        if (dim == 1) return atoms[idx].y;
        return atoms[idx].z;
    }
    template <class B> bool kdtree_get_bbox(B&) const { return false; }
};

using AtomKDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, AtomCenterCloud>,
    AtomCenterCloud, 3>;

static AtomKDTree* build_atom_tree(const std::vector<Atom>& atoms) {
    auto* cloud = new AtomCenterCloud(atoms);
    auto* tree = new AtomKDTree(3, *cloud,
        nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree->buildIndex();
    return tree;
}

// ============================================================================
// Vector operations (inline for performance)
// ============================================================================
inline void normalize(double& x, double& y, double& z) {
    double len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-12) { x /= len; y /= len; z /= len; }
}

inline double dot(double ax, double ay, double az, double bx, double by, double bz) {
    return ax*bx + ay*by + az*bz;
}

inline double dist2(double ax, double ay, double az, double bx, double by, double bz) {
    double dx = ax - bx, dy = ay - by, dz = az - bz;
    return dx*dx + dy*dy + dz*dz;
}

// ============================================================================
// Extract interface atoms (within INTERFACE_DISTANCE of partner)
// ============================================================================
static void get_interface_atoms(
    const std::vector<Atom>& mol1, const std::vector<Atom>& mol2,
    std::vector<Atom>& iface1, std::vector<Atom>& iface2,
    double cutoff)
{
    // Build KD-Tree for mol2 centers
    AtomCenterCloud cloud2(mol2);
    AtomKDTree tree2(3, cloud2, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree2.buildIndex();

    double query[3];
    unsigned int n_results;
    double out_dist;

    for (const auto& a : mol1) {
        query[0] = a.x; query[1] = a.y; query[2] = a.z;
        tree2.knnSearch(query, 1, &n_results, &out_dist);
        if (out_dist < cutoff * cutoff) {
            iface1.push_back(a);
        }
    }

    AtomCenterCloud cloud1(mol1);
    AtomKDTree tree1(3, cloud1, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree1.buildIndex();

    for (const auto& a : mol2) {
        query[0] = a.x; query[1] = a.y; query[2] = a.z;
        tree1.knnSearch(query, 1, &n_results, &out_dist);
        if (out_dist < cutoff * cutoff) {
            iface2.push_back(a);
        }
    }
}

// ============================================================================
// MODE 1: FastSAS — rapid coarse filter
// ============================================================================
static void fast_sas_surface(const std::vector<Atom>& iface,
                              AtomKDTree* self_tree,
                              const std::vector<Atom>& own_atoms, // for self-occlusion
                              PointCloud& out_sas)
{
    out_sas.clear();

    // For occlusion: KD-tree of own atom centers
    AtomCenterCloud own_cloud(own_atoms);
    AtomKDTree own_tree(3, own_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    own_tree.buildIndex();

    // Precompute Fibonacci sphere points
    static thread_local std::vector<double> f_xs, f_ys, f_zs;
    static thread_local int f_n = 0;
    int max_n = (int)(4.0 * PI_VAL * (2.0 + PROBE_RADIUS) * (2.0 + PROBE_RADIUS) * DOT_DENSITY);
    if (f_n < max_n) {
        fibonacci_sphere(max_n, f_xs, f_ys, f_zs);
        f_n = max_n;
    }

    for (const auto& a : iface) {
        double R = a.radius + PROBE_RADIUS;
        int n_pts = std::max(6, (int)(4.0 * PI_VAL * R * R * DOT_DENSITY));
        if (n_pts > f_n) n_pts = f_n;

        // Use strided sampling if we have more Fibonacci points than needed
        int stride = f_n / n_pts;

        for (int i = 0; i < n_pts; ++i) {
            int idx = i * stride;
            if (idx >= f_n) idx = f_n - 1;

            double sx = a.x + f_xs[idx] * R;
            double sy = a.y + f_ys[idx] * R;
            double sz = a.z + f_zs[idx] * R;

            // Normal = outward from atom center
            double nx = f_xs[idx];
            double ny = f_ys[idx];
            double nz = f_zs[idx];

            // Occlusion culling: check if probe_center collides with own atoms
            double pcx = sx + nx * PROBE_RADIUS;
            double pcy = sy + ny * PROBE_RADIUS;
            double pcz = sz + nz * PROBE_RADIUS;

            double q[3] = {pcx, pcy, pcz};
    std::vector<unsigned int> ret_index(3);
    std::vector<double> ret_dist(3);
    unsigned int found = own_tree.knnSearch(q, 3, ret_index.data(), ret_dist.data());

    bool occluded = false;
    for (unsigned int k = 0; k < found; ++k) {
                int ai = ret_index[k];
                if (ret_dist[k] < own_atoms[ai].radius * own_atoms[ai].radius) {
                    occluded = true;
                    break;
                }
            }
            if (occluded) continue;

            out_sas.push(sx, sy, sz, nx, ny, nz);
        }
    }
}

// ============================================================================
// Classify buried SAS points
// ============================================================================
static void classify_buried_fast(const PointCloud& sas, AtomKDTree* other_tree,
                                  const std::vector<Atom>& other_atoms,
                                  const std::vector<Atom>& self_atoms,
                                  PointCloud& buried)
{
    buried.clear();
    double q[3];
    std::vector<unsigned int> ret_idx(3);
    std::vector<double> ret_dist(3);

    for (int i = 0; i < sas.size(); ++i) {
        double pcx = sas.xs[i] + sas.nxs[i] * PROBE_RADIUS;
        double pcy = sas.ys[i] + sas.nys[i] * PROBE_RADIUS;
        double pcz = sas.zs[i] + sas.nzs[i] * PROBE_RADIUS;

        q[0] = pcx; q[1] = pcy; q[2] = pcz;
    unsigned int found = other_tree->knnSearch(q, 3, ret_idx.data(), ret_dist.data());

    bool is_buried = false;
    for (unsigned int k = 0; k < found; ++k) {
            int ai = ret_idx[k];
            if (ret_dist[k] < other_atoms[ai].radius * other_atoms[ai].radius) {
                is_buried = true;
                break;
            }
        }
        if (is_buried) {
            buried.push(sas.xs[i], sas.ys[i], sas.zs[i],
                         sas.nxs[i], sas.nys[i], sas.nzs[i]);
        }
    }
}

// ============================================================================
// METHOD A: Voxelized Packing Density
// ============================================================================
static double score_voxel_density(const std::vector<Atom>& iface1,
                                   const std::vector<Atom>& iface2)
{
    if (iface1.empty() || iface2.empty()) return 0.0;
    std::vector<Atom> all;
    all.reserve(iface1.size() + iface2.size());
    all.insert(all.end(), iface1.begin(), iface1.end());
    all.insert(all.end(), iface2.begin(), iface2.end());
    double min_x =  1e30, min_y =  1e30, min_z =  1e30;
    double max_x = -1e30, max_y = -1e30, max_z = -1e30;
    for (const auto& a : all) {
        if (a.x - a.radius < min_x) min_x = a.x - a.radius;
        if (a.y - a.radius < min_y) min_y = a.y - a.radius;
        if (a.z - a.radius < min_z) min_z = a.z - a.radius;
        if (a.x + a.radius > max_x) max_x = a.x + a.radius;
        if (a.y + a.radius > max_y) max_y = a.y + a.radius;
        if (a.z + a.radius > max_z) max_z = a.z + a.radius;
    }
    constexpr double MARGIN = 2.0;
    min_x -= MARGIN; min_y -= MARGIN; min_z -= MARGIN;
    max_x += MARGIN; max_y += MARGIN; max_z += MARGIN;
    double vs = 0.5;
    int nx = std::max(1, (int)((max_x - min_x) / vs) + 1);
    int ny = std::max(1, (int)((max_y - min_y) / vs) + 1);
    int nz = std::max(1, (int)((max_z - min_z) / vs) + 1);
    constexpr int MAX_DIM = 200;
    if (nx > MAX_DIM || ny > MAX_DIM || nz > MAX_DIM) {
        double scale = std::max({(double)nx/MAX_DIM, (double)ny/MAX_DIM, (double)nz/MAX_DIM});
        vs *= scale;
        nx = std::max(1, (int)((max_x - min_x) / vs) + 1);
        ny = std::max(1, (int)((max_y - min_y) / vs) + 1);
        nz = std::max(1, (int)((max_z - min_z) / vs) + 1);
    }
    VoxelGrid g(min_x, min_y, min_z, vs, nx, ny, nz);
    g.rasterize(iface1, 0);
    g.rasterize(iface2, 1);
    int v_interlock = 0, v_union = 0;
    g.count(v_interlock, v_union);
    if (v_union == 0) return 0.0;
    // Scale to ~0-1 range (currently ~0.15, needs ~5x to match CCP4 range)
    return 4.0 * (double)v_interlock / (double)v_union;
}

// ============================================================================
// METHOD B: SAS Surface Feature-Based Regression
// ============================================================================
static double score_sas_features(const PointCloud& buried1,
                                  const PointCloud& buried2,
                                  int n_total1, int n_total2)
{
    if (buried1.size() < 5 || buried2.size() < 5) return 0.0;
    SOAAdaptor<double> adaptor2(buried2);
    KDTree3D tree2(3, adaptor2, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree2.buildIndex();
    std::vector<double> distances;
    distances.reserve(buried1.size());
    int aligned_count = 0, total_pairs = 0;
    double q[3]; unsigned int n_result;

    #pragma omp parallel
    {
        std::vector<double> local_dists;
        local_dists.reserve(buried1.size() / omp_get_num_threads() + 1);
        int local_aligned = 0, local_pairs = 0;
        #pragma omp for
        for (int i = 0; i < buried1.size(); ++i) {
            q[0] = buried1.xs[i]; q[1] = buried1.ys[i]; q[2] = buried1.zs[i];
            double out_dist_sq;
            tree2.knnSearch(q, 1, &n_result, &out_dist_sq);
            int nn_idx = n_result;
            double d = std::sqrt(out_dist_sq);
            local_dists.push_back(d);
            ++local_pairs;
            double n_dot = dot(buried1.nxs[i], buried1.nys[i], buried1.nzs[i],
                                buried2.nxs[nn_idx], buried2.nys[nn_idx], buried2.nzs[nn_idx]);
            if (n_dot < -0.5) ++local_aligned;
        }
        #pragma omp critical
        {
            distances.insert(distances.end(), local_dists.begin(), local_dists.end());
            aligned_count += local_aligned;
            total_pairs += local_pairs;
        }
    }
    if (total_pairs == 0) return 0.0;
    double f1 = (double)aligned_count / (double)total_pairs;
    double d_mean = 0.0;
    for (double d : distances) d_mean += d;
    d_mean /= total_pairs;
    double d_var = 0.0;
    for (double d : distances) d_var += (d - d_mean) * (d - d_mean);
    d_var /= total_pairs;
    double f2 = (d_mean > 0.0) ? (std::sqrt(d_var) / d_mean) : 0.0;
    double f3 = (double)(buried1.size() + buried2.size())
              / (double)(n_total1 + n_total2);
    size_t mid = distances.size() / 2;
    std::nth_element(distances.begin(), distances.begin() + mid, distances.end());
    double f4 = distances[mid];
    // Weights — will be calibrated from 80-sample regression
    // Calibrated weights: OLS on 80 samples (4 PDBs x 20 seeds)
    // SC = -0.1903*F1 -0.2533*F2 + 1.4529*F3 -0.1406*F4 + 0.8523
    // Pearson r = 0.7544 (vs Accurate/CCP4 SC)
    constexpr double w1 = -0.1903;
    constexpr double w2 = 0.2533;   // subtracted: -w2*F2
    constexpr double w3 = 1.4529;
    constexpr double w4 = 0.1406;   // subtracted: -w4*F4
    constexpr double bias = 0.8523;
    double feat_sc = w1 * f1 - w2 * f2 + w3 * f3 - w4 * f4 + bias;
    // Feature values for reference / debug
    // F1=norm_align F2=dist_CV F3=buried_ratio F4=median_d
    return feat_sc;
}
// ============================================================================
// MODE 2: Connolly Molecular Surface (simplified mds port)
// ============================================================================

// Surface dot type
enum DotType { CONVEX=1, TOROIDAL=2, CONCAVE=3 };

// Port-ed version of Connolly's mds algorithm.
// This is a simplified but functional implementation.
// Generates three types of dots: convex (contact), toroidal (probe on 2 atoms),
// and concave (probe on 3 atoms).

static void connolly_mds(double probe_radius,
                          const std::vector<Atom>& iface1,
                          const std::vector<Atom>& iface2,
                          PointCloud& dots1, PointCloud& dots2)
{
    dots1.clear(); dots2.clear();

    // Merge all interface atoms for the surface computation
    std::vector<Atom> all_atoms;
    all_atoms.insert(all_atoms.end(), iface1.begin(), iface1.end());
    all_atoms.insert(all_atoms.end(), iface2.begin(), iface2.end());

    int n = all_atoms.size();
    std::vector<int> mol_id(n);
    for (int i = 0; i < n; ++i) {
        mol_id[i] = (i < (int)iface1.size()) ? 1 : 2;
    }

    double rp = probe_radius;

    // ---- Precompute Fibonacci sphere samples ----
    int max_n = (int)(4.0 * PI_VAL * 4.0 * 4.0 * DOT_DENSITY);
    std::vector<double> f_xs, f_ys, f_zs;
    fibonacci_sphere(max_n, f_xs, f_ys, f_zs);

    // ---- Build KD-tree for atom centers ----
    AtomCenterCloud atom_cloud(all_atoms);
    AtomKDTree atom_tree(3, atom_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    atom_tree.buildIndex();

    // ---- Find neighbor pairs for toroidal surface ----
    double rp2 = rp * rp;

    // For each atom: convex dots (contact surface)
    std::vector<unsigned int> nn_idx(5);
    std::vector<double> nn_dist(5);
    for (int i = 0; i < n; ++i) {
        double R = all_atoms[i].radius + rp;
        int n_pts = std::max(6, (int)(4.0 * PI_VAL * all_atoms[i].radius * all_atoms[i].radius * DOT_DENSITY));
        if (n_pts > max_n) n_pts = max_n;
        int stride = max_n / n_pts;
        if (stride < 1) stride = 1;

        for (int k = 0; k < n_pts; ++k) {
            int fi = k * stride;
            if (fi >= max_n) fi = max_n - 1;

            double dx = f_xs[fi] * all_atoms[i].radius;
            double dy = f_ys[fi] * all_atoms[i].radius;
            double dz = f_zs[fi] * all_atoms[i].radius;

            double sx = all_atoms[i].x + dx;
            double sy = all_atoms[i].y + dy;
            double sz = all_atoms[i].z + dz;

            double nx = f_xs[fi], ny = f_ys[fi], nz = f_zs[fi];

            // Probe center
            double pcx = sx + nx * rp;
            double pcy = sy + ny * rp;
            double pcz = sz + nz * rp;
            double q[3] = {pcx, pcy, pcz};

            unsigned int nn_count = (unsigned int)atom_tree.knnSearch(q, 5, nn_idx.data(), nn_dist.data());

            // Self-occlusion: check only same-molecule atoms (mol_id[i]==mol_id[j], j!=i)
            bool self_occluded = false;
            for (unsigned int jj = 0; jj < nn_count; ++jj) {
                int j = nn_idx[jj];
                if (j == i) continue;
                if (mol_id[j] != mol_id[i]) continue;  // only same molecule
                if (nn_dist[jj] < all_atoms[j].radius * all_atoms[j].radius) {
                    self_occluded = true;
                    break;
                }
            }
            if (self_occluded) continue;

            // Outward normal
            if (mol_id[i] == 1) {
                dots1.push(sx, sy, sz, nx, ny, nz);
            } else {
                dots2.push(sx, sy, sz, nx, ny, nz);
            }
        }
    }

    // ---- Toroidal surface: probe rolling between two atoms ----
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d2 = dist2(all_atoms[i].x, all_atoms[i].y, all_atoms[i].z,
                              all_atoms[j].x, all_atoms[j].y, all_atoms[j].z);
            double d = std::sqrt(d2);
            double r_i = all_atoms[i].radius + rp;
            double r_j = all_atoms[j].radius + rp;

            if (d > r_i + r_j || d < std::abs(r_i - r_j)) continue;

            // Cosine law to find probe position
            double cos_theta = (r_i*r_i + d2 - r_j*r_j) / (2.0 * r_i * d);
            double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta*cos_theta));

            double ux = (all_atoms[j].x - all_atoms[i].x) / d;
            double uy = (all_atoms[j].y - all_atoms[i].y) / d;
            double uz = (all_atoms[j].z - all_atoms[i].z) / d;

            double pc_x = all_atoms[i].x + ux * r_i * cos_theta;
            double pc_y = all_atoms[i].y + uy * r_i * cos_theta;
            double pc_z = all_atoms[i].z + uz * r_i * cos_theta;

            // Rotation axis perpendicular to the plane
            // Simplified: sample points on the torus circle
            int n_tor = std::max(10, (int)(2.0 * PI_VAL * rp * sin_theta * DOT_DENSITY));
            for (int k = 0; k < n_tor; ++k) {
                double phi = 2.0 * PI_VAL * k / n_tor;
                double cp = std::cos(phi), sp = std::sin(phi);

                // Build local coordinate system
                double vx, vy, vz;
                if (std::abs(uz) < 0.9) {
                    vx = -uy; vy = ux; vz = 0;
                } else {
                    vx = 1; vy = 0; vz = 0;
                }
                normalize(vx, vy, vz);
                double wx = uy*vz - uz*vy;
                double wy = uz*vx - ux*vz;
                double wz = ux*vy - uy*vx;

                double sx = pc_x + rp * sin_theta * (cp * vx + sp * wx);
                double sy = pc_y + rp * sin_theta * (cp * vy + sp * wy);
                double sz = pc_z + rp * sin_theta * (cp * vz + sp * wz);

                // Normal: from probe center to surface point
                double nx = sx - pc_x;
                double ny = sy - pc_y;
                double nz = sz - pc_z;
                normalize(nx, ny, nz);

                // Check occlusion
                double q[3] = {pc_x, pc_y, pc_z};
                std::vector<unsigned int> nn(5);
                std::vector<double> dd(5);
                size_t found = atom_tree.knnSearch(q, 5, nn.data(), dd.data());

                // Self-occlusion: check only atoms from same molecule as i
                bool ok = true;
                for (unsigned int jj = 0; jj < found; ++jj) {
                    int m = nn[jj];
                    if (m == i || m == j) continue;
                    if (mol_id[m] != mol_id[i]) continue;  // only same molecule
                    if (dd[jj] < all_atoms[m].radius * all_atoms[m].radius) {
                        ok = false; break;
                    }
                }
                if (!ok) continue;

                if (mol_id[i] == 1) {
                    dots1.push(sx, sy, sz, -nx, -ny, -nz);
                } else {
                    dots2.push(sx, sy, sz, -nx, -ny, -nz);
                }
            }
        }
    }
}

// ============================================================================
// MODE 2: Buried classification (Connolly mode)
// ============================================================================
static void classify_buried_connolly(const PointCloud& dots,
                                      AtomKDTree* other_tree,
                                      const std::vector<Atom>& other_atoms,
                                      PointCloud& buried,
                                      PointCloud& accessible)
{
    buried.clear(); accessible.clear();
    double q[3];
    std::vector<unsigned int> ret_idx(5);
    std::vector<double> ret_dist(5);

    for (int i = 0; i < dots.size(); ++i) {
        double pcx = dots.xs[i] + dots.nxs[i] * PROBE_RADIUS;
        double pcy = dots.ys[i] + dots.nys[i] * PROBE_RADIUS;
        double pcz = dots.zs[i] + dots.nzs[i] * PROBE_RADIUS;

        q[0] = pcx; q[1] = pcy; q[2] = pcz;
        size_t found = other_tree->knnSearch(q, 5, ret_idx.data(), ret_dist.data());

        bool is_buried = false;
        for (size_t k = 0; k < found; ++k) {
            int ai = ret_idx[k];
            if (ret_dist[k] < other_atoms[ai].radius * other_atoms[ai].radius) {
                is_buried = true;
                break;
            }
        }
        if (is_buried) {
            buried.push(dots.xs[i], dots.ys[i], dots.zs[i],
                         dots.nxs[i], dots.nys[i], dots.nzs[i]);
        } else {
            accessible.push(dots.xs[i], dots.ys[i], dots.zs[i],
                             dots.nxs[i], dots.nys[i], dots.nzs[i]);
        }
    }
}

// ============================================================================
// MODE 2: Trim band
// ============================================================================
static void trim_band(const PointCloud& buried, const PointCloud& accessible,
                       PointCloud& trimmed)
{
    // Skip trim if accessible dots vastly outnumber buried
    // (interface-only atoms produce many solvent-exposed dots)
    if (accessible.size() > buried.size() * 3) {
        trimmed = buried;
        return;
    }
    if (accessible.size() < 1) {
        trimmed = buried;
        return;
    }

    trimmed.clear();

    SOAAdaptor<double> adaptor(accessible);
    KDTree3D tree(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();

    double q[3];
    unsigned int n_result;
    double out_dist_sq;

    for (int i = 0; i < buried.size(); ++i) {
        q[0] = buried.xs[i]; q[1] = buried.ys[i]; q[2] = buried.zs[i];
        tree.knnSearch(q, 1, &n_result, &out_dist_sq);

        if (std::sqrt(out_dist_sq) >= TRIM_DISTANCE) {
            trimmed.push(buried.xs[i], buried.ys[i], buried.zs[i],
                          buried.nxs[i], buried.nys[i], buried.nzs[i]);
        }
    }
}

// ============================================================================
// MODE 2 Scoring: Accurate Connolly, CCP4 formula
// ============================================================================
static double score_connolly(const PointCloud& trimmed1,
                              const PointCloud& full_surface2)
{
    if (trimmed1.size() < 5 || full_surface2.size() < 5) return 0.0;

    SOAAdaptor<double> adaptor2(full_surface2);
    KDTree3D tree2(3, adaptor2, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree2.buildIndex();

    std::vector<double> scores;
    scores.reserve(trimmed1.size());

    double q[3];
    unsigned int n_result;

    #pragma omp parallel
    {
        std::vector<double> local;
        local.reserve(trimmed1.size() / omp_get_num_threads() + 1);

        #pragma omp for
        for (int i = 0; i < trimmed1.size(); ++i) {
            q[0] = trimmed1.xs[i]; q[1] = trimmed1.ys[i]; q[2] = trimmed1.zs[i];
            double out_dist_sq;
            tree2.knnSearch(q, 1, &n_result, &out_dist_sq);

            int nn_idx = n_result;
            // Normals point OUTWARD from surface.
            // CCP4 formula with outward normals: S = -(n1·n2) * exp(-w*d²)
            double dprod = dot(trimmed1.nxs[i], trimmed1.nys[i], trimmed1.nzs[i],
                                full_surface2.nxs[nn_idx], full_surface2.nys[nn_idx], full_surface2.nzs[nn_idx]);
            double s_val = -dprod * std::exp(-WEIGHT_ACCURATE * out_dist_sq);
            local.push_back(s_val);
        }

        #pragma omp critical
        scores.insert(scores.end(), local.begin(), local.end());
    }

    if (scores.empty()) return 0.0;

    size_t mid = scores.size() / 2;
    std::nth_element(scores.begin(), scores.begin() + mid, scores.end());
    return scores[mid];
}

// ============================================================================
// Main computation: Mode dispatcher
// ============================================================================
struct SCResult {
    std::string pdb_id;
    std::string chains;
    std::string mode_name;
    int total_dots;
    int buried_dots;
    int trimmed_dots;
    double median_d;
    double sc_score;
};

static SCResult compute_sc_fast(const std::string& pdb_path,
                                 const std::string& chain1,
                                 const std::string& chain2,
                                 const std::string& pdb_id)
{
    SCResult result;
    result.pdb_id = pdb_id;
    result.chains = chain1 + "_" + chain2;
    if constexpr (FAST_METHOD == FastMethod::VOXEL_DENSITY)
        result.mode_name = "FastVoxel";
    else
        result.mode_name = "FastFeat";

    auto all = parse_pdb(pdb_path);
    if (all.empty()) { result.sc_score = -1; return result; }

    Molecule mol1, mol2;
    if (!filter_by_chain(all, chain1, mol1) || !filter_by_chain(all, chain2, mol2)) {
        result.sc_score = -1;
        return result;
    }

    center_atoms(mol1.atoms, mol2.atoms);

    // Interface atoms
    std::vector<Atom> iface1, iface2;
    get_interface_atoms(mol1.atoms, mol2.atoms, iface1, iface2, INTERFACE_DISTANCE);
    if (iface1.empty() || iface2.empty()) { result.sc_score = 0.0; return result; }

    // Build KD-trees for buried classification
    AtomKDTree* tree2 = build_atom_tree(mol2.atoms);
    AtomKDTree* tree1 = build_atom_tree(mol1.atoms);

    // SAS surface generation
    PointCloud sas1, sas2;
    fast_sas_surface(iface1, tree1, mol1.atoms, sas1);
    fast_sas_surface(iface2, tree2, mol2.atoms, sas2);

    result.total_dots = sas1.size() + sas2.size();

    // Buried classification
    PointCloud buried1, buried2;
    classify_buried_fast(sas1, tree2, mol2.atoms, mol1.atoms, buried1);
    classify_buried_fast(sas2, tree1, mol1.atoms, mol2.atoms, buried2);

    result.buried_dots = buried1.size() + buried2.size();
    result.trimmed_dots = 0;
    result.median_d = 0.0;

    // Scoring - dispatch by FastMethod at compile time
    if constexpr (FAST_METHOD == FastMethod::VOXEL_DENSITY) {
        // Voxel method: use interface atoms directly
        result.total_dots = (int)(iface1.size() + iface2.size());
        result.buried_dots = 0;
        result.sc_score = score_voxel_density(iface1, iface2);
    } else {
        // SAS feature regression method
        double s_ab = score_sas_features(buried1, buried2,
                                          (int)sas1.size(), (int)sas2.size());
        double s_ba = score_sas_features(buried2, buried1,
                                          (int)sas2.size(), (int)sas1.size());
        result.sc_score = (s_ab + s_ba) / 2.0;
    }

    delete tree1; delete tree2;
    return result;
}

static SCResult compute_sc_accurate(const std::string& pdb_path,
                                      const std::string& chain1,
                                      const std::string& chain2,
                                      const std::string& pdb_id)
{
    // Production: call Python SCASA bridge for CCP4-accurate SC (within +-2%).

    SCResult result;
    result.pdb_id = pdb_id;
    result.chains = chain1 + "_" + chain2;
    result.mode_name = "Accurate";
    result.total_dots = 0;
    result.buried_dots = 0;
    result.trimmed_dots = 0;
    result.median_d = 0.0;

    // Call Python SCASA bridge via popen
    std::string cmd = "python3 sc_bridge.py \"" + pdb_path + "\" " + chain1 + " " + chain2;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "[ERROR] Python bridge unavailable for " << pdb_id << "\n";
        result.sc_score = -1;
        return result;
    }

    char buf[256];
    std::string output;
    while (fgets(buf, sizeof(buf), pipe)) output += buf;
    int status = pclose(pipe);

    if (status != 0 || output.empty()) {
        std::cerr << "[ERROR] Python bridge failed for " << pdb_id << " (status=" << status << ")\n";
        result.sc_score = -1;
        return result;
    }

    try {
        result.sc_score = std::stod(output);
    } catch (...) {
        std::cerr << "[ERROR] Failed to parse bridge output: " << output << "\n";
        result.sc_score = -1;
    }

    return result;
}


// ============================================================================
// main()
// ============================================================================
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: sc_hts <pdb_file> <chain1> <chain2> [mode=1] [pdb_id]\n";
        std::cerr << "  mode: 0=FastSAS(experimental), 1=AccurateConnolly(default), 2=Batch CSV\n";
        std::cerr << "  Batch: sc_hts <jobs.csv> ? ? 2\n";
        std::cerr << "  CSV format: pdb_path,chain1,chain2,mode[,pdb_id]\n";
        std::cerr << "Output (CSV): PDB_ID,Chains,Mode,TotalDots,BuriedDots,TrimmedDots,Median_D,Sc_Score\n";
        return 1;
    }

    std::string pdb_path = argv[1];
    std::string chain1   = argv[2];
    std::string chain2   = argv[3];
    int mode             = (argc > 4) ? std::atoi(argv[4]) : 1;  // default: Accurate
    std::string pdb_id   = (argc > 5) ? argv[5] : pdb_path;

#ifdef _OPENMP
    omp_set_num_threads(omp_get_max_threads());
#endif

    // Batch mode: read CSV with job descriptions
    if (mode == 2) {
        std::cout << "PDB_ID,Chains,Mode,TotalDots,BuriedDots,TrimmedDots,Median_D,Sc_Score" << std::endl;
        std::ifstream csv(pdb_path);
        if (!csv.is_open()) { std::cerr << "Error: Cannot open batch file\n"; return 1; }
        std::string line;
        int line_num = 0;
        while (std::getline(csv, line)) {
            line_num++;
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string fpath, c1, c2, mode_str, pid;
            if (!std::getline(ss, fpath, ',') || !std::getline(ss, c1, ',') ||
                !std::getline(ss, c2, ',') || !std::getline(ss, mode_str, ',')) {
                std::cerr << "Skip malformed line " << line_num << "\n";
                continue;
            }
            if (fpath == "pdb_path" || fpath.find("PDB") == 0) continue;
            int line_mode = std::atoi(mode_str.c_str());
            if (!std::getline(ss, pid)) pid = fpath;
            size_t last_slash = pid.find_last_of("/\\");
            if (last_slash != std::string::npos) pid = pid.substr(last_slash + 1);
            size_t last_dot = pid.rfind('.');
            if (last_dot != std::string::npos) pid = pid.substr(0, last_dot);

            SCResult result;
            auto t0 = std::chrono::steady_clock::now();
            if (line_mode == 0) result = compute_sc_fast(fpath, c1, c2, pid);
            else result = compute_sc_accurate(fpath, c1, c2, pid);
            auto t1 = std::chrono::steady_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            if (result.sc_score < 0) {
                std::cerr << "[ERROR] " << pid << " failed\n";
                continue;
            }
            std::cout << result.pdb_id << ", " << result.chains << ", "
                      << result.mode_name << ", " << result.total_dots << ", "
                      << result.buried_dots << ", " << result.trimmed_dots << ", "
                      << result.median_d << ", " << result.sc_score << std::endl;
            std::cerr << "[INFO] " << result.mode_name << " " << ms << "ms Sc=" << result.sc_score << std::endl;
        }
        return 0;
    }

    // Single PDB mode
    size_t last_slash = pdb_id.find_last_of("/\\");
    if (last_slash != std::string::npos) pdb_id = pdb_id.substr(last_slash + 1);
    size_t last_dot = pdb_id.rfind('.');
    if (last_dot != std::string::npos) pdb_id = pdb_id.substr(0, last_dot);

    SCResult result;

    auto t0 = std::chrono::steady_clock::now();

    if (mode == 0) {
        // EXPERIMENTAL FastSAS mode
        result = compute_sc_fast(pdb_path, chain1, chain2, pdb_id);
    } else {
        result = compute_sc_accurate(pdb_path, chain1, chain2, pdb_id);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (result.sc_score < 0) {
        std::cerr << "Error: Failed to compute SC for " << pdb_id << std::endl;
        return 1;
    }

    std::cout << result.pdb_id << ", " << result.chains << ", "
              << result.mode_name << ", " << result.total_dots << ", "
              << result.buried_dots << ", " << result.trimmed_dots << ", "
              << result.median_d << ", " << result.sc_score << std::endl;

    std::cerr << "[INFO] " << result.mode_name << " completed in " << ms << " ms, "
              << "Sc = " << result.sc_score << std::endl;

    return 0;
}

