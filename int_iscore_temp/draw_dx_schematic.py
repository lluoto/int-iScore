import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path('/media/cuixi/data01/lluoto/int-iScore')
PQR = ROOT / 'test.pqr'
DX = ROOT / 'test.pqr.dx'
OUT = ROOT / 'int_iscore_temp' / 'test_dx_schematic.png'


def parse_pqr(path: Path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            name = line[12:16].strip()
            if name.startswith('H'):
                continue
            atoms.append((
                float(line[30:38].strip()),
                float(line[38:46].strip()),
                float(line[46:54].strip()),
                name,
            ))
    return np.array([(x, y, z) for x, y, z, _ in atoms], dtype=float), [n for _, _, _, n in atoms]


def parse_dx(path: Path):
    counts = None
    origin = None
    deltas = []
    values = []
    reading = False
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('object 1 class gridpositions counts'):
                counts = tuple(map(int, line.split()[-3:]))
            elif line.startswith('origin'):
                origin = np.array(list(map(float, line.split()[1:4])), dtype=float)
            elif line.startswith('delta'):
                deltas.append(np.array(list(map(float, line.split()[1:4])), dtype=float))
            elif 'data follows' in line:
                reading = True
            elif reading:
                if line.startswith('attribute') or line.startswith('object'):
                    break
                values.extend(map(float, line.split()))
    if counts is None or origin is None or len(deltas) != 3:
        raise RuntimeError('Failed to parse DX header')
    grid = np.array(values, dtype=float).reshape(counts, order='C')
    spacing = np.array([d[np.flatnonzero(np.abs(d) > 1e-12)[0]] for d in deltas], dtype=float)
    return grid, origin, spacing


def slice_extent(origin, spacing, shape, a, b):
    xmin = origin[a]
    xmax = origin[a] + spacing[a] * (shape[a] - 1)
    ymin = origin[b]
    ymax = origin[b] + spacing[b] * (shape[b] - 1)
    return [xmin, xmax, ymin, ymax]


coords, names = parse_pqr(PQR)
grid, origin, spacing = parse_dx(DX)

cx, cy, cz = [s // 2 for s in grid.shape]
slices = [
    ('XY @ mid-Z', grid[:, :, cz].T, (0, 1), coords[:, [0, 1]]),
    ('XZ @ mid-Y', grid[:, cy, :].T, (0, 2), coords[:, [0, 2]]),
    ('YZ @ mid-X', grid[cx, :, :].T, (1, 2), coords[:, [1, 2]]),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
vmax = np.percentile(np.abs(grid), 99)
vmin = -vmax

for ax, (title, data, dims, pts) in zip(axes, slices):
    extent = slice_extent(origin, spacing, grid.shape, dims[0], dims[1])
    im = ax.imshow(data, origin='lower', extent=extent, cmap='coolwarm', vmin=vmin, vmax=vmax, aspect='auto')
    ax.scatter(pts[:, 0], pts[:, 1], s=4, c='k', alpha=0.45)
    ax.set_title(title)
    ax.set_xlabel(['X', 'X', 'Y'][axes.tolist().index(ax)])
    ax.set_ylabel(['Y', 'Z', 'Z'][axes.tolist().index(ax)])

cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85)
cbar.set_label('Electrostatic potential (DX)')
fig.suptitle('Schematic DX grid slices with projected structure points')

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200)
print(str(OUT))
