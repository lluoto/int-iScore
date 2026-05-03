

import os

from sklearn.cluster import DBSCAN,AffinityPropagation
from mpi4py import MPI
import subprocess
from modeller import *
from modeller import soap_pp
import math
import itertools
import numpy as np
import csv
import ast
import subprocess
import time
from tqdm import tqdm
from pathlib import Path
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import warnings
warnings.filterwarnings('ignore')

# ================================================
# 1. First, diagnose your RMSD matrix
# ================================================
# def diagnose_rmsd_matrix(rmsd_matrix):
#     """Check if RMSD matrix has issues."""
    
#     # Check basic properties
#     print("=" * 60)
#     print("RMSD MATRIX DIAGNOSIS")
#     print("=" * 60)
    
#     print(f"Shape: {rmsd_matrix.shape}")
#     print(f"Diagonal all zeros: {np.all(np.diag(rmsd_matrix) == 0)}")
#     print(f"Symmetry: {np.allclose(rmsd_matrix, rmsd_matrix.T)}")
    
#     # Get off-diagonal values
#     mask = ~np.eye(rmsd_matrix.shape[0], dtype=bool)
#     off_diag = rmsd_matrix[mask]
    
#     print(f"\nRMSD Statistics:")
#     print(f"  Min: {off_diag.min():.4f}")
#     print(f"  Max: {off_diag.max():.4f}")
#     print(f"  Mean: {off_diag.mean():.4f}")
#     print(f"  Std: {off_diag.std():.4f}")
#     print(f"  Median: {np.median(off_diag):.4f}")
    
#     # Check for uniform distances
#     cv = off_diag.std() / off_diag.mean() if off_diag.mean() > 0 else 0
#     print(f"  Coefficient of Variation: {cv:.4f}")
    
#     if cv < 0.1:
#         print("⚠️  WARNING: Very low variation! All structures are similar.")
#     if off_diag.max() < 1.0:
#         print("⚠️  WARNING: Very small RMSD values (<1.0 Å)")
    
#     # Check distance distribution
#     plt.figure(figsize=(10, 4))
    
#     plt.subplot(1, 2, 1)
#     plt.hist(off_diag.flatten(), bins=30, edgecolor='black', alpha=0.7)
#     plt.xlabel('RMSD (Å)')
#     plt.ylabel('Frequency')
#     plt.title('RMSD Distribution')
#     plt.grid(True, alpha=0.3)
    
#     plt.subplot(1, 2, 2)
#     plt.boxplot(off_diag.flatten())
#     plt.title('RMSD Box Plot')
#     plt.ylabel('RMSD (Å)')
    
#     plt.tight_layout()
#     plt.show()
#     plt.figure(figsize=(8, 8))  # 正方形画布

# # 使用 'viridis' 或 'plasma' colormap。对于RMSD， sequential colormap 更合适。
#     im = plt.imshow(rmsd_matrix,
#                 cmap='plasma',  # 另一种选择：'viridis', 'inferno', 'magma', 'cividis'
#                 interpolation='nearest',  # 或 'none'
#                 origin='upper',
#                 vmin=0,  # 设置颜色范围的最小值
#                 vmax=np.max(rmsd_matrix)  # 设置颜色范围的最大值
#                 )

# # 添加颜色条
#     cbar = plt.colorbar(im, shrink=0.8, pad=0.02)
#     cbar.set_label('RMSD (Å)', fontsize=12)

#     # 标签和标题
#     plt.title('RMSD Matrix Heatmap', fontsize=14)
#     # plt.xlabel('Frame ', fontsize=12)
#     # plt.ylabel('Frame Index', fontsize=12)

#     # 关键步骤：强制热图显示为正方形
#     ax = plt.gca()
#     ax.set_aspect('equal')  # 这确保了热图的单元格是正方形的

#     # 可选：添加网格线或调整刻度
#     # plt.xticks(np.arange(0, n_frames, step=10))
#     # plt.yticks(np.arange(0, n_frames, step=10))

#     plt.tight_layout()
#     plt.savefig('rmsd_heatmap_square.png', dpi=300, bbox_inches='tight')
#     plt.show()

#     print(f"RMSD矩阵形状: {rmsd_matrix.shape}")
#     print(f"RMSD最小值: {np.min(rmsd_matrix):.2f} Å, 最大值: {np.max(rmsd_matrix):.2f} Å")
#     return off_diag
# rmsd_matrix = np.load('rmsd_matrix.npy')
# # Run diagnosis on your RMSD matrix
# diagnose_rmsd_matrix(rmsd_matrix)
import os
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from modeller import Environ, Model, Selection

all_ref = pd.read_csv("new_null_process_3_3.csv")


def calculate_dope(model_path):
    """
    对单个模型按链分别计算 normalized DOPE-HR，并返回各链得分之和。
    """
    env = Environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')
    model_path = os.path.join('/media/cuixi/data1/af3_benchmark', model_path+'.cif')
    mdl = Model(env, file=model_path)
    chains = mdl.chains

    dope_scores = []
    pid = os.getpid()

    for chain_id in range(len(chains)):
        temp_filename = f"{model_path}_pid{pid}_chain_{chain_id}.pdb"

        try:
            sel = Selection(mdl.chains[chain_id])
            sel.write(file=temp_filename)

            single_model = Model(env, file=temp_filename)
            dope = single_model.assess_normalized_dopehr()
            dope_scores.append(dope)

        finally:
            if os.path.isfile(temp_filename):
                try:
                    os.remove(temp_filename)
                except OSError:
                    pass

    return np.sum(dope_scores)


def main():
    model_list = all_ref["filename"].tolist()
    # n_proc = max(1, cpu_count() - 1)

    with Pool(processes=6) as pool:
        refine_dope_list = list(
            tqdm(
                pool.imap(calculate_dope, model_list),
                total=len(model_list),
                desc="Calculating DOPE"
            )
        )

    all_ref["refine_dope"] = refine_dope_list
    all_ref.to_csv("null_new.csv", index=False)


if __name__ == "__main__":
    main()