"""
Sequence trimming module.

Provides three methods for trimming protein sequences to their structured core:
1. PDB CA atom coordinates (actual resolved residues)
2. UniProt TRANSMEM annotations (curated TM boundaries)
3. Kyte-Doolittle hydropathy (fallback)
"""
import math
import json
import os

# Kyte-Doolittle hydropathy scale
KD_SCALE = {
    'A': 1.800, 'R': -4.500, 'N': -3.500, 'D': -3.500,
    'C': 2.500, 'Q': -3.500, 'E': -3.500, 'G': -0.400,
    'H': -3.200, 'I': 4.500, 'L': 3.800, 'K': -3.900,
    'M': 1.900, 'F': 2.800, 'P': -1.600, 'S': -0.800,
    'T': -0.700, 'W': -0.900, 'Y': -1.300, 'V': 4.200,
}

# Default parameters
WINDOW_SIZE = 19
TM_THRESHOLD_INIT = 1.2
TM_THRESHOLD_MIN = 0.6
MIN_TM_LENGTH = 12
MAX_TM_LENGTH = 42
EXPECTED_TMS = 7
TM_ZSCORE_THRESHOLD = 0.3


# ---------------------------------------------------------------------------
# Kyte-Doolittle hydropathy
# ---------------------------------------------------------------------------
def hydropathy_score(residue):
    return KD_SCALE.get(residue.upper(), 0.0)


def sliding_hydropathy(seq, window=WINDOW_SIZE):
    """Smoothed hydropathy scores via sliding window."""
    scores = []
    half = window // 2
    for i in range(len(seq)):
        ws = max(0, i - half)
        we = min(len(seq), i + half + 1)
        wseq = seq[ws:we]
        scores.append(sum(hydropathy_score(aa) for aa in wseq) / len(wseq) if wseq else 0.0)
    return scores


def compute_zscore(values):
    mean = sum(values) / len(values)
    var = sum((v - mean) ** 2 for v in values) / len(values)
    std = math.sqrt(var) if var > 0 else 1.0
    return [(v - mean) / std for v in values]


def predict_tm_regions(seq, min_len=MIN_TM_LENGTH, max_len=MAX_TM_LENGTH):
    """
    Predict TM regions using adaptive hydropathy + z-score.
    Returns (regions_list, scores).
    """
    scores = sliding_hydropathy(seq)
    z_scores = compute_zscore(scores)

    # Strategy 1: multiple absolute thresholds
    for thresh in [TM_THRESHOLD_INIT, 1.0, 0.8, TM_THRESHOLD_MIN]:
        regions = []
        i = 0
        binary = [1 if s >= thresh else 0 for s in scores]
        while i < len(binary):
            if binary[i]:
                start = i
                while i < len(binary) and binary[i]:
                    i += 1
                end = i - 1
                length = end - start + 1
                if min_len <= length <= max_len:
                    avg_hp = sum(scores[start:end+1]) / length
                    regions.append((start, end, avg_hp, length))
            else:
                i += 1
        if len(regions) >= EXPECTED_TMS:
            regions.sort(key=lambda x: x[2], reverse=True)
            tm = regions[:EXPECTED_TMS]
            tm.sort(key=lambda x: x[0])
            return tm, scores

    # Strategy 2: z-score based
    binary = [1 if z >= TM_ZSCORE_THRESHOLD else 0 for z in z_scores]
    i, all_reg = 0, []
    while i < len(binary):
        if binary[i]:
            start = i
            while i < len(binary) and binary[i]:
                i += 1
            end = i - 1
            length = end - start + 1
            if min_len <= length <= max_len:
                avg_z = sum(z_scores[start:end+1]) / length
                avg_hp = sum(scores[start:end+1]) / length
                all_reg.append((start, end, avg_hp, length, avg_z))
        else:
            i += 1
    if all_reg:
        all_reg.sort(key=lambda x: x[4], reverse=True)
        tm = all_reg[:EXPECTED_TMS]
        tm.sort(key=lambda x: x[0])
        return [(s, e, hp, ln) for s, e, hp, ln, _ in tm], scores

    # Strategy 3: lowest threshold fallback
    binary = [1 if s >= TM_THRESHOLD_MIN else 0 for s in scores]
    i, regions = 0, []
    while i < len(binary):
        if binary[i]:
            start = i
            while i < len(binary) and binary[i]:
                i += 1
            end = i - 1
            length = end - start + 1
            if min_len <= length <= max_len:
                regions.append((start, end, sum(scores[start:end+1])/length, length))
        else:
            i += 1
    regions.sort(key=lambda x: x[0])
    return regions, scores


def get_7tm_boundaries(seq):
    """
    Find the GPCR 7TM bundle cluster.
    Returns (start, end, tm_regions) 0-based inclusive.
    """
    tm_regions, scores = predict_tm_regions(seq)
    if len(tm_regions) < 3:
        return 0, len(seq) - 1, tm_regions

    MAX_SPAN = 500
    best_start = best_count = 0
    for i in range(len(tm_regions)):
        for j in range(i, len(tm_regions)):
            span = tm_regions[j][1] - tm_regions[i][0]
            count = j - i + 1
            if span <= MAX_SPAN and count >= best_count:
                best_start = i
                best_count = count

    if best_count >= 4:
        cluster = tm_regions[best_start:best_start + best_count]
        tm_start = cluster[0][0]
        tm_end = cluster[-1][1]
    else:
        tm_start = tm_regions[0][0]
        tm_end = tm_regions[-1][1]

    buf_n = max(0, tm_start - 20)
    buf_c = min(len(seq) - 1, tm_end + 20)
    return buf_n, buf_c, tm_regions


# ---------------------------------------------------------------------------
# PDB CA data
# ---------------------------------------------------------------------------
def load_pdb_ca_ranges(path=None):
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def load_uniprot_tm_data(path=None):
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


# ---------------------------------------------------------------------------
# Main trimming
# ---------------------------------------------------------------------------
def trim_sequence(seq, uniprot_id, pdb_ca_path=None, uniprot_tm_path=None):
    """
    Trim a sequence to its 7TM core domain.

    Priority:
      1. UniProt TRANSMEM annotations (curated)
      2. PDB CA atom coverage (validated against UniProt)
      3. Kyte-Doolittle hydropathy (fallback)

    Returns (trimmed_seq, start, end, method_string).
    """
    # 1. UniProt boundaries
    tm_data = load_uniprot_tm_data(uniprot_tm_path)
    up_info = tm_data.get(uniprot_id, {})
    up_count = up_info.get('tm_count', 0)
    up_start = up_info.get('trim_start')
    up_end = up_info.get('trim_end')
    up_len = up_info.get('trim_length', 0)

    if up_count >= 7 and up_start and up_end:
        uni_s = up_start - 1
        uni_e = up_end
        trimmed = seq[uni_s:uni_e]
        method = 'UniProt ({}-{})'.format(up_start, up_end)

        # 2. Check PDB CA consistency
        ca_data = load_pdb_ca_ranges(pdb_ca_path)
        ca_info = ca_data.get(uniprot_id, {})
        first_ca = ca_info.get('first_ca')
        last_ca = ca_info.get('last_ca')
        if first_ca and last_ca and last_ca <= len(seq):
            ca_len = last_ca - first_ca + 1
            diff = abs(ca_len - up_len)
            if diff <= 50:
                p_s = first_ca - 1
                p_e = last_ca - 1
                trimmed = seq[p_s:p_e + 1]
                pdb_id = ca_info.get('pdb_id', '?')
                method = 'CA:{}_{} ({}-{})'.format(pdb_id, ca_info.get('chain', '?'), first_ca, last_ca)
                return trimmed, p_s, p_e, method

        return trimmed, uni_s, uni_e - 1, method

    # 3. Kyte-Doolittle fallback
    start, end, tm_regions = get_7tm_boundaries(seq)
    trimmed = seq[start:end + 1]
    method = 'Kyte-Doolittle (TMs={})'.format(len(tm_regions))
    return trimmed, start, end, method
