import csv
import sys
from argparse import Namespace

import pandas as pd

sys.path.insert(0, "/mnt/sdb2/lluoto/int-iScore/src")

from int_iscore.modes.af3_nonref import process_nonref_file
from int_iscore.modes.af3_reference import process_single_file
from int_iscore.modes.md_snapshots import process_md_snapshot


def as_dict(headers, row):
    return {header: value for header, value in zip(headers, row)}


def normalize_filename(value):
    if value.endswith(".cif") or value.endswith(".pdb"):
        return value
    return value + ".cif"


def main():
    rows = []

    nonref_args = Namespace(
        chimera_path=None,
        vmd_path=None,
        pdb2pqr_path=None,
        apbs_path=None,
        ccp4_path=None,
        temp_dir="/mnt/sdb2/lluoto/int-iScore/int_iscore_temp",
    )
    nonref_headers = [
        "filename", "dope", "soap", "sc", "ec", "cpscore",
        "bsa", "interface_residues", "pairs", "ramachandran",
    ]
    nonref_row = process_nonref_file(
        "/mnt/sdb2/lluoto/5-HT2C-2CD28/pred.rank_0.cif",
        (["A"], ["B", "C"]),
        nonref_args,
    )
    nonref_result = as_dict(nonref_headers, nonref_row)
    nonref_legacy = pd.read_csv("/mnt/sdb2/lluoto/int-iScore/5HT2C_intiscore.csv")
    nonref_legacy_row = nonref_legacy[nonref_legacy["filename"] == "pred.rank_0.cif"].iloc[0].to_dict()
    rows.append({
        "mode": "af3-nonref",
        "sample": "pred.rank_0.cif",
        "legacy_dope": nonref_legacy_row["dope"],
        "new_dope": nonref_result["dope"],
        "legacy_soap": nonref_legacy_row["soap"],
        "new_soap": nonref_result["soap"],
        "legacy_cpscore": nonref_legacy_row["cpscore"],
        "new_cpscore": nonref_result["cpscore"],
        "legacy_bsa": nonref_legacy_row["bsa"],
        "new_bsa": nonref_result["bsa"],
        "legacy_ec": nonref_legacy_row["ec"],
        "new_ec": nonref_result["ec"],
        "new_sc": nonref_result["sc"],
    })

    ref_args = Namespace(
        input="/mnt/sdb2/lluoto/6A6I",
        reference=None,
        reference_dir="/mnt/sdb2/lluoto/int-iScore",
        server_output=False,
        no_clean_pdb=True,
        chimera_path=None,
        vmd_path=None,
        pdb2pqr_path=None,
        apbs_path=None,
        ccp4_path=None,
    )
    ref_headers = [
        "filename", "average_plddt", "ranking_score", "iptm", "ptm", "dope", "soap",
        "dockq", "sc", "cpscore", "bsa", "interface_plddt", "notification",
        "interface_residues", "pairs", "ramachandran", "ec",
    ]
    ref_row = process_single_file(
        "/mnt/sdb2/lluoto/6A6I/seed-614460_sample-1/summary_confidences.json",
        (["A"], ["B"]),
        ref_args,
    )
    ref_result = as_dict(ref_headers, ref_row)
    ref_legacy = pd.read_csv("/mnt/sdb2/lluoto/int-iScore/all_3_6.csv")
    ref_legacy_row = ref_legacy[ref_legacy["filename"] == "6A6I_seed-614460_sample-1model"].iloc[0].to_dict()
    rows.append({
        "mode": "af3-ref",
        "sample": "6A6I_seed-614460_sample-1model",
        "legacy_dope": ref_legacy_row["dope"],
        "new_dope": ref_result["dope"],
        "legacy_soap": ref_legacy_row["soap"],
        "new_soap": ref_result["soap"],
        "legacy_cpscore": ref_legacy_row["cpscore"],
        "new_cpscore": ref_result["cpscore"],
        "legacy_bsa": ref_legacy_row["bsa"],
        "new_bsa": ref_result["bsa"],
        "legacy_ec": ref_legacy_row["ec"],
        "new_ec": ref_result["ec"],
        "legacy_sc": ref_legacy_row["sc"],
        "new_sc": ref_result["sc"],
        "legacy_dockq": ref_legacy_row["dockq"],
        "new_dockq": ref_result["dockq"],
    })

    md_args = Namespace(
        chimera_path=None,
        vmd_path=None,
        pdb2pqr_path=None,
        apbs_path=None,
        ccp4_path=None,
    )
    md_headers = [
        "filename", "dope", "soap", "frustration_score", "cpscore",
        "bsa", "interface_residues", "ramachandran", "pairs", "sc", "ec",
    ]
    md_row = process_md_snapshot(
        "/mnt/sdb2/lluoto/split_frames/3070_0082-repeat_1_0731.pdb",
        (["A", "B", "C", "D"], ["E"]),
        md_args,
    )
    md_result = as_dict(md_headers, md_row)
    with open("/mnt/sdb2/lluoto/int-iScore/md_3070_0082_hi_3_3.csv", newline="") as handle:
        reader = csv.reader(handle)
        first_md_row = next(reader)
    md_legacy_row = {
        "filename": first_md_row[0],
        "dope": float(first_md_row[1]),
        "soap": float(first_md_row[2]),
        "cpscore": float(first_md_row[4]),
        "bsa": float(first_md_row[5]),
    }
    rows.append({
        "mode": "md",
        "sample": "3070_0082-repeat_1_0731.pdb",
        "legacy_dope": md_legacy_row["dope"],
        "new_dope": md_result["dope"],
        "legacy_soap": md_legacy_row["soap"],
        "new_soap": md_result["soap"],
        "legacy_cpscore": md_legacy_row["cpscore"],
        "new_cpscore": md_result["cpscore"],
        "legacy_bsa": md_legacy_row["bsa"],
        "new_bsa": md_result["bsa"],
        "new_ec": md_result["ec"],
        "new_sc": md_result["sc"],
    })

    df = pd.DataFrame(rows)
    output = "/mnt/sdb2/lluoto/int-iScore/int_iscore_temp/substitution_test_results.csv"
    df.to_csv(output, index=False)
    print(df.to_string(index=False))
    print(f"WROTE {output}")


if __name__ == "__main__":
    main()
