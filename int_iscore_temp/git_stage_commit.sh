#!/bin/bash
set -euo pipefail

cd /media/cuixi/data01/lluoto/int-iScore

git add \
  README.md \
  CURRENT_STATUS_SC_BSA.md \
  pyproject.toml \
  setup.py \
  src/int_iscore/utils/sc_calculator.py \
  src/int_iscore/utils/sc_backend.cpp \
  src/int_iscore/utils/sc_ec.py \
  src/int_iscore/core/metrics.py \
  src/int_iscore/tools/__init__.py \
  src/int_iscore/tools/frustration_legacy.py \
  src/int_iscore/tools/input_data_generate_legacy.py \
  src/int_iscore/modes/md_snapshots.py \
  src/int_iscore/modes/af3_nonref.py \
  src/int_iscore/modes/af3_reference.py \
  src/int_iscore/cli.py \
  calculate_md.py \
  calculate_non_ref_af3.py \
  calculate_af3.py

git diff --cached --numstat
git commit -m "integrate legacy wrappers and document current SC/BSA status"
