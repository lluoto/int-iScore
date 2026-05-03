#!/bin/bash
set -euo pipefail

PDB_PATH="$1"
CHAIN1="$2"
shift 2
CHAINS2=("$@")

export CCP4_MASTER=/media/cuixi/data01/Softwares/ccp4-9.0.015-linux64
export CCP4=$CCP4_MASTER/ccp4-9
export CCP4_SCR=/tmp/ccp4_sc_$(whoami)
export CBIN=$CCP4/bin
export CLIB=$CCP4/lib
export CLIBD=$CCP4/lib/data
export CETC=$CCP4/etc
export CINCL=$CCP4/include
export CHTML=$CCP4/html
export CEXAM=$CCP4/examples
export CCP4I_TOP=$CCP4/share/ccp4i
export MMCIFDIC=$CLIB/ccp4/cif_mmdic.lib
export CLIBD_MON=$CCP4/lib/data/monomers/
export CCP4_HELPDIR=$CCP4/help/
export CCP4_OPEN=UNKNOWN
export PATH=$CCP4/etc:$CCP4/bin:$PATH

mkdir -p "$CCP4_SCR"

{
  printf 'MOLECULE 1\n'
  printf 'CHAIN %s\n' "$CHAIN1"
  printf 'MOLECULE 2\n'
  for chain in "${CHAINS2[@]}"; do
    printf 'CHAIN %s\n' "$chain"
  done
  printf 'END\n'
} | sc XYZIN "$PDB_PATH"
