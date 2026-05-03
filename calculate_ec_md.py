import sys
import numpy as np
import ast
import argparse
import csv
import os
import subprocess
import tempfile
import pandas as pd
# --- Function Definitions ---
def get_args():
    parser = argparse.ArgumentParser(description='Calculate electrostatic complementarity')
    parser.add_argument('--model','-m', type=str, default='', help='PQR file from APBS')
    parser.add_argument('--sc','-q', type=str, default='', help='SC score file')
    parser.add_argument('--output','-o', type=str, default='', help='Output directory')
    parser.add_argument('--input','-i', type=str, default='', help='Input directory')
    parser.add_argument('--vmd_path','-v', type=str, default='vmd', help='Path to VMD executable')
    return parser.parse_args()

def calculate_residue_scores(residue_potentials_dict):
    """Calculate electrostatic complementarity scores for residues"""
    residue_score = {}
    for residue, potentials in residue_potentials_dict.items():
        if potentials:  # Only calculate if potentials exist
            residue_score[residue] = np.mean(potentials)
    return residue_score


def map_atoms_vmd(pqr_file, interface_residues, dx_file, vmd_exec='vmd'):
    """Map electrostatic potential values to atoms/residues using VMD"""
    # Generate TCL script content with stderr redirection
    tcl_script_content = """
package require volutil

set M_PI [expr {4 * atan(1)}]
proc generate_sphere_points {center radius npoints} {
    global M_PI
    set points {}
    set phi [expr {$M_PI * (3.0 - sqrt(5.0))}]
    for {set i 0} {$i < $npoints} {incr i} {
        set y [expr {1.0 - (double($i) / ($npoints-1)) * 2.0}]
        set r [expr {sqrt(1 - $y*$y)}]
        set theta [expr {$phi * $i}]
        set x [expr {cos($theta) * $r}]
        set z [expr {sin($theta) * $r}]
        set px [expr {$x * $radius + [lindex $center 0]}]
        set py [expr {$y * $radius + [lindex $center 1]}]
        set pz [expr {$z * $radius + [lindex $center 2]}]
        lappend points [list $px $py $pz]
    }
    return $points
}

set pqr_file [lindex $argv 0]
set dx_file [lindex $argv 1]

mol new $pqr_file type pqr
set vol [volutil read dx $dx_file]
set sel [atomselect top 'not name H*']
set all_data {}

foreach index [$sel list] {
    set atom [atomselect top "index $index"]
    set res_id [lindex [$atom get resid] 0]
    set chain [lindex [$atom get chain] 0]
    set atom_name [lindex [$atom get name] 0]
    set res_name [lindex [$atom get resname] 0]
    set xyz [lindex [$atom get {x y z}] 0]
    set radius [lindex [$atom get radius] 0]
    
    if {$radius == ""} { 
        puts stderr "Warning: Missing radius for atom $atom_name in residue $res_id$chain"
        continue 
    }
    
    set shell_radius [expr {$radius + 1.4}]
    set points [generate_sphere_points $xyz $shell_radius 100]
    set potentials {}
    
    foreach point $points {
        if {[catch {volutil interpolate $vol $point} pot]} {
            puts stderr "Error interpolating at $point: $pot"
            continue
        }
        if {![string is double $pot]} {
            puts stderr "Warning: Invalid potential value '$pot' at point: $point"
            continue
        }
        lappend potentials $pot
    }
    
    if {[llength $potentials] == 0} { 
        puts stderr "Warning: No potentials for atom $atom_name in residue $res_id$chain"
        continue 
    }
    
    set avg_pot [expr [tcl::mathop::+ {*}$potentials] / [llength $potentials]]
    puts "[format "%d_%s" $res_id $chain] [format "%.4f" $avg_pot]"
}
exit
"""

    # Write TCL script to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tcl', delete=False) as f_tcl:
        f_tcl.write(tcl_script_content)
        tcl_script_path = f_tcl.name
    with open('tcl_text.tcl', 'w') as f:
        f.write(tcl_script_content)
    # Prepare residue dictionary
    residue_potentials_dict = {res: [] for res in interface_residues}
    
    # Run VMD command with text-only rendering
    cmd = [vmd_exec, '-dispdev', 'text', '-eofexit', '-e', tcl_script_path, 
           '-args', pqr_file, dx_file]
    print(' '.join(cmd))
    try:
        # Enable software rendering fallback
        env = os.environ.copy()
        env['VMD_USE_SOFTWARE_RENDERER'] = '1'
        os.system('LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/libGL.so.1:$LD_LIBRARY_PATH')
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env
        )
        stdout, stderr = process.communicate()
        
        # Log any errors from stderr
        if stderr:
            print("VMD stderr output:", file=sys.stderr)
            print(stderr, file=sys.stderr)
        
        # Parse stdout for valid numeric output
        for line in stdout.splitlines():
            parts = line.strip().split()
            if len(parts) != 2:
                continue
                
            residue_key, potential_str = parts
            try:
                potential = float(potential_str)
            except ValueError:
                continue
                
            if residue_key in residue_potentials_dict:
                residue_potentials_dict[residue_key].append(potential)
                
    finally:
        # Cleanup temp file
        os.unlink(tcl_script_path)
    
    return residue_potentials_dict
# --- Main Execution ---
if __name__ == "__main__":
    args = get_args()
    pqr_file = args.model
    dx_file = pqr_file + '.dx'  # APBS output file
    vmd_exec = args.vmd_path
    output_dir = args.output
    input_dir = args.input
    sc_score = args.sc
    prefix=pqr_file.split('/')[-1].split('.')[0]
    pdb_file=os.path.basename(pqr_file.replace('pqr','pdb'))
    # Create output directory if needed
    df=pd.read_csv(input_dir,sep=',')
    mask=df['filename']==pdb_file
    
    pdb_refine=pdb_file.replace('.pdb','_refine.pdb')
    in_file=pqr_file.replace('pqr','in')
    command=f'pdbfixer {pdb_file} --output {pdb_refine} &wait & pdb2pqr --noopt --ff AMBER --apbs-input {pdb_refine} {pqr_file} & wait & apbs {in_file}'
    # Initialize output file

    
    # Process each PDB entry
 # Check filename match
    res=df.loc[mask,'interface_residues'].astype(str).iloc[0]
    # interface_res = ast.literal_eval(res)
    interface_res =res.strip('[] \n').split(' ')
    # Map potentials to residues using VMD
    residue_potentials_dict = map_atoms_vmd(
        pqr_file, 
        interface_res, 
        dx_file, 
        vmd_exec
    )
    
    # Calculate residue scores
    score_dict = calculate_residue_scores(residue_potentials_dict)
    score = np.mean(list(score_dict.values())) if score_dict else 0.0
    
    # Update output line
    df.loc[mask,'sc'] = sc_score
    df.loc[mask,'ec'] = score
    print(df.loc[mask])
    
    