
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



set sel [atomselect top 'not name H*']
set all_data {}

foreach index [$sel list] {
    set atom [atomselect top "index $index"]
    set res_id [lindex [$atom get resid] 0]
    set chain
     [lindex [$atom get chain] 0]
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
