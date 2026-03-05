source findcontactsV2-3.tcl
puts done source normall
set current_read [file rootname [molinfo top get name]]
calculate_interface_sasa "$current_read" "index 5256 to 46602" "index 1 to 5255"
exit
