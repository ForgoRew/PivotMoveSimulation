#logfile vmd.log

proc addbond { p q } {
    global b
    set a [ lreplace $b $p $p [ lreplace [ lindex $b $p ] 0 -1 $q ] ]
    lassign [ list $a ] b
    set a [ lreplace $b $q $q [ lreplace [ lindex $b $q ] 0 -1 $p ] ]
    lassign [ list $a ] b
}
proc remove_long_bonds { max_length } {
    for { set i 0 } { $i < [ molinfo top get numatoms ] } { incr i } {
        set bead [ atomselect top "index $i" ]
        set bonds [ lindex [$bead getbonds] 0 ]
        if { [ llength bonds ] > 0 } {
            set bonds_new {}
            set xyz [ lindex [$bead get {x y z}] 0 ]
            foreach j $bonds {
                set bead_to [ atomselect top "index $j" ]
                set xyz_to [ lindex [$bead_to get {x y z}] 0 ]
                if { [ vecdist $xyz $xyz_to ] < $max_length } {
                    lappend bonds_new $j
                }
            }
        $bead setbonds [ list $bonds_new ]
        }
    }
}
set s [atomselect top "all"]
set b {}
for {set i 0} {$i < [$s num]} {incr i} {
    lappend b {}
}

#for {set i 0} {$i < [$s num]} {incr i} {
#    addbond $i $i+1
#    }

addbond 0 1
addbond 1 2
addbond 2 3
addbond 3 4
addbond 4 5
addbond 5 6
addbond 6 7
addbond 7 8
addbond 8 9
addbond 9 10
addbond 10 11
addbond 11 12
addbond 12 13
addbond 13 14
addbond 14 15
addbond 15 16
addbond 16 17
addbond 17 18
addbond 18 19
$s setbonds $b
remove_long_bonds 5.0
pbc set {1000 1000 1000 90.0 90.0 90.0} -all
# pbc box_draw
# pbc wrap -all
axes location off
color Display Background white
#color Display Background black
pbc wrap -center com

#COLOR definitions
#------------------------------------------
#red
#ColorID 1
#silver
#color change rgb 6 0.300000 0.300000 0.40000
#blue
#color change rgb 0 0.100000 0.100000 1.000000
#green
#color change rgb 7 0.200000 0.550000 0.300000
#cyan
#color change rgb 10 0.250000 0.500000 0.750000
#orange
#color change rgb 3 1.000000 0.300000 0.000000

#PARTICLE REPRESENTATIONS
#----------------------------------------
#all polymer segments
mol selection { name 'H:' or 'P:'}
mol material Opaque
mol color ColorID 6
mol representation Licorice 0.5 12 12
mol addrep top

## Hydrofilní (H)
#mol selection { name 'H:' }
#mol color ColorID 0
#mol material Opaque
#mol representation CPK 3.81 1 20 1
#mol addrep top
#
#    # Hydrofóbní (P)
#mol selection { name 'P:' }
#mol color ColorID 1
#mol material Opaque
#mol representation CPK 3.81 1 20 1
#mol addrep top

# FASTA aminoacid code
# A (Alanine)
mol selection { name 'A:' }
mol color ColorID 0
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# C (Cysteine)
mol selection { name 'C:' }
mol color ColorID 1
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# D (Aspartic acid)
mol selection { name 'D:' }
mol color ColorID 2
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# E (Glutamic acid)
mol selection { name 'E:' }
mol color ColorID 3
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# F (Phenylalanine)
mol selection { name 'F:' }
mol color ColorID 4
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# G (Glycine)
mol selection { name 'G:' }
mol color ColorID 5
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# H (Histidine)
mol selection { name 'H:' }
mol color ColorID 6
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# I (Isoleucine)
mol selection { name 'I:' }
mol color ColorID 7
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# K (Lysine)
mol selection { name 'K:' }
mol color ColorID 8
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# L (Leucine)
mol selection { name 'L:' }
mol color ColorID 9
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# M (Methionine)
mol selection { name 'M:' }
mol color ColorID 10
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# N (Asaparagine)
mol selection { name 'N:' }
mol color ColorID 11
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# P (Proline)
mol selection { name 'P:' }
mol color ColorID 12
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# Q (Glutamine)
mol selection { name 'Q:' }
mol color ColorID 13
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# R (Arginine)
mol selection { name 'R:' }
mol color ColorID 14
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# S (Serine)
mol selection { name 'S:' }
mol color ColorID 15
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# T (Threonine)
mol selection { name 'T:' }
mol color ColorID 16
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# V (Valine)
mol selection { name 'V:' }
mol color ColorID 17
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# W (Tryptophane)
mol selection { name 'W:' }
mol color ColorID 18
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
# Y (Tyrosine)
mol selection { name 'Y:' }
mol color ColorID 19
mol material Opaque
mol representation CPK 3.81 1 20 1
mol addrep top
