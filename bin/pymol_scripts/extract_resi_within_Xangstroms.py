# Script to extract residues within a certain number of angstroms from the active site
# import sys
# sys.stdout = open('8angstroms.txt', 'w')

# Fetch model
fetch 4KU5_A
hide everything, 4KU5_A
show cartoon, 4KU5_A
set seq_view, 1
color marine, 4KU5_A

# Select catalytic residues
select S143, resi 143

# Select all residues within 8 angstroms
select all_res, *
select 8angstroms, all_res within 8 of S143
remove resn hoh
iterate 8angstroms (name ca), print(resn, resi)

