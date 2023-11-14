# This code extracts all the hb interactions between the substrat and residues
It shows the distance and the counting of the same atom pairs involving hb.

command (in the input folder): python3 hb_2.py -probe 2cht_h_modify.probe -pdb 2cht_h_modify.pdb -s A:203

# This code shows the coordinate shift for each atom in the reaction
command (in the input folder): python3 shift_1.py -pdb template.pdb,reactant_all_compare_43.pdb,TS_all_compare_16.pdb,product_all_compare_69.pdb