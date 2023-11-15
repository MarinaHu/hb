import numpy as np
import argparse
from scipy.spatial import distance
pdborder = ['template','reactant','TS','product']
A = ['F','O','N','S','P','Cl']
top = ['H06','H01','C9','C8','O7','C11','O4','O3']
def gen_pdb(pdb):
    pdb = pdb.replace(',',' ')
    pdb = pdb.split()
    return pdb

def gen_xyz(pdb):
    top_dict = {}
    bottom_dict = {}
    for n in range(len(pdb)):
        with open(pdb[n], 'r') as fp:
            data = fp.readlines()
            for line in data:
                chain = line[21:22].replace(' ','')
                res = int(line[23:26].strip())
                atom = line[11:17].replace(' ','')
                x = float((line[29:38]).replace(' ',''))
                y = float((line[38:46]).replace(' ',''))
                z = float((line[46:54]).replace(' ',''))
                value = [x,y,z]
                if chain == 'A' and res == 203:
                    if atom in top:
                        if atom not in top_dict.keys():
                            top_dict[atom] = [[],[],[],[]]
                            top_dict[atom][n] = value
                        else:
                            top_dict[atom][n] = value
                    else:
                        if atom not in bottom_dict.keys():
                            bottom_dict[atom] = [[],[],[],[]]
                            bottom_dict[atom][n] = value
                        else:
                            bottom_dict[atom][n] = value
    return top_dict, bottom_dict

def gen_bond(top_dict, bottom_dict):
    bond_dict = {}
    for kt,vt in top_dict.items():
        for kb,vb in bottom_dict.items():
            bond_dict[kt,kb] = []
            for n in range(4):
                top = np.array([[vt[n][0],vt[n][1],vt[n][2]]])
                bottom = np.array([[vb[n][0],vb[n][1],vb[n][2]]])
                bond = distance.cdist(top, bottom,  'euclidean')
                bond = bond[0,0]
                bond_dict[kt,kb].append(bond)
            l = bond_dict[kt,kb]
            #print(l[1])
            r = l[1] - l[0]
            ts = l[2] - l[0]
            p = l[3] - l[0]
            #print(r)
            bond_dict[kt,kb].append(r)
            bond_dict[kt,kb].append(ts)
            bond_dict[kt,kb].append(p)
            #print(bond_dict)
            bond_dict[kt,kb].append(ts-r)
            bond_dict[kt,kb].append(p-ts)
            bond_dict[kt,kb].append(p-r)      
            #print(bond_dict)   
    return bond_dict
                  
def write_shift_probe(bond_dict):                                   
    with open('bond.dat', 'w') as f:
        for atomt,atomb in bond_dict.keys():
            v = bond_dict[atomt,atomb]
            f.write('%-4s'%atomt + '%3s'%atomb)
            #print(v)
            for i in range(len(v)): 
                f.write('%6.1f'%v[i])
            f.write('\n')
    return
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract all the HB interactions between the substrat and residues')
    parser.add_argument('-pdb', dest='pdb', default = 'None', help='protonated pdbfile')
    #parser.add_argument('-s', dest='seed', default = 'None', help='Chain:Resid1,Resid2;Chain:Resid1,Resid2')
    args = parser.parse_args()


    pdb = args.pdb
    #seed = args.seed
    
    pdb = gen_pdb(pdb)
    top_dict, bottom_dict = gen_xyz(pdb) 
    bond_dict= gen_bond(top_dict, bottom_dict) 
    write_shift_probe(bond_dict)


# python3 bond_1.py -pdb template.pdb,reactant_all_compare_43.pdb,TS_all_compare_16.pdb,product_all_compare_69.pdb