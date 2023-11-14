import numpy as np
import argparse
from scipy.spatial import distance
pdborder = ['template','reactant','TS','product']
A = ['F','O','N','S','P','Cl']
def gen_pdb(pdb):
    #print(pdb)
    #pdb = [pdb1,pdb2,pdb3,pdb4]
    pdb = pdb.replace(',',' ')
    pdb = pdb.split()
    return pdb

def gen_shift_dict(pdb):
    shift_dict = {}
    for n in range(len(pdb)):
        with open(pdb[n], 'r') as fp:
            data = fp.readlines()
            for line in data:
                chain = line[21:22].replace(' ','')
                res = int(line[23:26].strip())
                #code = line[17:21].replace(' ','')
                atom = line[11:17].replace(' ','')
                x = float((line[29:38]).replace(' ',''))
                y = float((line[38:46]).replace(' ',''))
                z = float((line[46:54]).replace(' ',''))
                value = [x,y,z]
                if (chain,res,atom) not in shift_dict.keys():
                    shift_dict[chain,res,atom] = [[],[],[],[]]
                    shift_dict[chain,res,atom][n] = value
                else:
                    shift_dict[chain,res,atom][n] = value 
    #print(shift_dict)                   
    return shift_dict

def gen_move(shift_dict):
    shift = {}
    for chain,res,atom in shift_dict.keys():
        shift[chain,res,atom] = []
        v = shift_dict[chain,res,atom] 
        ref = np.array([[v[0][0],v[0][1],v[0][2]]])
        r = np.array([[v[1][0],v[1][1],v[1][2]]])
        ts = np.array([[v[2][0],v[2][1],v[2][2]]])
        p = np.array([[v[3][0],v[3][1],v[3][2]]])
        sr = distance.cdist(ref, r,  'euclidean')
        sts = distance.cdist(ref, ts,  'euclidean') 
        sp = distance.cdist(ref, p,  'euclidean')
        sr = sr[0,0]
        sts = sts[0,0]
        sp = sp[0,0]
        shift[chain,res,atom] = [sr,sts,sp,sts-sr,sp-sts,sp-sr]
    return shift
                  
# present hb in template,reactant,TS and product, and diff based on the template
def write_shift_probe(shift):                                   
    with open('shift.dat', 'w') as f:
        for chain,res,atom in shift.keys():
            v = shift[chain,res,atom] 
            f.write('%-4s'%chain + '%-4s'%res + '%-6s'%atom + '%8.1f'%v[0]+ '%8.1f'%v[1]+ '%8.1f'%v[2]+ '%8.1f'%v[3]+ '%8.1f'%v[4]+ '%8.1f'%v[5])
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
    shift_dict = gen_shift_dict(pdb) 
    shift = gen_move(shift_dict) 
    write_shift_probe(shift)


# python3 shift_1.py -pdb template.pdb,reactant_all_compare_43.pdb,TS_all_compare_16.pdb,product_all_compare_69.pdb