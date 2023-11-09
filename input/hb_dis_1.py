import numpy as np
import argparse
from scipy.spatial import distance

A = ['F','O','N','S','P','Cl']

def gen_hb_dis(pdb):
    with open(pdb, 'r') as fp:
        data = fp.readlines()
    cd_seed = {}
    cd_r = {}
    dis_dict = {}
    for line in data:
        chain = line[21:22].replace(' ','')
        res = int(line[23:26].strip())
        code = line[17:21].replace(' ','')
        atom = line[11:17].replace(' ','')
        x = float((line[29:38]).replace(' ',''))
        y = float((line[38:46]).replace(' ',''))
        z = float((line[46:54]).replace(' ',''))
        if chain+str(res) == 'A203':
            cd_seed['A',203,atom] = [x,y,z]
        elif chain+str(res) != 'A203':
           cd_r[chain,res,atom] = [x,y,z] 

    for chain,res,atom in cd_seed.keys():
        if 'H' == atom[0]:
            seed_x,seed_y,seed_z = cd_seed['A',203,atom][0],cd_seed['A',203,atom][1],cd_seed['A',203,atom][2]
            chain_seed,res_seed,atom_seed = chain,res,atom
            for chain,res,atom in cd_r.keys(): 
                for i in A:
                    if i == atom[0]:
                        rx,ry,rz = cd_r[chain,res,atom][0],cd_r[chain,res,atom][1],cd_r[chain,res,atom][2]
                        chain_r,res_r,atom_r = chain,res,atom 
                        X = np.array([[seed_x,seed_y,seed_z ]])
                        Y = np.array([[rx,ry,rz]])
                        dis = distance.cdist(X, Y,  'euclidean')
                        dis = dis[0,0]
                        dis = format(dis,".2f")
                        dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r] = dis
        else:
            for i in A:
                if i == atom[0]:
                    seed_x,seed_y,seed_z = cd_seed['A',203,atom][0],cd_seed['A',203,atom][1],cd_seed['A',203,atom][2]
                    chain_seed,res_seed,atom_seed = chain,res,atom
                    for chain,res,atom in cd_r.keys(): 
                       if 'H' == atom[0]:
                            rx,ry,rz = cd_r[chain,res,atom][0],cd_r[chain,res,atom][1],cd_r[chain,res,atom][2]
                            chain_r,res_r,atom_r = chain,res,atom 
                            X = np.array([[seed_x,seed_y,seed_z ]])
                            Y = np.array([[rx,ry,rz]])
                            dis = distance.cdist(X, Y,  'euclidean')
                            dis = dis[0,0]
                            dis = format(dis,".2f")
                            dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r] = dis
    return dis_dict
                            

def write_hb_probe(dis_dict):                                   
    with open('hb_dist.dat', 'w') as f:
        for chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r in dis_dict.keys():
            dis = dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r]
            if float(dis) < 3.15:            
                f.write('%-4s'%chain_seed + '%-4s'%res_seed + '%-4s'%atom_seed+ '%-4s'%chain_r + '%-4s'%res_r + '%-10s'%atom_r + dis)
                f.write('\n')
    return
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract all the HB interactions between the substrat and residues')
    parser.add_argument('-pdb', dest='pdb', default = 'None', help='protonated pdbfile')
    #parser.add_argument('-s', dest='seed', default = 'None', help='Chain:Resid1,Resid2;Chain:Resid1,Resid2')
    args = parser.parse_args()


    pdb = args.pdb
    #seed = args.seed

    #seed = gen_seed(seed) 
    dis_dict = gen_hb_dis(pdb)
    write_hb_probe(dis_dict)


# python3 hb_2.py -probe 2cht_h_modify.probe -pdb 2cht_h_modify.pdb -s A:203