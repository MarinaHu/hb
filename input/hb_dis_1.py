import numpy as np
import argparse
from scipy.spatial import distance
pdborder = ['template','reactant','TS','product']
A = ['F','O','N','S','P','Cl']
def gen_pdb(pdb):
    pdb = pdb.replace(',',' ')
    pdb = pdb.split()
    return pdb

def gen_hb_dis(pdb):
    dis_dict = {}
    for n in range(len(pdb)):
        with open(pdb[n], 'r') as fp:
            data = fp.readlines()
            cd_seed = {}
            cd_r = {}
            for line in data:
                chain = line[21:22].replace(' ','')
                res = int(line[23:26].strip())
                #code = line[17:21].replace(' ','')
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
                    #print(cd_seed['A',203,atom])
                    chain_seed,res_seed,atom_seed = chain,res,atom
                    for chain,res,atom in cd_r.keys(): 
                        for i in A:
                            if i == atom[0]:
                                rx,ry,rz = cd_r[chain,res,atom][0],cd_r[chain,res,atom][1],cd_r[chain,res,atom][2]
                                chain_r,res_r,atom_r = chain,res,atom
                                if (chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r) not in dis_dict.keys():
                                    dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r] = [0.0,0.0,0.0,0.0]
                                
                                X = np.array([[seed_x,seed_y,seed_z]])
                                Y = np.array([[rx,ry,rz]])
                                dis = distance.cdist(X, Y,  'euclidean')
                                dis = dis[0,0]
                                #dis = format(dis,".3f")
                                dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r][n] = dis
                                #print(dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r],n)
                else:
                    for i in A:
                        if i == atom[0]:
                            seed_x,seed_y,seed_z = cd_seed['A',203,atom][0],cd_seed['A',203,atom][1],cd_seed['A',203,atom][2]
                            chain_seed,res_seed,atom_seed = chain,res,atom
                            for chain,res,atom in cd_r.keys(): 
                                if 'H' == atom[0]:
                                    rx,ry,rz = cd_r[chain,res,atom][0],cd_r[chain,res,atom][1],cd_r[chain,res,atom][2]
                                    chain_r,res_r,atom_r = chain,res,atom
                                    if (chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r) not in dis_dict.keys():
                                        dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r] = [0.0,0.0,0.0,0.0] 
                                    X = np.array([[seed_x,seed_y,seed_z]])
                                    Y = np.array([[rx,ry,rz]])
                                    dis = distance.cdist(X, Y,  'euclidean')
                                    dis = dis[0,0]
                                    #dis = format(dis,".3f")
                                    dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r][n] = dis
    return dis_dict
                                
# present hb in template,reactant,TS and product, and diff based on the template
def write_hb_probe(dis_dict):                                   
    with open('hb_dist.dat', 'w') as f:
        #print(dis_dict)
        for chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r in dis_dict.keys():
            dis = dis_dict[chain_seed,res_seed,atom_seed,chain_r,res_r,atom_r]
            r = dis[1]-dis[0]
            ts = dis[2]-dis[0]
            p = dis[3]-dis[0]
            dts = ts - r
            dpts = p-ts
            dpr = p-r

            """""
            r = format(r,".1f")
            ts = format(ts,".1f")
            p = format(p,".1f")
            dts = format(dts,".1f")
            dpts = format(dpts,".1f")
            dpr = format(dpr,".1f")
"""
            if (dis[0] < 5) and (dis[1] < 5) and (dis[2] < 5) and (dis[3] < 5):
                #dis[0] = format(dis[0],".3f")
                #dis[1] = format(dis[1],".3f")
                #dis[2] = format(dis[2],".3f")
                #dis[3] = format(dis[3],".3f")
                f.write('%-4s'%chain_seed + '%-4s'%res_seed + '%-4s'%atom_seed+ '%-4s'%chain_r + '%-4s'%res_r + '%-10s'%atom_r + '%8.3f'%dis[0]+ '%8.3f'%dis[1]+ '%8.3f'%dis[2]+ '%8.3f'%dis[3]+ '%8.1f'%r + '%8.1f'%ts+ '%8.1f'%p+ '%8.1f'%dts + '%8.1f'%dpts+ '%8.1f'%dpr)
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
    dis_dict = gen_hb_dis(pdb)
    write_hb_probe(dis_dict)


# python3 hb_dis_1.py -pdb template.pdb,reactant_all_compare_43.pdb,TS_all_compare_16.pdb,product_all_compare_69.pdb