import numpy as np
import argparse
from scipy.spatial import distance

A = ['F','O','N','S','P','Cl']

def gen_seed(seed_raw):
    seed_raw = seed_raw.replace(';',' ')
    seed_raw = seed_raw.replace(':',' ')
    seed_raw = seed_raw.split()
    seed = []
    for i in range(len(seed_raw)//2):
        seed_id = []
        seed_id = seed_raw[2*i+1].replace(',',' ')
        seed_id = seed_id.split()
        for j in range(len(seed_id)):
            seed.append(seed_raw[2*i] + '-' + seed_id[j])
    return seed

def gen_xyz(a,pdb):
    X = np.empty(shape=(3,1))
    with open(pdb, 'r') as fp:  
        data = fp.readlines() 
        for line in data:
            atom = line[11:17].replace(' ','')
            code = line[17:21].replace(' ','')
            chain = line[21:22].replace(' ','')
            res = line[23:26].replace(' ','')
            x = line[31:38].replace(' ','')
            y = line[38:46].replace(' ','')
            z = line[46:54].replace(' ','')
            pdb = '%-4s'%chain + '%-4s'%res + '%-5s'%code + '%-5s'%atom
            if pdb == a:
                x = float(x)
                y = float(y)
                z = float(z)
                X = np.array([[x,y,z]])  
    return X

def calc_dis(a,b,pdb):
    X = gen_xyz(a,pdb)
    Y = gen_xyz(b,pdb)
    dis = distance.cdist(X, Y,  'euclidean')
    dis = dis[0,0]
    dis = format(dis,".2f")    
    return dis

def H_loc(atom1,atom2,a,b,hb,intype,pdb):
    if 'H' == atom1[0]:   
        for i in A:
            if i == atom2[0]:
                if (a + '%-4s'%'<->' + b + '(' + intype + ')') not in hb.keys():
                    dis = calc_dis(a,b,pdb) 
                    hb[a + '%-4s'%'<->' + b+ '(' + intype + ')'] = [dis, 1]
                else:
                    hb[a + '%-4s'%'<->' + b+ '(' + intype + ')'][1] += 1
    return hb

def situation(hb,atom1,atom2,intype,a,b,pdb,a_res,b_res,seed):
    if (a_res in seed) and (b_res not in seed):  
        H_loc(atom1,atom2,a,b,hb,intype,pdb)
        H_loc(atom2,atom1,a,b,hb,intype,pdb)
    return hb

def gen_hb(probe,pdb,seed):
    with open(probe, 'r') as fp:
        data = fp.readlines()
        hb = {}
        for line in data:
            intype = line[6:8].replace(' ','')
            chain1 = line[9:11].replace(' ','')
            res1 = int(line[11:16].strip())
            code1 = line[16:20].replace(' ','')
            atom1 = line[20:25].replace(' ','')
            chain2 = line[25:28].replace(' ','')
            chain2 = chain2.replace(':','')
            res2 = line[28:33].replace(' ','')
            res2 = res2.replace(chain2,'')
            res2 = int(res2)
            code2 = line[32:36].replace(' ','')
            atom2 = line[37:41].replace(' ','')
            a = '%-4s'%chain1 + '%-4s'%str(res1) + '%-5s'%code1 + '%-5s'%atom1
            b = '%-4s'%chain2 + '%-4s'%str(res2) + '%-5s'%code2 + '%-5s'%atom2
            a_res = chain1+'-'+ str(res1)
            b_res = chain2+'-'+ str(res2)            
            situation(hb,atom1,atom2,intype,a,b,pdb,a_res,b_res,seed)
            situation(hb,atom1,atom2,intype,b,a,pdb,b_res,a_res,seed)               
    return hb

def write_hb_probe(hb):                                   
    with open('hb.probe', 'w') as f:
        for k,v in hb.items():
            f.write('%-40s'%k + '%-4s'%':' + '%-2s'%v[0]+ '%-2s'%',' + '%-4s'%v[1])
            f.write('\n')
    return
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract all the HB interactions between the substrat and residues')
    parser.add_argument('-probe', dest='probe', default = 'None', help='protonated probefile')
    parser.add_argument('-pdb', dest='pdb', default = 'None', help='protonated pdbfile')
    parser.add_argument('-s', dest='seed', default = 'None', help='Chain:Resid1,Resid2;Chain:Resid1,Resid2')
    args = parser.parse_args()

    probe = args.probe
    pdb = args.pdb
    seed = args.seed

    seed = gen_seed(seed) 
    hb = gen_hb(probe,pdb,seed)
    write_hb_probe(hb)


# python3 hb_2.py -probe 2cht_h_modify.probe -pdb 2cht_h_modify.pdb -s A:203