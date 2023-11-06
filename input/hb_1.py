import numpy as np
import os
import argparse
from scipy.spatial import distance

A = ['F','O','N','S','P','Cl']
def gen_hb(probe,pdb):
    with open(probe, 'r') as fa:
        with open(pdb, 'r') as fb:  
            data = fa.readlines()
            datab = fb.readlines()
  
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
                if (res1 == 203) and (res2 != 203):  
                    if 'H' == atom1[0]:   
                        for i in A:
                            if i == atom2[0]:
                                if (a + '%-4s'%'<->' + b) not in hb.keys():
                                    for lineb in datab:
                                        atom = lineb[11:17].replace(' ','')
                                        code = lineb[17:21].replace(' ','')
                                        chain = lineb[21:22].replace(' ','')
                                        res = lineb[23:26].replace(' ','')
                                        x = lineb[31:38].replace(' ','')
                                        y = lineb[38:46].replace(' ','')
                                        z = lineb[46:54].replace(' ','')
                                        pdb = '%-4s'%chain + '%-4s'%res + '%-5s'%code + '%-5s'%atom
                                        if pdb == a:
                                            chain1 = chain
                                            res1 = res
                                            code1 = code
                                            atom1 = atom
                                            x1 = float(x)
                                            y1 = float(y)
                                            z1 = float(z)
                                        if pdb == b:
                                            chain2 = chain
                                            res2 = res
                                            code2 = code
                                            atom2 = atom
                                            x2 = float(x)
                                            y2 = float(y)
                                            z2 = float(z)
                                    X = np.array([[x1,y1,z1]])
                                    Y = np.array([[x2,y2,z2]])
                                    dis = distance.cdist(X, Y,  'euclidean')
                                    dis = dis[0,0]
                                    dis = format(dis,".2f")
                                    hb[a + '%-4s'%'<->' + b] = [dis, 1, intype]
                                else:
                                    hb[a + '%-4s'%'<->' + b][1] += 1
                    elif 'H' == atom2[0]:   
                        for i in A:
                            if i == atom1[0]:
                                if (a + '%-4s'%'<->' + b) not in hb.keys():
                                    for lineb in datab:
                                        atom = lineb[11:17].replace(' ','')
                                        code = lineb[17:21].replace(' ','')
                                        chain = lineb[21:22].replace(' ','')
                                        res = lineb[23:26].replace(' ','')
                                        x = lineb[31:38].replace(' ','')
                                        y = lineb[38:46].replace(' ','')
                                        z = lineb[46:54].replace(' ','')
                                        pdb = '%-4s'%chain + '%-4s'%res + '%-5s'%code + '%-5s'%atom
                                        if pdb == a:
                                            chain1 = chain
                                            res1 = res
                                            code1 = code
                                            atom1 = atom
                                            x1 = float(x)
                                            y1 = float(y)
                                            z1 = float(z)
                                        if pdb == b:
                                            chain2 = chain
                                            res2 = res
                                            code2 = code
                                            atom2 = atom
                                            x2 = float(x)
                                            y2 = float(y)
                                            z2 = float(z)
                                    X = np.array([[x1,y1,z1]])
                                    Y = np.array([[x2,y2,z2]])
                                    dis = distance.cdist(X, Y,  'euclidean')
                                    dis = dis[0,0]
                                    dis = format(dis,".2f")
                                    hb[a + '%-4s'%'<->' + b] = [dis, 1, intype]
                                else:
                                    hb[a + '%-4s'%'<->' + b][1] += 1
                elif (res2 == 203) and (res1 != 203):
                    if 'H' == atom1[0]:   
                        for i in A:
                            if i == atom2[0]:
                                
                                if (b + '%-4s'%'<->' + a) not in hb.keys():
                                    for lineb in datab:
                                        atom = lineb[11:17].replace(' ','')
                                        code = lineb[17:21].replace(' ','')
                                        chain = lineb[21:22].replace(' ','')
                                        res = lineb[23:26].replace(' ','')
                                        x = lineb[31:38].replace(' ','')
                                        y = lineb[38:46].replace(' ','')
                                        z = lineb[46:54].replace(' ','')
                                        pdb = '%-4s'%chain + '%-4s'%res + '%-5s'%code + '%-5s'%atom
                                        if pdb == a:
                                            chain1 = chain
                                            res1 = res
                                            code1 = code
                                            atom1 = atom
                                            x1 = float(x)
                                            y1 = float(y)
                                            z1 = float(z)
                                        if pdb == b:
                                            chain2 = chain
                                            res2 = res
                                            code2 = code
                                            atom2 = atom
                                            x2 = float(x)
                                            y2 = float(y)
                                            z2 = float(z)
                                    X = np.array([[x1,y1,z1]])
                                    Y = np.array([[x2,y2,z2]])
                                    dis = distance.cdist(X, Y,  'euclidean')
                                    dis = dis[0,0]
                                    dis = format(dis,".2f")
                                    hb[b + '%-4s'%'<->' + a] = [dis, 1, intype]
                                else:
                                    hb[b + '%-4s'%'<->' + a][1] += 1
                    elif 'H' == atom2[0]:   
                        for i in A:
                            if i == atom1[0]:
                                if (b + '%-4s'%'<->' + a) not in hb.keys():
                                    for lineb in datab:
                                        atom = lineb[11:17].replace(' ','')
                                        code = lineb[17:21].replace(' ','')
                                        chain = lineb[21:22].replace(' ','')
                                        res = lineb[23:26].replace(' ','')
                                        x = lineb[31:38].replace(' ','')
                                        y = lineb[38:46].replace(' ','')
                                        z = lineb[46:54].replace(' ','')
                                        pdb = '%-4s'%chain + '%-4s'%res + '%-5s'%code + '%-5s'%atom
                                        if pdb == a:
                                            chain1 = chain
                                            res1 = res
                                            code1 = code
                                            atom1 = atom
                                            x1 = float(x)
                                            y1 = float(y)
                                            z1 = float(z)
                                        if pdb == b:
                                            chain2 = chain
                                            res2 = res
                                            code2 = code
                                            atom2 = atom
                                            x2 = float(x)
                                            y2 = float(y)
                                            z2 = float(z)
                                    X = np.array([[x1,y1,z1]])
                                    Y = np.array([[x2,y2,z2]])
                                    dis = distance.cdist(X, Y,  'euclidean')
                                    dis = dis[0,0]
                                    dis = format(dis,".2f")
                                    hb[b + '%-4s'%'<->' + a] = [dis, 1, intype]
                                else:
                                    hb[b + '%-4s'%'<->' + a][1] += 1
                                    
    with open('hb.probe', 'w') as f:
        for k,v in hb.items():
            f.write('%-40s'%k + '%-4s'%':' + '%-2s'%v[0]+ '%-2s'%',' + '%-4s'%v[1] + '%-4s'%v[2])
            f.write('\n')

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Extract all the HB interactions between the substrat and residues')
    parser.add_argument('-probe', dest='probe', default = 'None', help='protonated probefile')
    parser.add_argument('-pdb', dest='pdb', default = 'None', help='protonated pdbfile')
    args = parser.parse_args()

    probe = args.probe
    pdb = args.pdb

    gen_hb(probe,pdb)

