from vina import Vina
import os
import math 
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from scipy.stats import epps_singleton_2samp
from scipy.stats import anderson_ksamp
import matplotlib.pyplot as plt
import random
from rdkit import Chem
from rdkit.Chem import Lipinski



def bepred(path='',rcpqt=''):
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    if os.path.exists(path+'Progress.st'):
        with open(path+'Progress.st')as fin:
            for line in fin:
                if 'PDBPath' in line:
                    line=line.rstrip().split(':')
                    line=line[-1].split('/')
                    global cwdir
                    global pdb_nm
                    cwdir='/'.join(line[0:-1])+'/'
                    pdb_nm=line[-1].split('.')[0]
                elif 'LigPath' in line:
                    global LigPath
                    line=line.rstrip().split(':')
                    LigPath=line[-1]
                elif 'cavity_grid' in line:
                    line=line.rstrip().split(':')
                    crds=line[-1]
                    global cords
                    cords=[]
                    if len(crds.split('\t'))>1:
                        for i in range(len(crds.split('\t'))):
                            if crds.split('\t')[i]!='':
                                cords.append(float(crds.split('\t')[i]))
                    else:
                        for i in range(len(crds.split(' '))):
                            if crds.split(' ')[i]!='':
                                cords.append(float(crds.split(' ')[i]))

    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try :LigPath
    except: raise NameError('Error: LigBuilder path not found: Please run Synthesizer module first')
    try: cwdir
    except: raise NameError('Error: Error: Could not find deafult output directory: Please run Synthesizer module first')
    try: cords
    except: raise NameError('Error: Necessary cavity information not found: Please run Synthesizer module first')
    os.chdir(cwdir)
    if os.path.exists('Synth.csv'):   
        synth=pd.read_csv('Synth.csv')
    else:
        raise NameError('Error: Synthesized compund data not found: Please run Synthesizer module first')
    if len(rcpqt)==0:
        if os.path.exists(cwdir+pdb_nm+'.pdb'):
            os.system('obabel -ipdb '+cwdir+pdb_nm+'.pdb -xr -O'+cwdir+pdb_nm+'.pdbqt')
            rcpqt=pdb_nm+'.pdbqt'
    else:
        if os.path.exists(rcpqt):
            pass
        else:
            raise NameError('Error:',rcpqt,'not found in the default output directory')
    
    os.chdir(cwdir)
    ot=[]
    with open(pdb_nm+'.pdbqt')as fin:
        for line in fin:
            if line.startswith('ENDROOT'):
                ot.append('TER\n')
            elif line.startswith('TORSDOF'):
                pass
            elif line.startswith('ROOT'):
                pass
            else:
                ot.append(line)
    with open(pdb_nm+'.pdbqt','w')as fout:
        for line in ot:
            fout.write(line)
    #os.system('mkdir pdbqts')
    #os.system('mv '+cwdir+'*.pdbqts pdbqts/')
    xc=(cords[0]+cords[3])/2
    yc=(cords[1]+cords[4])/2
    zc=(cords[2]+cords[5])/2
    bs=int(math.sqrt(pow((cords[5]-cords[2]),2)+pow((cords[4]-cords[1]),2)+pow((cords[3]-cords[0]),2))/math.sqrt(3))+2
    all_be=[]
    fnd_lig=[]
    for lig in synth['ID'].tolist():
        if os.path.exists(lig+'.pdbqt'):
            v = Vina(sf_name='vina')
            v.set_receptor(cwdir+rcpqt)
            v.set_ligand_from_file(lig+'.pdbqt')
            v.compute_vina_maps(center=[xc,yc,zc], box_size=[bs,bs,bs])
            v.dock(exhaustiveness=9, n_poses=20)
            all_be.append(min(v.score()))
            fnd_lig.append(lig)
            print(lig,min(v.score()))
        else:
            #print('warning: some synthesized molecules seems to be missing')
            print('Missing:'+lig)
    if len(all_be)==len(synth['ID'].tolist()):
        synth["BE"] = all_be
        synth.to_csv(cwdir+'Synth_BE.csv')
        f = open(cwdir+"Progress.st", "a")
        f.write('7. BE output:'+cwdir+'Synth_BE.csv\n')
        f.close()
        return synth,cwdir
    else:
        df=pd.DataFrame()
        df['BE']=all_be
        df['ID']=fnd_lig
        result = pd.merge(synth, df, on="ID")
        result.to_csv(cwdir+'Synth_BE.csv')
        f = open(cwdir+"Progress.st", "a")
        f.write('7. BE output:'+cwdir+'Synth_BE.csv\n')
        f.close()
        return result,cwdir

def ecdf(data):
    n = len(data)
    x = np.sort(data)
    y = np.arange(1, n+1) / n
    return x, y

def synthetic_interaction(path='',pdbqt='',method='KS-test',p_val=0.05,plot='Density'):
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    res,cwdir = bepred(path,pdbqt)
    cutoffs=[]
    for i in range(int(min(res['BE'].tolist())),int(max(res['BE'].tolist()))):
        ls1=[]
        ls2=[]
        for j in res['BE'].tolist():
            if j<=i:
                ls1.append(j)
            else:
                ls2.append(j)
        if len(ls1)<=5 or len(ls2)<=5:
            continue
        if method=='KS-test':
                testr = ks_2samp(ls1, ls2)
        if method=='ES-test':
                testr = epps_singleton_2samp(ls1, ls2)
        if method=='AD-test':
                testr = anderson_ksamp([np.array(ls1), np.array(ls2)])
        if len(testr) > 2:
                if testr[2] < p_val:
                    if len(ls1)/len(res['BE'].tolist()) <=0.7 and len(ls1)/len(res['BE'].tolist()) >=0.5:
                        cutoffs.append(i)
                        if plot=='ECDF':
                            x1, y1 = ecdf(ls1)
                            x2, y2 = ecdf(ls2)
                            fig = plt.figure(figsize=(7, 7),linewidth = 4)
                            plt.plot(x1, y1, marker='.')
                            plt.plot(x2, y2, marker='.')
                            plt.title('Cutoff:'+str(i)+' '+method+' Statistics='+str(testr[0])+' p-value='+str(testr[2]))
                            plt.xlabel("Binding Affinity")
                            plt.ylabel('ECDF')
                            plt.legend(['HAB', 'LAB'])
                            plt.savefig(cwdir+'ECDF_at_cutoff('+str(i)+')_with_'+method+'.pdf')
                            plt.show()
                        else:
                            ls1_lb=['HAB']*len(ls1)
                            ls2_lb=['LAB']*len(ls2)
                            ls1.extend(ls2)
                            ls1_lb.extend(ls2_lb)
                            df=pd.DataFrame()
                            df['BE']=ls1
                            df['status']=ls1_lb
                            data_wide = df.pivot(columns='status',values='BE')
                            data_wide.plot.density(figsize = (7, 7),linewidth = 4)
                            plt.title('Cutoff:'+str(i)+' '+method+' Statistics='+str(testr[0])+' p-value='+str(testr[2]))
                            plt.xlabel("Binding Affinity")
                            plt.savefig(cwdir+'Density_at_cutoff('+str(i)+')_with_'+method+'.pdf')
                            plt.show()
        elif testr[1] < p_val:
            if len(ls1)/len(res['BE'].tolist()) <=0.7 and len(ls1)/len(res['BE'].tolist()) >=0.5:
                cutoffs.append(i)
                if plot=='ECDF':
                    x1, y1 = ecdf(ls1)
                    x2, y2 = ecdf(ls2)
                    fig = plt.figure(figsize=(7, 7),linewidth = 4)
                    plt.plot(x1, y1, marker='.')
                    plt.plot(x2, y2, marker='.')
                    plt.title('Cutoff:'+str(i)+' '+method+' Statistics='+str(testr[0])+' p-value='+str(testr[1]))
                    plt.xlabel("Binding Affinity")
                    plt.ylabel('ECDF')
                    plt.legend(['HAB', 'LAB'])
                    plt.savefig(cwdir+'ECDF_at_cutoff('+str(i)+')_with_'+method+'.pdf')
                    plt.show()
                else:
                    ls1_lb=['HAB']*len(ls1)
                    ls2_lb=['LAB']*len(ls2)
                    ls1.extend(ls2)
                    ls1_lb.extend(ls2_lb)
                    df=pd.DataFrame()
                    df['BE']=ls1
                    df['status']=ls1_lb
                    data_wide = df.pivot(columns='status',values='BE')
                    data_wide.plot.density(figsize = (7, 7),linewidth = 4)
                    plt.title('Cutoff:'+str(i)+' '+method+' Statistics='+str(testr[0])+' p-value='+str(testr[1]))
                    plt.xlabel("Binding Affinity")
                    plt.savefig(cwdir+'Density_at_cutoff('+str(i)+')_with_'+method+'.pdf')
                    plt.show()
                
    if len(cutoffs)==0:
        ls1=[]
        ls2=[]
        for j in res['BE'].tolist():
            if j<=-7:
                ls1.append(j)
            else:
                ls2.append(j)
        if len(ls1)<250:
            print('Minimal requirement of HAB does not statisfy\nTarget cavity is not optimal for drug target!!')
        else:
            print('The synthesized molecules are not calssifiable on the basis of molecular interaction\nPlease proceed  with HAB Decoys.')
    else:
        print('Suggested cutoffs with ',method,' method:\n',cutoffs)
    return res

def create_decoy_smiles(smiles_list):
    decoy_smiles_list = []
    elements =['[H]','[He]','[Li]','[Be]','B','C','N','O','F','[Ne]','[Na+]','[Mg+2]','[Al+3]','[Si]','P','S','Cl','[Ar]','[K+]','[Ca+2]','[Sc+3]','[Ti+2]','[V+3]','[Cr+3]','[Mn+2]','[Fe+2]','[Co+2]','[Ni+2]','[Cu+2]','[Zn+2]','[Ga+3]','[Ge]','[As]','[Se]','Br','[Kr]','[Rb+]','[Sr+2]','[Y+3]','[Zr+4]','[Nb+3]','[Mo]','[Tc]','[Ru+3]','[Rh+3]','[Pd+2]','[Ag+]','[Cd+2]','[In+3]','[Sn]','[Sb]','[Te]','I','[Xe]','[Cs+]','[Ba+2]','[La+3]','[Ce+3]','[Pr+3]','[Nd+3]','[Pm+3]','[Sm+3]','[Eu+2]','[Gd+3]','[Tb+3]','[Dy+3]','[Ho+3]','[Er+3]','[Tm+3]','[Yb+3]','[Lu+3]','[Hf+4]','[Ta+5]','[W+6]','[Re]','[Os]','[Ir+3]','[Pt+2]','[Au+]','[Hg+2]','[Tl+]','[Pb]','[Bi+3]','[Po]','[At]','[Rn]','[Fr+]','[Ra+2]','[Ac+3]','[Th+4]','[Pa+5]','[U+6]','[Np+6]','[Pu+6]','[Am+3]','[Cm+3]','[Bk+3]','[Cf+3]','[Es+3]','[Fm+3]','[Md+3]','[No+2]','[Lr+3]','[Rf]','[Db]','[Sg]','[Bh]','[Hs]','[Mt]','[Ds]','[Rg]','[Cn]','[Nh]','[Fl]','[Mc]','[Lv]','[Ts]','[Og]']
    query_smiles_list=[]
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        atoms = mol.GetAtoms()
        atom = random.choice(atoms)
        sub_ele=[]
        for ele in elements :
            elmol=Chem.MolFromSmiles(ele)
            eltoms = elmol.GetAtoms()
            eltom = random.choice(eltoms)
            if atom.GetTotalValence() == eltom.GetTotalValence():
                if ele != atom.GetSymbol():
                    sub_ele.append(ele)
        if len(sub_ele)==0:
            continue
        query_smiles_list.append(smiles)
        element = random.choice(sub_ele)
        atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(element))
        decoy_smiles = Chem.MolToSmiles(mol)
        decoy_smiles_list.append(decoy_smiles)
    return decoy_smiles_list

def synthetic_decoys(path='',cf=-7,decoy_csv=''):
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    if os.path.exists(path+'Progress.st'):
        with open(path+'Progress.st')as fin:
            for line in fin:
                if 'PDBPath' in line:
                    line=line.rstrip().split(':')
                    line=line[-1].split('/')
                    global cwdir
                    global pdb_nm
                    cwdir='/'.join(line[0:-1])+'/'
                    pdb_nm=line[-1].split('.')[0]
                elif 'BE output' in line:
                    line=line.rstrip().split(':')[-1]
                    be_dat=pd.read_csv(line)
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try: cwdir
    except: raise NameError('Error: Error: Could not find deafult output directory: Please run Synthesizer module first')
    if len(decoy_csv)>0:
        be_dat=pd.read_csv(decoy_csv)
        lb=[]
        st=[]
        anot=list(be_dat['Annotation'].tolist())
        for j in range(len(anot)):
            if anot[j]=='HAB':
                lb.append('HAB')
                st.append(1)
            else:
                lb.append('Decoy')
                st.append(0)
                
        if lb.count('HAB')<250 or lb.count('Decoy')<150:
            raise NameError('Minimal requirement of HAB/Decoys does not statisfy\nDataset provided either too small wrongly annotatted: Please recheck the cloumn names and annotaion type')
        else:
            df=pd.DataFrame()
            df['SMILES']=list(be_dat['SMILES'].tolist())
            df['status']=st
            df['Label']=lb
            df.to_csv(cwdir+'Labeled_cmps.csv')
            f = open(cwdir+"Progress.st", "a")
            f.write('8. CMP labeled:'+cwdir+'Labeled_cmps.csv\n')
            f.close()
            return df
    else:
        try: be_dat
        except: raise NameError('Error: Synthetic comlund interaction data not found: Please run synthetic_interaction function first')
    ls1=[]
    ls2=[]
    sms=list(be_dat['SMILES'].tolist())
    bes=list(be_dat['BE'].tolist())
    for j in range(len(bes)):
        if bes[j]<=cf:
            ls1.append(sms[j])
        else:
            ls2.append(sms[j])
    if len(ls1)<250:
        print('Minimal requirement of HAB does not statisfy\nTarget cavity is not optimal for drug target!!')
    else:
        dcys=create_decoy_smiles(ls1)
        ls1_lb=['HAB']*len(ls1)
        ls2_lb=['LAB']*len(dcys)
        ls1_st=[1]*len(ls1)
        ls2_st=[0]*len(dcys)
        ls1.extend(dcys)
        ls1_lb.extend(ls2_lb)
        ls1_st.extend(ls2_st)
        df=pd.DataFrame()
        df['SMILES']=ls1
        df['status']=ls1_st
        df['Label']=ls1_lb
        df.to_csv(cwdir+'Labeled_cmps.csv')
        f = open(cwdir+"Progress.st", "a")
        f.write('8. CMP labeled:'+cwdir+'Labeled_cmps.csv\n')
        f.close()
        return df


def synthetic_classify(path='',cf=-7):
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    if os.path.exists(path+'Progress.st'):
        with open(path+'Progress.st')as fin:
            for line in fin:
                if 'PDBPath' in line:
                    line=line.rstrip().split(':')
                    line=line[-1].split('/')
                    global cwdir
                    global pdb_nm
                    cwdir='/'.join(line[0:-1])+'/'
                    pdb_nm=line[-1].split('.')[0]
                elif 'BE output' in line:
                    line=line.rstrip().split(':')[-1]
                    dat=pd.read_csv(line)
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try: cwdir
    except: raise NameError('Error: Error: Could not find deafult output directory: Please run Synthesizer module first')
    try: dat
    except: raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    ls1=[]
    ls1_sm=[]
    ls2=[]
    ls2_sm=[]
    BEs=dat['BE'].tolist()
    SMs=dat['SMILES'].tolist()
    for j in range(len(BEs)):
        if BEs[j]<=cf:
            ls1.append(BEs[j])
            ls1_sm.append(SMs[j])
        else:
            ls2.append(BEs[j])
            ls2_sm.append(SMs[j])
    ls1_lb=['HAB']*len(ls1)
    ls2_lb=['LAB']*len(ls2)
    ls1_st=[1]*len(ls1)
    ls2_st=[0]*len(ls2)
    ls1.extend(ls2)
    ls1_lb.extend(ls2_lb)
    ls1_st.extend(ls2_st)
    ls1_sm.extend(ls2_sm)
    df=pd.DataFrame()
    df['SMILES']=ls1_sm
    df['BE']=ls1
    df['status']=ls1_st
    df['Label']=ls1_lb
    df.to_csv(cwdir+'Labeled_cmps.csv')
    f = open(cwdir+"Progress.st", "a")
    f.write('8. CMP labeled:'+cwdir+'Labeled_cmps.csv\n')
    f.close()
    return df
