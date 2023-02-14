import os
import shutil
import pandas as pd
import py3Dmol
from datetime import datetime

def set_LigBuilder(path_to_ligb,Output_dir):
        global LigB_path
        if os.path.exists(path_to_ligb+'/build/'):
                LigB_path=path_to_ligb+'/build/'
        elif os.path.exists(path_to_ligb+'/LigBuilderV3.0/build/'):
                LigB_path=path_to_ligb+'/LigBuilderV3.0/build/'
        elif os.path.exists(path_to_ligb+'build/'):
                LigB_path=path_to_ligb+'build/'
        elif os.path.exists(path_to_ligb+'LigBuilderV3.0/build/'):
                LigB_path=path_to_ligb+'LigBuilderV3.0/build/'
        else:
                raise NameError('Error: LigBuilderV3 not found: No such file or directory')
        f = open(Output_dir+"Progress.st", "a")
        f.write('1. LigPath:'+LigB_path+'\n')
        f.close()
        
def set_libiomp5(lpath,Output_dir):
        if lpath.endswith('/'):
                lpath=lpath[:-1]
        global lib_path
        if 'libiomp5.so' in lpath:
                lpath=lpath.split('/')
                npath=lpath[:-1]
                if os.path.exists('/'.join(npath)):
                        lib_path='/'.join(npath)
                else:
                        raise NameError('Error: libiomp5.so not found: No such file or directory')
        else:
                if os.path.exists(lpath+'/libiomp5.so'):
                        lib_path=lpath
                else:
                        raise NameError('Error: libiomp5.so not found: No such file or directory')
        '''if os.path.exists(lpath+'/libiomp5.so'):
                lib_path=lpath
        else:
                raise NameError('Error: libiomp5.so not found: No such file or directory')'''
        f = open(Output_dir+"Progress.st", "a")
        f.write('2. libiomp5Path:'+lib_path+'\n')
        f.close()
        
def Set_paths(LigBuilder_path='',Output_dir='',libiomp5_file=''):
        if len(Ouput_dir)==0:
                Output_dir=os.getcwd()+'/'
                print('Warning:',Output_dir,' is set as the default output directory')
        else:
                if Output_dir.endswith('/'):
                        if os.path.exists(Output_dir):
                                print(Output_dir,' is set as the default output directory')
                        else:
                                raise NameError('Error:',Output_dir,'not found: No such file or directory')
                else:
                        if os.path.exists(Output_dir+'/'):
                                Output_dir=Output_dir+'/'
                                print(Output_dir,' is set as the default output directory')  
                        else:
                                raise NameError('Error:',Output_dir,'not found: No such file or directory')      
        global Out_path
        Out_path=Output_dir
        f = open(Output_dir+"Progress.st", "w")
        f.write('0. OutDir:'+Out_path+'\n')
        f.close()
        if len(LigBuilder_path)==0:
                raise NameError('Error: Please provide path to the LigBuilder installation directory')
        else:
                set_LigBuilder(LigBuilder_path,Output_dir)
        if len(libiomp5_file)==0:
                raise NameError('Error: Please provide path to the libiomp5.so file')
        else:
                set_libiomp5(libiomp5_file,Output_dir)
                
        
def unique(list1):
  
    # initialize a null list
    unique_list = []
  
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def input_structure(pdb=''):
        if os.path.exists('Progress.st'):
                with open('Progress.st')as fin:
                        for line in fin:
                                if 'LigPath' in line:
                                        global LigPath
                                        line=line.rstrip().split(':')
                                        LigPath=line[-1]
                                if 'libiomp5_path' in line:
                                        global lib_path
                                        line=line.rstrip().split(':')
                                        lib_path=line[-1]
                                if 'OutDir' in line:
                                        global Out_path
                                        line=line.rstrip().split(':')
                                        Out_path=line[-1]
        else:
                try: Out_path
                except: raise NameError('Error: Currect directory is not default output directory: please change your working directory to default output dircetory')
                try: LigPath
                except: raise NameError('Error: LigBuilderV3 path not defind: please run Set_paths function first')
                try: lib_path
                except: raise NameError('Error: libiomp5.so path not defind: please run Set_paths function first')
        global Out_dir
        global pdb_nm
        if len(pdb)==0:
                raise NameError('Error: Please provide a protein file in PDB format as input')
        else:
                if '/' in  pdb:
                        pathdata=pdb.split('/')
                        if pathdata[-1].endswith('pdb'):
                                pdb_nm=pathdata[-1].split('.')[0]
                        else:
                                raise NameError('Error: file provided is in incorrect format/extenstion')
                        shutil.copy2(pdb,Out_path)
                        pdb_path=Out_path+pathdata[-1]
                        Out_dir=Out_path
                else:
                        if os.path.exists(Out_path+pdb):
                                Out_dir=Out_path
                                pdb_path=Out_path+pdb
                                pdb_nm=pdb.split('.')[0]
                        else:
                                raise NameError('Error: ',pdb,' not found: Please provide the full path')
        if os.path.exists(LigB_path+pdb_nm):
                if os.path.exists(LigB_path+pdb_nm+'cavity.log'):
                        with open(LigB_path+pdb_nm+'cavity.log')as f:
                                lines=f.read().splitlines()
                        if len(lines)<5:
                                os.system('rm '+LigB_path+pdb_nm+'cavity.log')
                        else:
                                if lines[-4].startswith('Cavity  has done the job successfully'):
                                        f = open(Out_path+"Progress.st", "a")
                                        f.write('3. PDBPath:'+pdb_path+'\n')
                                        f.close()
                                        print(lines[-5])
                                        return
                                else:
                                        os.system('rm '+LigB_path+pdb_nm+'cavity.log')
        else:
                os.mkdir(LigB_path+pdb_nm)
        if os.path.exists(LigB_path+pdb_nm):
                shutil.copy2(pdb_path,LigB_path+pdb_nm)
        with open(LigB_path+pdb_nm+'.input','w')as fout:
                with open(LigB_path+'cavity.input')as fin:
                        for line in fin:
                                if line.startswith("RECEPTOR_FILE"):
                                        fout.write('RECEPTOR_FILE                       '+pdb_nm+'/'+pdb_nm+'.pdb\n')
                                else:
                                        fout.write(line)
        os.system('chmod 775 '+LigB_path+'cavity')
        os.system('chmod 775 '+LigB_path+'build')
        os.chdir(LigB_path)
        os.system('nohup ./cavity '+LigB_path+pdb_nm+'.input > '+pdb_nm+'cavity.log &')
        while(True):
                if os.path.exists(pdb_nm+'cavity.log'):
                        with open(pdb_nm+'cavity.log')as f:
                                lines=f.read().splitlines()
                        if len(lines)<5:
                                for d in lines:
                                        if 'rror' in d:
                                                raise NameError('\n'.join(lines))
                                        elif 'ulimit -s unlimited' in d:
                                                raise NameError(d)
                                        else:
                                                continue
                        else:
                                if lines[-4].startswith('Cavity  has done the job successfully'):
                                        break
                                else:
                                        for d in lines:
                                                if 'rror' in d:
                                                        raise NameError(d)
                                                else:
                                                        continue
        if os.path.exists(pdb_nm+'cavity.log'):
                with open(pdb_nm+'cavity.log')as f:
                        for line in f:
                                if line.startswith('Cavity will output'):
                                        print(line)
        os.chdir(Out_path)
        f = open(Out_path+"Progress.st", "a")
        f.write('3. PDBPath:'+pdb_path+'\n')
        f.close()

def cavity_view(CvID):
	global Cav_detect_ID
	Cav_detect_ID=CvID
	if os.path.exists('Progress.st'):
		with open('Progress.st')as fin:
			for line in fin:
				if 'LigPath' in line:
					global LigPath
					line=line.rstrip().split(':')
					LigPath=line[-1]
				if 'libiomp5_path' in line:
					global lib_path
					line=line.rstrip().split(':')
					lib_path=line[-1]
				if 'PDBPath' in line:
					line=line.rstrip().split(':')
					line=line[-1].split('/')
					global cwdir
					global pdb_nm
					cwdir='/'.join(line[0:-1])+'/'
					pdb_nm=line[-1].split('.')[0]
	else:
		try: LigPath
		except: raise NameError('Error: LigBuilderV3 path not defind: please run Set_paths function first')
		try: lib_path
		except: raise NameError('Error: libiomp5.so path not defind: please run Set_paths function first')
		try: pdb_nm
		except: raise NameError('Error: pdb file not provided: please run input_structure function first')
		try: Out_path
		except: 
			raise NameError('Error: Currect directory is not default output directory: please change your working directory to default output dircetory')
		else:
			cwdir=Out_path
	
	cv_ln=[]
	if os.path.exists(LigB_path+pdb_nm+"/"+pdb_nm+"_surface_"+str(CvID)+".pdb"):
		with open(LigB_path+pdb_nm+"/"+pdb_nm+"_surface_"+str(CvID)+".pdb")as fin:
			for line in fin:
				if line.startswith('HETATM'):
					cv_ln.append(line.rstrip())
	else:
		raise NameError('Cavity'+str(CvID)+' pdb file not found : please run cavity_detect function again')
	if os.path.exists(LigB_path+pdb_nm+"/"+pdb_nm+".pdb"):
		with open(cwdir+pdb_nm+"/"+pdb_nm+"_cmp.pdb",'w')as fout:
			with open(LigB_path+pdb_nm+"/"+pdb_nm+".pdb")as fin:
				for line in fin:
					if line.startswith('END'):
						fout.write('\n'.join(cv_ln)+'\nEND\n')
					else:
						fout.write(line)
	else:
		if os.path.exists(cwdir+pdb_nm+"/"+pdb_nm+".pdb"):
			with open(cwdir+pdb_nm+"/"+pdb_nm+"_cmp.pdb",'w')as fout:
				with open(cwdir+pdb_nm+"/"+pdb_nm+".pdb")as fin:
					for line in fin:
						if line.startswith('END'):
							fout.write('\n'.join(cv_ln)+'\nEND\n')
						else:
							fout.write(line)
		else:
			raise NameError('Recptor pdb file not found in the current directory')
	with open(cwdir+pdb_nm+"/"+pdb_nm+"_cmp.pdb") as ifile:
		system = "".join([x for x in ifile])
	view = py3Dmol.view(width=800, height=600)
	view.addModelsAsFrames(system)
	view.setStyle({'model': -1}, {"stick": {'color': 'spectrum'}})
	view.zoomTo()
	view.show()
	
def compund_synthesis(CavID,ligand_count=500,detected='F'):
        if os.path.exists('Progress.st'):
                with open('Progress.st')as fin:
                        for line in fin:
                                if 'LigPath' in line:
                                        global LigPath
                                        line=line.rstrip().split(':')
                                        LigPath=line[-1]
                                if 'libiomp5_path' in line:
                                        global lib_path
                                        line=line.rstrip().split(':')
                                        lib_path=line[-1]
                                if 'PDBPath' in line:
                                        line=line.rstrip().split(':')
                                        line=line[-1].split('/')
                                        global cwdir
                                        global pdb_nm
                                        cwdir='/'.join(line[0:-1])+'/'
                                        pdb_nm=line[-1].split('.')[0]
        else:
                try: LigPath
                except: raise NameError('Error: LigBuilderV3 path not defind: please run set_LigBuilder function first')
                try: lib_path
                except: raise NameError('Error: libiomp5.so path not defind: please run libiomp5_path function first')
                try: pdb_nm
                except: raise NameError('Error: pdb file not provided: please run input_structure function first')
                try: Out_path
                except: raise NameError('Error: Currect directory is not default output directory: please change your working directory to default output dircetory')
                else: cwdir=Out_path

        if detected=='F':
                f = open(cwdir+"Progress.st", "a")
                f.write('4. Cavity:'+str(CavID)+'\n')
                f.close()
        with open(LigB_path+pdb_nm+'_build.input','w')as fout:
                with open(LigB_path[:-6]+'example/build-just-example.input')as fin:
                        for line in fin:
                                if line.startswith('$TASK_NAME$'):
                                        fout.write('$TASK_NAME$                         '+pdb_nm+'\n')
                                elif line.startswith('POCKET_ATOM_FILE'):
                                        fout.write('POCKET_ATOM_FILE[1]                 '+pdb_nm+'/'+pdb_nm+'_pocket_'+str(CavID)+'.txt\n')
                                elif line.startswith('POCKET_GRID_FILE'):
                                        fout.write('POCKET_GRID_FILE[1]                 '+pdb_nm+'/'+pdb_nm+'_grid_'+str(CavID)+'.txt\n')
                                elif line.startswith('PROTEIN_FILE'):
                                        fout.write('PROTEIN_FILE[1]                     '+pdb_nm+'/'+pdb_nm+'.pdb\n')
                                elif line.startswith('SEED_LIGAND_LIST'):
                                        fout.write('SEED_LIGAND_LIST                    '+pdb_nm+'/extract/INDEX\n')
                                elif line.startswith('REPORT_FIGURE_TYPE'):
                                        fout.write('#REPORT_FIGURE_TYPE                 2\n')
                                elif line.startswith('USER_SEED_LIST'):
                                        fout.write('USER_SEED_LIST                          '+pdb_nm+'/user.list\n')
                                elif line.startswith('SEED_OUTPUT_DIRECTORY'):
                                        fout.write('SEED_OUTPUT_DIRECTORY                   '+pdb_nm+'/extract\n')
                                elif line.startswith('EXTRACT_INTERACTION_LIST'):
                                        fout.write('EXTRACT_INTERACTION_LIST                '+pdb_nm+'/inhibitor.list\n')
                                elif line.startswith('OUTPUT_EXTRACT_INTERACTION'):
                                        fout.write('OUTPUT_EXTRACT_INTERACTION              '+pdb_nm+'/key_grid_$TASK_NAME$\n')
                                else:
                                        fout.write(line)

        smi_list=[]
        with open(cwdir+'Synth.csv','w')as fout:
                fout.write('ID,SMILES\n')
        strt_ct=datetime.now()
        strt_ts = strt_ct.timestamp()
        start_ts= datetime.fromtimestamp(strt_ts)
        while(True):
                os.chdir(LigB_path)
                if os.path.exists(LigB_path+'result/'):
                        os.system('rm -r '+LigB_path+'result/')
                #os.system('export LD_LIBRARY_PATH=\"'+lib_path+':$PATH\"')
                #os.system('export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'+lib_path+'/')
                os.system('export LD_LIBRARY_PATH=\"'+lib_path+':$PATH\"\nnohup ./build -Automatic '+pdb_nm+'_build.input > build.log &')
                while(True):
                        with open('build.log', 'r') as fin:
                                for line in fin:
                                        if 'libiomp5.so: cannot open' in line:
                                                raise NameError(line+'\nPlease rerun the libiomp5_path function again !!')
                                        elif 'error' in line:
                                                raise NameError(line)
                        if os.path.exists(LigB_path+'run_'+pdb_nm+'.list'):
                                os.system('chmod 775 run_'+pdb_nm+'.list')
                                os.system('export LD_LIBRARY_PATH=\"'+lib_path+':$PATH\"\nnohup ./run_'+pdb_nm+'.list > run_'+pdb_nm+'.log &')
                                break
                while(True):
                        with open('build.log', 'r') as f:
                                lines=f.read().splitlines()
                        if len(lines)<3:
                                continue
                        if lines[-2].startswith('LigBuilder/Automatic has done the job successfully'):
                                break
                        end_ct=datetime.now()
                        end_ts = end_ct.timestamp()
                        ende_ts = datetime.fromtimestamp(end_ts)
                        delta = ende_ts - start_ts
                        if delta.total_seconds() > 14400:
                                raise NameError('Error: Cavity:'+str(CavID)+' is not an optimal target site')
                if os.path.exists(LigB_path+'result/cluster_'+pdb_nm+'/cluster_000001.mol2'):
                        cmp_c=0
                        for fis in os.listdir(LigB_path+'result/cluster_'+pdb_nm+'/'):
                                if fis.startswith('cluster_'):
                                        cmp_c+=1
                else:
                        continue
                with open(cwdir+'Synth.csv','a')as fout:
                        for i in range(1,cmp_c+1):
                                os.system('obabel -imol2 '+LigB_path+'result/cluster_'+pdb_nm+'/cluster_00000'+str(i)+'.mol2 -osmi -O'+LigB_path+'result/cluster_'+pdb_nm+'/cmp.smi')
                                with open(LigB_path+'result/cluster_'+pdb_nm+'/cmp.smi')as fin:
                                        for line in fin:
                                                line=line.split('\t')
                                                if line[0] not in smi_list:
                                                        try:
                                                                os.system('obabel -imol2 '+LigB_path+'result/cluster_'+pdb_nm+'/cluster_00000'+str(i)+'.mol2 -opdbqt -O'+cwdir+'cmp'+str(1+len(smi_list))+'.pdbqt')
                                                        except:
                                                                continue
                                                        fout.write('cmp'+str(1+len(smi_list))+','+line[0]+'\n')
                                                        smi_list.append(line[0])
                if len(smi_list)>ligand_count:
                        os.chdir(cwdir)
                        f = open("Progress.st", "a")
                        f.write('5. cavity_synth:'+str(len(smi_list))+'\n')
                        f.close()
                        if os.path.exists(LigB_path+pdb_nm):
                                st='No'
                                crds=''
                                with open(LigB_path+pdb_nm+"/"+pdb_nm+"_surface_"+str(CavID)+".pdb")as fin:
                                        for line in fin:
                                                if 'REMARK   3 Area :' in line:
                                                        st='Go'
                                                        continue
                                                else:
                                                        if st=='Go':
                                                                if line.endswith('Z\n'):
                                                                        continue
                                                                elif line.endswith('REMARK   3 \n'):
                                                                        st='No'
                                                                        continue
                                                                else:
                                                                        line=line.rstrip().split(' ')[-1]
                                                                        crds+=line+'\t'
                        os.chdir(cwdir)
                        f = open("Progress.st", "a")
                        f.write('6. cavity_grid:'+crds+'\n')
                        f.close()
                        break
        res_data=pd.read_csv(cwdir+'Synth.csv')
        return res_data

def cavity_detect(res_list):
        if os.path.exists('Progress.st'):
                with open('Progress.st')as fin:
                        for line in fin:
                                if 'LigPath' in line:
                                        global LigPath
                                        line=line.rstrip().split(':')
                                        LigPath=line[-1]
                                if 'libiomp5_path' in line:
                                        global lib_path
                                        line=line.rstrip().split(':')
                                        lib_path=line[-1]
                                if 'PDBPath' in line:
                                        line=line.rstrip().split(':')
                                        line=line[-1].split('/')
                                        global cwdir
                                        global pdb_nm
                                        cwdir='/'.join(line[0:-1])+'/'
                                        pdb_nm=line[-1].split('.')[0]
        else:
                try: LigPath
                except: raise NameError('Error: LigBuilderV3 path not defind: please run set_LigBuilder function first')
                try: lib_path
                except: raise NameError('Error: libiomp5.so path not defind: please run libiomp5_path function first')
                try: pdb_nm
                except: raise NameError('Error: pdb file not provided: please run input_structure function first')
                try: Out_path
                except: raise NameError('Error: Currect directory is not default output directory: please change your working directory to default output dircetory')
                else: cwdir=Out_path
                
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N','GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W','ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        alphasite=[]
        if os.path.exists(cwdir+res_list):
                with open(cwdir+res_list)as fin:
                        for line in fin:
                                line=line.rstrip().split("\t")
                                alphasite.append(str(int(line[-1]))+","+line[0].strip())
        else:
                raise NameError('Error: ',res_list,' not found: No such file or directory')
        cav=1
        res_percent={}
        #print("C"+str(k)+"_ymdb"+str(i)+"/FM"+str(j)+"/gpcr_ymdb"+str(i)+"a_conf"+str(j)+"th_cavity_"+str(cav)+".pdb")
        #print(alphasite)
        if os.path.exists(LigB_path+pdb_nm):
                while(os.path.exists(LigB_path+pdb_nm+"/"+pdb_nm+"_cavity_"+str(cav)+".pdb")):
                        with open(LigB_path+pdb_nm+"/"+pdb_nm+"_cavity_"+str(cav)+".pdb")as fin:
                                res_pos=[]
                                for line in fin:
                                        line=line.rstrip()
                                        if line.startswith("ATOM"):
                                                res=d[line[17:20]]
                                                pos=int(line[22:26])
                                                res_pos.append(str(pos)+","+res)
                        unique_res_pos=unique(res_pos)
                        #print(unique_res_pos)
                        ALcount=0
                        for aRES in alphasite:
                                if aRES in unique_res_pos:
                                        ALcount+=1
                        #print(ALcount)
                        res_percent[cav]=((ALcount/len(alphasite))*100)
                        cav+=1
                #print(max(res_percent.values()))
                #print(res_percent)
                #print(list(res_percent.values()).count(max(res_percent.values())))
                if list(res_percent.values()).count(max(res_percent.values())) > 1:
                        Dscore={}
                        for ind in list_duplicates_of(list(res_percent.values()),max(res_percent.values())):
                                CAVITY=list(res_percent.keys())[ind]
                                with open(LigB_path+pdb_nm+"/"+pdb_nm+"_surface_"+str(CAVITY)+".pdb")as fin:
                                        for line in fin:
                                                if 'DrugScore' in line:
                                                        line=line.rstrip()
                                                        Dscore[CAVITY]=float(line.split(' ')[-1])
                        f = open(cwdir+"Progress.st", "a")
                        f.write('4. Cavity:'+str(list(Dscore.keys())[list(Dscore.values()).index(max(Dscore.values()))])+'\n')
                        f.close()
                        os.system('echo Cav%:'+str(max(res_percent.values())))
                        return compund_synthesis(list(Dscore.keys())[list(Dscore.values()).index(max(Dscore.values()))],detected='T')
                else:
                        CAVITY=list(res_percent.keys())[list(res_percent.values()).index(max(res_percent.values()))]
                        f = open(cwdir+"Progress.st", "a")
                        f.write('4. Cavity:'+str(CAVITY)+'\n')
                        f.close()
                        os.system('echo Cav%:'+str(max(res_percent.values())))
                        return compund_synthesis(CAVITY,detected='T')
