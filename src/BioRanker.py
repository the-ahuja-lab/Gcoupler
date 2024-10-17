from signaturizer import Signaturizer
from tqdm import tqdm
import os
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import networkx as nx
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt
import re





def ends_with_integer(string):
    return bool(re.search(r'\d+$', string))

def biofet(smi,Chs):
    sig_df=pd.DataFrame()
    sig_df['smiles']=smi
    for chi in Chs:
        if ends_with_integer(chi):
            sign = Signaturizer(chi)
            results = sign.predict(smi)
            df=pd.DataFrame(results.signature)
            for clm in list(df.columns):
                df=df.rename(columns={clm:chi+'_'+str(clm)})
            sig_df=pd.concat([sig_df,df],axis = 1)
        else:
            for i in range(1,6):
                sign = Signaturizer(chi+str(i))
                results = sign.predict(smi)
                df=pd.DataFrame(results.signature)
                for clm in list(df.columns):
                    df=df.rename(columns={clm:chi+str(i)+'_'+str(clm)}) 
                sig_df=pd.concat([sig_df,df],axis = 1)
    return sig_df

def hetplot(df,D):
    df1=df.drop(['QID','smiles'],axis=1)
    df1.index=df['QID'].tolist()
    df_transposed = df1.T
    scaler = StandardScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df_transposed), index=df_transposed.index, columns=df_transposed.columns)
    df_scaled = df_scaled.T
    sns.clustermap(df_scaled, 
                   cmap='viridis',  
                   figsize=(10, 8), 
                   standard_scale=0,  # Scale the data for clustering by rows
                   method='average',  # Linkage method for clustering
                   metric='euclidean' # Distance metric
                  )
    plt.suptitle('Top Hits', y=1.02)
    plt.savefig(D+'TopHit_heatmap.pdf')
    plt.show()


def rankMatrix(smilist,HABlist,Prt,t_hit):
    fetDict={'Chem':'A','FP2D':'A1','FP3D':'A2','Sfld':'A3','StrKeys':'A4','PhysChem':'A5','Tgt':'B','MoA':'B1','MetaGns':'B2','Cry':'B3','Bnd':'B4','HTSBio':'B5','Ntwk':'C','SMRoles':'C1','SMPaths':'C2','SigPaths':'C3','BioProc':'C4','Intome':'C5','Cls':'D','Trans':'D1','CancCL':'D2','ChemGen':'D3','Morph':'D4','CellBio':'D5','Clncs':'E','ThrpA':'E1','Indctns':'E2','SEff':'E3','DisTox':'E4','DDI':'E5'}
    if isinstance(Prt, list):
        desc=[fetDict[fts] for fts in Prt]
        tempdict={k: v for k, v in zip(desc, Prt)}
    elif isinstance(Prt, int):
        if Prt==5:
            desc=['A','B','C','D','E']
        elif Prt==25:
            desc=[]
            for x in ['A','B','C','D','E']:
                for y in range(1,6):
                    desc.append(x+str(y))
            tempPrt=[]
            for value in desc:
                keys = [key for key, v in fetDict.items() if v == value]
                tempPrt.extend(keys)
            tempdict={k: v for k, v in zip(desc, tempPrt)}
        else:
            raise NameError('Error: batch calculation is only possible for 5 or 25 BioAtlas Space!')
    else:
        desc=['A','B','C','D','E']
        tempPrt=[]
        for value in desc:
            keys = [key for key, v in fetDict.items() if v == value]
            tempPrt.extend(keys)
        tempdict={k: v for k, v in zip(desc, tempPrt)}
    
    query_fet=biofet(smilist,desc)
    HAB_fet=biofet(HABlist,desc)
    
    Res5_df=pd.DataFrame()
    Res5_df['smiles']=query_fet['smiles'].tolist()
    for dsc in tqdm(desc):
        Qsub_df=query_fet.filter(regex='^'+dsc)
        Qsub_df.index=query_fet['smiles']
        Hsub_df=HAB_fet.filter(regex='^'+dsc)
        Hsub_df.index=HAB_fet['smiles']
        cos_sim = cosine_similarity(Hsub_df.values, Qsub_df.values)
        similarity_df = pd.DataFrame(cos_sim, index=[f"{i}" for i in list(Hsub_df.index)],
                                     columns=[f"{i}" for i in list(Qsub_df.index)])

        similarity_matrix = np.array(similarity_df)
        transition_matrix = similarity_matrix / similarity_matrix.sum(axis=0, keepdims=True)
        B = nx.DiGraph()

        row_labels = list(similarity_df.index)
        col_labels = list(similarity_df.columns)

        B.add_nodes_from(row_labels, bipartite=0)  # Add row nodes
        B.add_nodes_from(col_labels, bipartite=1)  # Add column nodes

        for i, row in enumerate(row_labels):
            for j, col in enumerate(col_labels):
                B.add_edge(row, col, weight=transition_matrix[i, j])

        # Calculate PageRank
        pagerank_scores = nx.pagerank(B, weight='weight')

        col_scores = {node: pagerank_scores[node] for node in col_labels}

        df = pd.DataFrame({
            'smiles': list(col_scores.keys()),
            dsc: list(col_scores.values())
        })
        Res5_df=Res5_df.merge(df,on='smiles')
        
    #top % hit
    lists=[]
    for i in desc:
        lists.append(Res5_df.sort_values(by=i, ascending=False)['smiles'].tolist()[:round(len(Res5_df)*t_hit/100)])
        
    #union
    result_list = []
    seen = set()
    for lst in lists:
        for item in lst:
            if item not in seen:
                result_list.append(item)
                seen.add(item)
                
    Sub_Res5_df=Res5_df.loc[Res5_df['smiles'].isin(result_list)]
    Sub_Res5_df.loc[:,'QID']=['Hit'+str(string+1) for string in range(len(Sub_Res5_df))]
    for clm in list(Sub_Res5_df.columns):
        if clm in tempdict.keys():
            Sub_Res5_df.rename(columns={clm:tempdict[clm]},inplace=True)
    return Sub_Res5_df

def threshold_optimize(path='',method='G-mean',mdl=''):
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
                    global TO
                    global pdb_nm
                    TO='/'.join(line[0:-1])+'/'
                    pdb_nm=line[-1].split('.')[0]
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery>
    try: TO
    except: raise NameError('Error: Could not find default output directory: Please run Synthesizer module first')
    else:
        if os.path.exists(TO):
            pass
        else:
            raise NameError('Error: Could not locate default output directory: Set directory either removed or does not exist')
    if len(mdl)==0:
        raise NameError('Error: no model selected: Please select a GNN model of choice')
    elif mdl=='GCM':
        if os.path.exists(TO+'GraphConv_train.csv'):
                pass
        else:
                raise NameError('Error: pre-trained model not found: Please run multi_model_test function first')
        tdf=pd.read_csv(TO+'GraphConv_train.csv')
        fpr, tpr, thresholds = roc_curve(tdf['GT_label'], tdf['PD_prob'])
        print(mdl+' Training AUC Score - ' + str(auc(fpr, tpr)))
        if method=='G-mean':
            gmeans = sqrt(tpr * (1-fpr))
            ix = argmax(gmeans)
            print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
        else:
            J = tpr - fpr
            ix = argmax(J)
            print('Best Threshold=%f, Youden’s J statistic=%.3f' % (thresholds[ix], J[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
    elif mdl=='AFP':
        if os.path.exists(TO+'AttentiveFP_train.csv'):
                pass
        else:
                raise NameError('Error: pre-trained model not found: Please run multi_model_test function first')
        tdf=pd.read_csv(TO+'AttentiveFP_train.csv')
        fpr, tpr, thresholds = roc_curve(tdf['GT_label'], tdf['PD_prob'])
        print(mdl+' Training AUC Score - ' + str(auc(fpr, tpr)))
        if method=='G-mean':
            gmeans = sqrt(tpr * (1-fpr))
            ix = argmax(gmeans)
            print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
        else:
            J = tpr - fpr
            ix = argmax(J)
            print('Best Threshold=%f, Youden’s J statistic=%.3f' % (thresholds[ix], J[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
    elif mdl=='GCN':
        if os.path.exists(TO+'GCNetwork_train.csv'):
                pass
        else:
                raise NameError('Error: pre-trained model not found: Please run multi_model_test function first')
        tdf=pd.read_csv(TO+'GCNetwork_train.csv')
        fpr, tpr, thresholds = roc_curve(tdf['GT_label'], tdf['PD_prob'])
        print(mdl+' Training AUC Score - ' + str(auc(fpr, tpr)))
        if method=='G-mean':
            gmeans = sqrt(tpr * (1-fpr))
            ix = argmax(gmeans)
            print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
        else:
            J = tpr - fpr
            ix = argmax(J)
            print('Best Threshold=%f, Youden’s J statistic=%.3f' % (thresholds[ix], J[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
    elif mdl=='GAT':
        if os.path.exists(TO+'GAT_train.csv'):
                pass
        else:
                raise NameError('Error: pre-trained model not found: Please run multi_model_test function first')
        tdf=pd.read_csv(TO+'GAT_train.csv')
        fpr, tpr, thresholds = roc_curve(tdf['GT_label'], tdf['PD_prob'])
        print(mdl+' Training AUC Score - ' + str(auc(fpr, tpr)))
        if method=='G-mean':
            gmeans = sqrt(tpr * (1-fpr))
            ix = argmax(gmeans)
            print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))
        else:
            J = tpr - fpr
            ix = argmax(J)
            print('Best Threshold=%f, Youden’s J statistic=%.3f' % (thresholds[ix], J[ix]))
            #print('Threshold False Positive Rate - %f, Threshold True Positive Rate - %f' % (fpr[ix], tpr[ix]))

def analyse(path='',fi='',Qdf='',threshold=0.9,Property=5,top_hit=15,out_matrix=False):
        if path.endswith('/'):
                pass
        else:
                path=path+'/'
        global TO
        global fnm
        if os.path.exists(path+'Progress.st'):
                with open(path+'Progress.st')as fin:
                        for line in fin:
                                if 'PDBPath' in line:
                                        line=line.rstrip().split(':')
                                        line=line[-1].split('/')
                                        TO='/'.join(line[0:-1])+'/'
                                        pdb_nm=line[-1].split('.')[0]
                                elif 'CMP labeled' in line:
                                        line=line.rstrip().split(':')
                                        fnm=line[-1]
        else:
                raise NameError('Path to default output folder either not provided or wrong')
                #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
        if len(fi)==0:
                try: TO
                except: raise NameError('Error: Could not find default output directory: Please run Synthesizer module first')
                if os.path.exists(TO):
                        pass
                else:
                        raise NameError('Error: Could not locate default output directory: Set directory either removed or does not exist')
                try: fnm
                except: raise NameError('Error: synthetic compound file not found: please provide the file path to the output of vina_pred function !!')
                if os.path.exists(fnm):
                        Data=pd.read_csv(fnm)
                        fnm='Labeled_cmps.csv'
                else:
                        raise NameError('Error: synthetic compound file not found: please provide the file path to the output of vina_pred function !!')
        elif os.path.exists(fi):
                if '/' in fi:
                        if fi.endswith('.csv'):
                                Data=pd.read_csv(fi)
                                fi=fi.split('/')
                                TO='/'.join(fi[0:-1])+'/'
                                fnm=fi[-1]
                                os.chdir(TO)
                        else:
                                if fi.endswith('/'):
                                        if os.path.exists(fi+'Labeled_cmps.csv'):
                                                Data=pd.read_csv(fi+'Labeled_cmps.csv')
                                                TO=fi
                                                fnm='Labeled_cmps.csv'
                                                os.chdir(fi)
                                        else:
                                                raise NameError('Error: synthetic compound file not found: please provide the file path to the output of vina_pred function !!')
                                else:
                                        if os.path.exists(fi+'/Labeled_cmps.csv'):
                                                Data=pd.read_csv(fi+'/Labeled_cmps.csv')
                                                TO=fi+'/'
                                                fnm='Labeled_cmps.csv'
                                                os.chdir(TO)
                                        else:
                                                raise NameError('Error: synthetic compound file not found: please provide the file path to the output of vina_pred function !!')
        else:
                raise NameError('Error: synthetic compound file not found: please provide the file path to the output of vina_pred function !!')

        global Dat
        Dat=Data.rename(columns={'SMILES':'smiles','Status':'status'})
        HABlist=Dat[Dat['status']==1]['smiles'].tolist()
        
        if isinstance(Qdf, list):
                Qsmilist=Qdf
        elif isinstance(Qdf, pd.DataFrame):
                if 'Probability' in list(Qdf.columns):
                        Qsmilist=Qdf[Qdf['Probability']>=threshold]['SMILES'].tolist()
                        if len(Qsmilist)>1:
                                print(len(Qsmilist),'Molecules qualified with the probability threshold',threshold)
                        else:
                                raise NameError('Error: No moclue qualified, please relax the probability threshold !')
                else:
                        raise NameError('Error: \"Probability\" column could not be found !')
        else:
                raise NameError('Error: No query molcule is provided !!')
        res1 = rankMatrix(Qsmilist,HABlist,Prt=Property,t_hit=top_hit)
        hetplot(res1,TO)
        if out_matrix:
                return res1
        else:
                return res1[['QID','smiles']]
