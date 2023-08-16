import deepchem as dc
import numpy as np
import pandas as pd
import tensorflow as tf
from deepchem.feat import MolGraphConvFeaturizer
from deepchem.feat.graph_data import GraphData
from deepchem.feat.mol_graphs import ConvMol
from deepchem.data import NumpyDataset
from imblearn.over_sampling import SMOTE
from sklearn import preprocessing
import warnings
warnings.filterwarnings("ignore")
from collections import Counter
from sklearn.metrics import accuracy_score,roc_auc_score,cohen_kappa_score,f1_score,precision_score,recall_score,roc_curve,auc
import os
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold
import imblearn
from sklearn.metrics import matthews_corrcoef
from rdkit import Chem
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split
import sys
import dgl
import matplotlib.pyplot as plt
import sklearn
import seaborn as sns
from numpy import sqrt
from numpy import argmax

np.random.seed(123)
tf.random.set_seed(123)

def MolGraphConv(df):
    featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
    mols = [Chem.MolFromSmiles(smiles) for smiles in df["smiles"]]
    features = featurizer.featurize(mols)
    dataset = dc.data.NumpyDataset(X=features, y=df["status"], ids=df["smiles"]) 
    return dataset
def ConvMol(df):
    featurizer = dc.feat.ConvMolFeaturizer()
    mols = [Chem.MolFromSmiles(smiles) for smiles in df["smiles"]]
    features = featurizer.featurize(mols)
    dataset = dc.data.NumpyDataset(X=features, y=df["status"], ids=df["smiles"])  
    return dataset
def get_labels_2(pred_test): #Getting discrete labels from probability values    
    test_pred = [] 
    test_probs= []
    for i in range(pred_test.shape[0]):
        test_probs.append(pred_test[i][0][1])
        if(pred_test[i][0][0]>pred_test[i][0][1]):
            test_pred.append(0)
        else:
            test_pred.append(1)
    return test_pred,test_probs
def get_labels(pred_test): #Getting discrete labels from probability values    
    test_pred = [] 
    test_probs= []
    for i in range(pred_test.shape[0]):
        test_probs.append(pred_test[i][1])
        if(pred_test[i][0]>pred_test[i][1]):
            test_pred.append(0)
        else:
            test_pred.append(1)
    return test_pred,test_probs
    
    
def generate_dataset_GraphData(data,feat):  
    smiles = data['smiles']
    status=data['status']
    featurizer = feat
    features = featurizer.featurize(smiles)  
    indices = [ i for i, data in enumerate(features) if type(data) is not GraphData] #indices which failed to featurize  
    return indices
def generate_dataset_ConvMol(data,feat):  
    smiles = data['smiles']
    featurizer = feat
    features = featurizer.featurize(smiles)  
    indices = [ i for i, data in enumerate(features) if type(data) is not ConvMol]   
    return indices


def model_evaluation(truth,pred,label,D):
        score={}
        print(D+" Accuracy:",metrics.accuracy_score(truth, label))
        print(D+" MCC Score:",matthews_corrcoef(truth, label))
        print(D+" F1 Score:",f1_score(truth, label))
        fpr_rf, tpr_rf, _ = roc_curve(truth,pred)
        roc_auc_rf = auc(fpr_rf, tpr_rf)
        print(D+" AUC VALUE:",roc_auc_rf)
        kappa_rf=sklearn.metrics.cohen_kappa_score(truth,label)
        print(D+" kappa Score:",kappa_rf)
        print(D+" Precision:",metrics.precision_score(truth, label))
        print(D+" Recall:",metrics.recall_score(truth, label))
        score["Accuracy:"] = metrics.accuracy_score(truth, label)
        score["MCC Score:"] = matthews_corrcoef(truth, label)
        score["F1 Score:"] = f1_score(truth, label)
        score["AUC VALUE:"] = roc_auc_rf
        score["kappa Score:"] = kappa_rf
        score["Precision:"] = metrics.precision_score(truth, label)
        score["Recall:"] = metrics.recall_score(truth, label)
        display = metrics.RocCurveDisplay(fpr=fpr_rf, tpr=tpr_rf, roc_auc=roc_auc_rf,estimator_name=D)
        display.plot()
        plt.savefig(D+'AUC_ROC.pdf')
        plt.show()
        return pd.DataFrame(score,index=[D])

def multi_model_test(path='',fi=''):
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
        Dat=Dat[['smiles','status']]

        feat = dc.feat.MolGraphConvFeaturizer() 
        while True:# remove all error producing samples
                indices = generate_dataset_GraphData(Dat,feat)
                if len(indices) > 0: #keep removing unless there are no error prone samples
                        Dat.drop(index=indices,inplace=True)
                else:
                        break
                print('After MolGraphConv screening',len(Dat))

        FROM=TO
        Dat.to_csv(FROM+'1st-level-filetred.csv',index=False)#After 1st level of filtering
        f = open(TO+"Progress.st", "a")
        f.write('9. Filtered SMILES:'+TO+'1st-level-filetred.csv\n')
        f.close()

        splitter = dc.splits.ScaffoldSplitter()
        MGFeatures=MolGraphConv(Dat)
        MGtrain, MGtest = splitter.train_test_split(MGFeatures,seed=10)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGtrain)
        MGtrain = transformer.transform(MGtrain)
        CMFeatures=ConvMol(Dat)
        CMtrain, CMtest = splitter.train_test_split(CMFeatures,seed=10)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=CMtrain)
        CMtrain = transformer.transform(CMtrain)

        All_score=pd.DataFrame()

        deep_model='GraphConv'
        CMmodel = dc.models.GraphConvModel(mode='classification',n_tasks=1,model_dir=TO+deep_model)
        CMmodel.fit(CMtrain, nb_epoch=50)
        train_probs = CMmodel.predict(CMtrain)
        train_labels,train_preds = get_labels_2(train_probs)
        df1=model_evaluation(CMtrain.y,train_preds,train_labels,deep_model+' Train')
        test_probs = CMmodel.predict(CMtest)
        test_labels,test_preds = get_labels_2(test_probs)
        df2=model_evaluation(CMtest.y,test_preds,test_labels,deep_model+' Test')
        tdf=pd.DataFrame()
        tdf['GT_label']=CMtrain.y
        tdf['PD_label']=train_labels
        tdf['PD_prob']=train_preds
        tdf.to_csv(deep_model+'_train.csv')
        
        
        deep_model='AttentiveFP'
        AtFPmodel = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',model_dir=TO+deep_model)
        AtFPmodel.fit(MGtrain,nb_epoch=100)
        train_probs = AtFPmodel.predict(MGtrain)
        train_labels,train_preds = get_labels(train_probs)
        df3=model_evaluation(MGtrain.y,train_preds,train_labels,deep_model+' Train')
        test_probs = AtFPmodel.predict(MGtest)
        test_labels,test_preds = get_labels(test_probs)
        df4=model_evaluation(MGtest.y,test_preds,test_labels,deep_model+' Test')
        tdf=pd.DataFrame()
        tdf['GT_label']=MGtrain.y
        tdf['PD_label']=train_labels
        tdf['PD_prob']=train_preds
        tdf.to_csv(deep_model+'_train.csv')

        deep_model='GCNetwork'
        GCNmodel = dc.models.GCNModel(n_tasks=1,mode='classification',model_dir=TO+deep_model)
        GCNmodel.fit(MGtrain,nb_epoch=100)
        train_probs = GCNmodel.predict(MGtrain)
        train_labels,train_preds = get_labels(train_probs)
        df5=model_evaluation(MGtrain.y,train_preds,train_labels,deep_model+' Train')
        test_probs = GCNmodel.predict(MGtest)
        test_labels,test_preds = get_labels(test_probs)
        df6=model_evaluation(MGtest.y,test_preds,test_labels,deep_model+' Test')
        tdf=pd.DataFrame()
        tdf['GT_label']=MGtrain.y
        tdf['PD_label']=train_labels
        tdf['PD_prob']=train_preds
        tdf.to_csv(deep_model+'_train.csv')

        deep_model='GAT'
        GATmodel = dc.models.GATModel(mode='classification', n_tasks=1,model_dir=TO+deep_model)
        GATmodel.fit(MGtrain, nb_epoch=5)
        train_probs = GATmodel.predict(MGtrain)
        train_labels,train_preds = get_labels(train_probs)
        df7=model_evaluation(MGtrain.y,train_preds,train_labels,deep_model+' Train')
        test_probs = GATmodel.predict(MGtest)
        test_labels,test_preds = get_labels(test_probs)
        df8=model_evaluation(MGtest.y,test_preds,test_labels,deep_model+' Test')
        tdf=pd.DataFrame()
        tdf['GT_label']=MGtrain.y
        tdf['PD_label']=train_labels
        tdf['PD_prob']=train_preds
        tdf.to_csv(deep_model+'_train.csv')
        plt.clf()

        df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8])
        ax=sns.heatmap(df, annot=True)
        plt.savefig('Base_model_heatmap.pdf', format='pdf')
        df.to_csv(TO+'Base_model_stats.csv')
        f = open(TO+"Progress.st", "a")
        f.write('10. Base mds:'+TO+'Base_model_stats.csv\n')
        f.close()
        return df
        
def model_hpt(path='',mdl='',params={}):
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
                elif 'Filtered SMILES' in line:
                    line=line.rstrip().split(':')[-1]
                    Dat=pd.read_csv(line)
                elif 'CMP labeled' in line:
                    line=line.rstrip().split(':')
                    fnm=line[-1]
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try: TO
    except: raise NameError('Error: Could not find default output directory: Please run Synthesizer module first')
    else:
        if os.path.exists(TO):
            pass
        else:
            raise NameError('Error: Could not locate default output directory: Set directory either removed or does not exist')
        try: Dat
        except:
            try: fnm
            except: raise NameError('Error: cannot find interaction information: please run synth_classify function first')
            Data=pd.read_csv(fnm)
            Dat=Data.rename(columns={'SMILES':'smiles','Status':'status'})
            Dat=Dat[['smiles','status']]
            feat = dc.feat.MolGraphConvFeaturizer() 
            while True:# remove all error producing samples
                indices = generate_dataset_GraphData(Dat,feat)
                if len(indices) > 0: #keep removing unless there are no error prone samples
                    Dat.drop(index=indices,inplace=True)
                else:
                    break

            FROM=TO
            Dat.to_csv(FROM+'1st-level-filetred.csv',index=False)#After 1st level of filtering
            f = open(TO+"Progress.st", "a")
            f.write('9. Filtered SMILES:'+TO+'1st-level-filetred.csv\n')
            f.close()
    splitter = dc.splits.ScaffoldSplitter()
    MGFeatures=MolGraphConv(Dat)
    MGtrain, MGtest = splitter.train_test_split(MGFeatures,seed=10)
    transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGtrain)
    MGtrain = transformer.transform(MGtrain)
    CMFeatures=ConvMol(Dat)
    CMtrain, CMtest = splitter.train_test_split(CMFeatures,seed=10)
    transformer = dc.trans.DuplicateBalancingTransformer(dataset=CMtrain)
    CMtrain = transformer.transform(CMtrain)
    if os.path.exists(TO+'HPT'):
        pass
    else:
        os.mkdir(TO+'HPT')
    default_GraphConv={
        'graph_conv_layers':[64, 64],
        'dense_layer_size':128,
        'dropout':0.0,
        'number_atom_features':75,
        'batch_size':100
    }
    default_AttentiveFP={
        'num_layers':2,
        'num_timesteps':2,
        'graph_feat_size':200,
        'dropout':0
    }
    default_GCNetwork={
        'batch_size' : 16, 
        'graph_conv_layers' : [64, 64], 
        'predictor_hidden_feats' : 128, 
        'learning_rate' : 0.01,
        'predictor_droput' : 0
    }
    default_GAT={
        'n_attention_heads':8,
        'dropout':0,
        'alpha':0.2
    }
    if len(mdl)==0:
        raise NameError('Error: No model selected: please provide the name of the model to proceed with hyper parameter tuning')
    elif mdl=='GCM':
        deep_model='HPT/GraphConv'
        AtFPmetric = dc.metrics.Metric(dc.metrics.roc_auc_score)
        # defining `model_builder`
        def model_builder(**model_params):
            number_atom_features=model_params['number_atom_features']
            graph_conv_layers=model_params['graph_conv_layers']
            dropout=model_params['dropout']
            batch_size=model_params['batch_size']
            dense_layer_size=model_params['dense_layer_size']
            return dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=number_atom_features,graph_conv_layers=graph_conv_layers,dropout=dropout,batch_size=batch_size,dense_layer_size=dense_layer_size,model_dir=TO+deep_model)
        # the parameters which are to be optimized
        if len(params)==0:
                params = {
                'number_atom_features':list(range(0,100,25)),
                'graph_conv_layers':[[32,32],[64,64],[128,128]],
                'dropout':[0, 0.1, 0.5],
                'batch_size':list(range(100,200,50)),
                'dense_layer_size':list(range(120,140,2))
                }
        else:
                pass
        # Creating optimizer and searching over hyperparameters
        optimizer = dc.hyper.GridHyperparamOpt(model_builder)
        best_model, best_hyperparams, all_results =   optimizer.hyperparam_search(params, CMtrain, CMtest, AtFPmetric)
        CMmodel = dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=best_hyperparams['number_atom_features'],graph_conv_layers=best_hyperparams['graph_conv_layers'],batch_size=best_hyperparams['batch_size'],dropout=best_hyperparams['dropout'],dense_layer_size=best_hyperparams['dense_layer_size'],model_dir=TO+deep_model)
        CMmodel.fit(CMtrain,nb_epoch=100)
        CMmodel2 = dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=default_GraphConv['number_atom_features'],graph_conv_layers=default_GraphConv['graph_conv_layers'],batch_size=default_GraphConv['batch_size'],dropout=default_GraphConv['dropout'],dense_layer_size=default_GraphConv['dense_layer_size'],model_dir=TO+deep_model)
        CMmodel2.fit(CMtrain,nb_epoch=100)
        if CMmodel.evaluate(CMtest, [AtFPmetric])['roc_auc_score']> CMmodel2.evaluate(CMtest, [AtFPmetric])['roc_auc_score']:
            print(best_hyperparams)
            return best_hyperparams
        else:
            print(default_GraphConv)
            return default_GraphConv
    elif mdl=='AFP':
        deep_model='HPT/AttentiveFP'
        AtFPmetric = dc.metrics.Metric(dc.metrics.roc_auc_score)
        # defining `model_builder`
        def model_builder(**model_params):
            num_layers=model_params['num_layers']
            num_timesteps=model_params['num_timesteps']
            graph_feat_size=model_params['graph_feat_size']
            dropout=model_params['dropout']
            return dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=num_layers,num_timesteps=num_timesteps,graph_feat_size=graph_feat_size,dropout=dropout,model_dir=TO+deep_model)
        # the parameters which are to be optimized
        if len(params)==0:
                params = {
                'num_layers':list(range(5,10)),
                'num_timesteps':list(range(5,10)),
                'graph_feat_size':list(range(150,250,25)),
                'dropout':[0, 0.1, 0.5]
                }
        else:
                pass
        # Creating optimizer and searching over hyperparameters
        optimizer = dc.hyper.GridHyperparamOpt(model_builder)
        best_model, best_hyperparams, all_results =   optimizer.hyperparam_search(params, MGtrain, MGtest, AtFPmetric)
        AtFPmodel = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=best_hyperparams['num_layers'],num_timesteps=best_hyperparams['num_timesteps'],graph_feat_size=best_hyperparams['graph_feat_size'],dropout=best_hyperparams['dropout'],model_dir=TO+deep_model)
        AtFPmodel.fit(MGtrain,nb_epoch=100)
        AtFPmodel2 = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=default_AttentiveFP['num_layers'],num_timesteps=default_AttentiveFP['num_timesteps'],graph_feat_size=default_AttentiveFP['graph_feat_size'],dropout=default_AttentiveFP['dropout'],model_dir=TO+deep_model)
        AtFPmodel2.fit(MGtrain,nb_epoch=100)
        if AtFPmodel.evaluate(MGtest, [AtFPmetric])['roc_auc_score']> AtFPmodel2.evaluate(MGtest, [AtFPmetric])['roc_auc_score']:
            print(best_hyperparams)
            return best_hyperparams
        else:
            print(default_AttentiveFP)
            return default_AttentiveFP
            
    elif mdl=='GCN':
        deep_model='HPT/GCNetwork'
        AtFPmetric = dc.metrics.Metric(dc.metrics.roc_auc_score)
        # defining `model_builder`
        def model_builder(**model_params):
            batch_size=model_params['batch_size']
            graph_conv_layers=model_params['graph_conv_layers']
            predictor_hidden_feats=model_params['predictor_hidden_feats']
            learning_rate=model_params['learning_rate']
            predictor_droput=model_params['predictor_droput']
            return dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=batch_size,graph_conv_layers=graph_conv_layers,predictor_hidden_feats=predictor_hidden_feats,learning_rate=learning_rate,predictor_droput=predictor_droput,model_dir=TO+deep_model)
        # the parameters which are to be optimized
        if len(params)==0:
                params = {
                'batch_size' : list(range(10,30,2)), 
                'graph_conv_layers' : [[32,32],[64, 64],[128,128]], 
                'predictor_hidden_feats' : list(range(100,300,50)), 
                'learning_rate' : [0.01,0.1,1.0],
                'predictor_droput' : [0]
                }
        else:
                pass
        # Creating optimizer and searching over hyperparameters
        optimizer = dc.hyper.GridHyperparamOpt(model_builder)
        best_model, best_hyperparams, all_results =   optimizer.hyperparam_search(params, MGtrain, MGtest, AtFPmetric)
        GCNmodel = dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=best_hyperparams['batch_size'],graph_conv_layers=best_hyperparams['graph_conv_layers'],predictor_hidden_feats=best_hyperparams['predictor_hidden_feats'],learning_rate=best_hyperparams['learning_rate'],predictor_droput=best_hyperparams['predictor_droput'],model_dir=TO+deep_model)
        GCNmodel.fit(MGtrain,nb_epoch=100)
        GCNmodel2 = dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=default_GCNetwork['batch_size'],graph_conv_layers=default_GCNetwork['graph_conv_layers'],predictor_hidden_feats=default_GCNetwork['predictor_hidden_feats'],learning_rate=default_GCNetwork['learning_rate'],predictor_droput=default_GCNetwork['predictor_droput'],model_dir=TO+deep_model)
        GCNmodel2.fit(MGtrain,nb_epoch=100)
        if GCNmodel.evaluate(MGtest, [AtFPmetric])['roc_auc_score']> GCNmodel2.evaluate(MGtest, [AtFPmetric])['roc_auc_score']:
            print(best_hyperparams)
            return best_hyperparams
        else:
            print(default_GCNetwork)
            return default_GCNetwork
    elif mdl=='GAT':
        deep_model='HPT/GAT'
        AtFPmetric = dc.metrics.Metric(dc.metrics.roc_auc_score)
        # defining `model_builder`
        def model_builder(**model_params):
            alpha=model_params['alpha']
            n_attention_heads=model_params['n_attention_heads']
            dropout=model_params['dropout']
            return dc.models.GATModel(n_tasks=1,mode='classification',alpha=alpha,n_attention_heads=n_attention_heads,dropout=dropout,model_dir=TO+deep_model)
        # the parameters which are to be optimized
        if len(params)==0:
                params = {
                'alpha':[0.1,0.2,0.4],
                'n_attention_heads':list(range(2,20,2)),
                'dropout':[0, 0.1, 0.5]
                }
        else:
                pass
        # Creating optimizer and searching over hyperparameters
        optimizer = dc.hyper.GridHyperparamOpt(model_builder)
        best_model, best_hyperparams, all_results =   optimizer.hyperparam_search(params, MGtrain, MGtest, AtFPmetric)
        GATmodel = dc.models.GATModel(n_tasks=1,mode='classification',alpha=best_hyperparams['alpha'],n_attention_heads=best_hyperparams['n_attention_heads'],dropout=best_hyperparams['dropout'],model_dir=TO+deep_model)
        GATmodel.fit(MGtrain,nb_epoch=100)
        GATmodel2 = dc.models.GATModel(n_tasks=1,mode='classification',alpha=default_GAT['alpha'],n_attention_heads=default_GAT['n_attention_heads'],dropout=default_GAT['dropout'],model_dir=TO+deep_model)
        GATmodel2.fit(MGtrain,nb_epoch=100)
        if GATmodel.evaluate(MGtest, [AtFPmetric])['roc_auc_score']> GATmodel2.evaluate(MGtest, [AtFPmetric])['roc_auc_score']:
            print(best_hyperparams)
            return best_hyperparams
        else:
            print(default_GAT)
            return default_GAT
    else:
        raise NameError('Error: No model selected: please provide a valid model name to proceed with hyper parameter tuning')
        

        
def MD_kfold(path='',mdl='',params={},k=3):
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
                elif 'Filtered SMILES' in line:
                    line=line.rstrip().split(':')[-1]
                    Dat=pd.read_csv(line)
                elif 'CMP labeled' in line:
                    line=line.rstrip().split(':')
                    fnm=line[-1]
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try: TO
    except: raise NameError('Error: Could not find default output directory: Please run Synthesizer module first')
    else:
        if os.path.exists(TO):
            pass
        else:
            raise NameError('Error: Could not locate default output directory: Set directory either removed or does not exist')
        try: Dat
        except:
            try: fnm
            except: raise NameError('Error: cannot find interaction information: please run synth_classify function first')
            Data=pd.read_csv(fnm)
            Dat=Data.rename(columns={'SMILES':'smiles','Status':'status'})
            Dat=Dat[['smiles','status']]
            feat = dc.feat.MolGraphConvFeaturizer() 
            while True:# remove all error producing samples
                indices = generate_dataset_GraphData(Dat,feat)
                if len(indices) > 0: #keep removing unless there are no error prone samples
                    Dat.drop(index=indices,inplace=True)
                else:
                    break

            FROM=TO
            Dat.to_csv(FROM+'1st-level-filetred.csv',index=False)#After 1st level of filtering
            f = open(path+"Progress.st", "a")
            f.write('9. Filtered SMILES:'+TO+'1st-level-filetred.csv\n')
            f.close()
    best_hpt=model_hpt(path,mdl,params=params)
    f = open(TO+"Progress.st", "a")
    f.write('11. Mdl HPTs:'+mdl+'%'+str(k)+'%'+','.join(list(best_hpt.keys()))+'%'+str(list(best_hpt.values()))[1:-1]+'\n')
    f.close()
    MGFeatures=MolGraphConv(Dat)
    CMFeatures=ConvMol(Dat)
    if len(mdl)==0:
        raise NameError('Error: No model selected: please provide the name of the model to proceed with hyper parameter tuning')
    elif mdl=='GCM':
        if os.path.exists(str(k)+'FoldCV'):
            pass
        else:
            os.mkdir(str(k)+'FoldCV')
        splitter = dc.splits.RandomSplitter()
        cv=splitter.k_fold_split(CMFeatures,k,seed=10)
        scrs_trn=[]
        scrs_tst=[]
        for i in range(k):
            CMtrain, CMtest = cv[i][0],cv[i][1]
            transformer = dc.trans.DuplicateBalancingTransformer(dataset=CMtrain)
            CMtrain = transformer.transform(CMtrain)
            deep_model=str(k)+'FoldCV/GraphConv_CV_Fold'+str(i)
            CMmodel = dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=best_hpt['number_atom_features'],graph_conv_layers=best_hpt['graph_conv_layers'],batch_size=best_hpt['batch_size'],dropout=best_hpt['dropout'],dense_layer_size=best_hpt['dense_layer_size'],model_dir=TO+deep_model)
            CMmodel.fit(CMtrain,nb_epoch=100)
            train_probs = CMmodel.predict(CMtrain)
            train_labels,train_preds = get_labels_2(train_probs)
            scrs_trn.append(model_evaluation(CMtrain.y,train_preds,train_labels,deep_model))
            test_probs = CMmodel.predict(CMtest)
            test_labels,test_preds = get_labels_2(test_probs)
            scrs_tst.append(model_evaluation(CMtest.y,test_preds,test_labels,deep_model))
        FdfTrain=pd.concat(scrs_trn)
        FdfTrain.to_csv(TO+mdl+"_"+str(k)+"F_CV_Train_scores.csv")
        long_FdfTrain=pd.melt(FdfTrain)
        long_FdfTrain=long_FdfTrain.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTrain)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Train.pdf')
        plt.show()
        FdfTest=pd.concat(scrs_tst)
        FdfTest.to_csv(TO+mdl+"_"+str(k)+"F_CV_Test_scores.csv")
        long_FdfTest=pd.melt(FdfTest)
        long_FdfTest=long_FdfTest.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTest)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Test.pdf')
        plt.show()
        deep_model='GraphConv_100'
        CMmodel100 = dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=best_hpt['number_atom_features'],graph_conv_layers=best_hpt['graph_conv_layers'],batch_size=best_hpt['batch_size'],dropout=best_hpt['dropout'],dense_layer_size=best_hpt['dense_layer_size'],model_dir=TO+deep_model)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=CMFeatures)
        CMFeatures = transformer.transform(CMFeatures)
        CMmodel100.fit(CMFeatures,nb_epoch=100)
        f = open(TO+"Progress.st", "a")
        f.write('12. Final Mdl:'+TO+deep_model+'\n')
        f.close()
        return CMmodel100
    elif mdl=='AFP':
        if os.path.exists(str(k)+'FoldCV'):
            pass
        else:
            os.mkdir(str(k)+'FoldCV')
        splitter = dc.splits.RandomSplitter()
        cv=splitter.k_fold_split(MGFeatures,k,seed=10)
        scrs_trn=[]
        scrs_tst=[]
        for i in range(k):
            MGtrain, MGtest = cv[i][0],cv[i][1]
            transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGtrain)
            MGtrain = transformer.transform(MGtrain)
            deep_model=str(k)+'FoldCV/AttentiveFP_CV_Fold'+str(i)
            AtFPmodel = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=best_hpt['num_layers'],num_timesteps=best_hpt['num_timesteps'],graph_feat_size=best_hpt['graph_feat_size'],dropout=best_hpt['dropout'],model_dir=TO+deep_model)
            AtFPmodel.fit(MGtrain,nb_epoch=100)
            train_probs = AtFPmodel.predict(MGtrain)
            train_labels,train_preds = get_labels(train_probs)
            scrs_trn.append(model_evaluation(MGtrain.y,train_preds,train_labels,deep_model))
            test_probs = AtFPmodel.predict(MGtest)
            test_labels,test_preds = get_labels(test_probs)
            scrs_tst.append(model_evaluation(MGtest.y,test_preds,test_labels,deep_model))
        FdfTrain=pd.concat(scrs_trn)
        FdfTrain.to_csv(TO+mdl+"_"+str(k)+"F_CV_Train_scores.csv")
        long_FdfTrain=pd.melt(FdfTrain)
        long_FdfTrain=long_FdfTrain.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTrain)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Train.pdf')
        plt.show()
        FdfTest=pd.concat(scrs_tst)
        FdfTest.to_csv(TO+mdl+"_"+str(k)+"F_CV_Test_scores.csv")
        long_FdfTest=pd.melt(FdfTest)
        long_FdfTest=long_FdfTest.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTest)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Test.pdf')
        plt.show()
        deep_model='AttentiveFP_100'
        AtFPmodel100 = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=best_hpt['num_layers'],num_timesteps=best_hpt['num_timesteps'],graph_feat_size=best_hpt['graph_feat_size'],dropout=best_hpt['dropout'],model_dir=TO+deep_model)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGFeatures)
        MGFeatures = transformer.transform(MGFeatures)
        AtFPmodel100.fit(MGFeatures,nb_epoch=100)
        f = open(TO+"Progress.st", "a")
        f.write('12. Final Mdl:'+TO+deep_model+'\n')
        f.close()
        return AtFPmodel100
    elif mdl=='GCN':
        if os.path.exists(str(k)+'FoldCV'):
            pass
        else:
            os.mkdir(str(k)+'FoldCV')
        splitter = dc.splits.RandomSplitter()
        cv=splitter.k_fold_split(MGFeatures,k,seed=10)
        scrs_trn=[]
        scrs_tst=[]
        for i in range(k):
            MGtrain, MGtest = cv[i][0],cv[i][1]
            transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGtrain)
            MGtrain = transformer.transform(MGtrain)
            deep_model=str(k)+'FoldCV/GCNetwork_CV_Fold'+str(i)
            GCNmodel = dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=best_hpt['batch_size'],graph_conv_layers=best_hpt['graph_conv_layers'],predictor_hidden_feats=best_hpt['predictor_hidden_feats'],learning_rate=best_hpt['learning_rate'],predictor_droput=best_hpt['predictor_droput'],model_dir=TO+deep_model)
            GCNmodel.fit(MGtrain,nb_epoch=100)
            train_probs = GCNmodel.predict(MGtrain)
            train_labels,train_preds = get_labels(train_probs)
            scrs_trn.append(model_evaluation(MGtrain.y,train_preds,train_labels,deep_model))
            test_probs = GCNmodel.predict(MGtest)
            test_labels,test_preds = get_labels(test_probs)
            scrs_tst.append(model_evaluation(MGtest.y,test_preds,test_labels,deep_model))
        FdfTrain=pd.concat(scrs_trn)
        FdfTrain.to_csv(TO+mdl+"_"+str(k)+"F_CV_Train_scores.csv")
        long_FdfTrain=pd.melt(FdfTrain)
        long_FdfTrain=long_FdfTrain.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTrain)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Train.pdf')
        plt.show()
        FdfTest=pd.concat(scrs_tst)
        FdfTest.to_csv(TO+mdl+"_"+str(k)+"F_CV_Test_scores.csv")
        long_FdfTest=pd.melt(FdfTest)
        long_FdfTest=long_FdfTest.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTest)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Test.pdf')
        plt.show()
        deep_model='GCNetwork_100'
        GCNmodel100 = dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=best_hpt['batch_size'],graph_conv_layers=best_hpt['graph_conv_layers'],predictor_hidden_feats=best_hpt['predictor_hidden_feats'],learning_rate=best_hpt['learning_rate'],predictor_droput=best_hpt['predictor_droput'],model_dir=TO+deep_model)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGFeatures)
        MGFeatures = transformer.transform(MGFeatures)
        GCNmodel100.fit(MGFeatures,nb_epoch=100)
        f = open(TO+"Progress.st", "a")
        f.write('12. Final Mdl:'+TO+deep_model+'\n')
        f.close()
        return GCNmodel100
    elif mdl=='GAT':
        if os.path.exists(str(k)+'FoldCV'):
            pass
        else:
            os.mkdir(str(k)+'FoldCV')
        splitter = dc.splits.RandomSplitter()
        cv=splitter.k_fold_split(MGFeatures,k,seed=10)
        scrs_trn=[]
        scrs_tst=[]
        for i in range(k):
            MGtrain, MGtest = cv[i][0],cv[i][1]
            transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGtrain)
            MGtrain = transformer.transform(MGtrain)
            deep_model=str(k)+'FoldCV/GAT_CV_Fold'+str(i)
            GATmodel = dc.models.GATModel(n_tasks=1,mode='classification',alpha=best_hpt['alpha'],n_attention_heads=best_hpt['n_attention_heads'],dropout=best_hpt['dropout'],model_dir=TO+deep_model)
            GATmodel.fit(MGtrain,nb_epoch=100)
            train_probs = GATmodel.predict(MGtrain)
            train_labels,train_preds = get_labels(train_probs)
            scrs_trn.append(model_evaluation(MGtrain.y,train_preds,train_labels,deep_model))
            test_probs = GATmodel.predict(MGtest)
            test_labels,test_preds = get_labels(test_probs)
            scrs_tst.append(model_evaluation(MGtest.y,test_preds,test_labels,deep_model))
        FdfTrain=pd.concat(scrs_trn)
        FdfTrain.to_csv(TO+mdl+"_"+str(k)+"F_CV_Train_scores.csv")
        long_FdfTrain=pd.melt(FdfTrain)
        long_FdfTrain=long_FdfTrain.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTrain)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Train.pdf')
        plt.show()
        FdfTest=pd.concat(scrs_tst)
        FdfTest.to_csv(TO+mdl+"_"+str(k)+"F_CV_Test_scores.csv")
        long_FdfTest=pd.melt(FdfTest)
        long_FdfTest=long_FdfTest.rename(columns={'variable':'Model parameters','value':'Probability'})
        ax = sns.boxplot(x="Model parameters", y="Probability", data=long_FdfTest)
        ax.set_ylim([0, 1]) 
        ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
        plt.savefig(TO+mdl+'_'+str(k)+'Fold_CV_Test.pdf')
        plt.show()
        deep_model='GAT_100'
        GATmodel100 = dc.models.GATModel(n_tasks=1,mode='classification',alpha=best_hpt['alpha'],n_attention_heads=best_hpt['n_attention_heads'],dropout=best_hpt['dropout'],model_dir=TO+deep_model)
        transformer = dc.trans.DuplicateBalancingTransformer(dataset=MGFeatures)
        MGFeatures = transformer.transform(MGFeatures)
        GATmodel100.fit(MGFeatures,nb_epoch=100)
        f = open(TO+"Progress.st", "a")
        f.write('12. Final Mdl:'+TO+deep_model+'\n')
        f.close()
        return GATmodel100
    else:
        raise NameError('Error: No model selected: please provide a valid model name to proceed with hyper parameter tuning')
        
        
def MD_pred(path='',smi_list=[],):
    if path.endswith('/'):
        pass
    else:
        path=path+'/'
    if len(smi_list)==0:
        raise NameError('Error: no molecules provided for prediction')
    else:
        query=pd.DataFrame()
        query['smiles']=smi_list
        query['status']=['NA']*len(smi_list)
        feat = dc.feat.MolGraphConvFeaturizer() 
        while True:# remove all error producing samples
            indices = generate_dataset_GraphData(query,feat)
            if len(indices) > 0: #keep removing unless there are no error prone samples
                query.drop(index=indices,inplace=True)
            else:
                break
        print(len(query),'query molecules passed the GraphConversion screening')
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
                elif 'Final Mdl' in line:
                    line=line.rstrip().split(':')
                    restore_path=line[-1]
                elif 'Mdl HPTs' in line:
                    line=line.rstrip().split(':')[-1]
                    data=line.split('%')
                    mdl=data[0]
                    vals=data[3]
                    kes=data[2].split(',')
                    vals=[]
                    for i in data[3].split(','):
                        if '[' in i:
                            sub_l=[]
                            sub_l.append(int(i.split('[')[-1]))
                        elif ']' in i:
                            sub_l.append(int(i.split(']')[0]))
                            vals.append(sub_l)
                            sub_l=[]
                        elif '.' in i:
                            vals.append(float(i))
                        else:
                            vals.append(int(i))
                    hpts={}
                    if len(kes)==len(vals):
                        for i in range(len(kes)):
                            hpts[kes[i]]=vals[i]
                    else:
                        raise NameError('Model Hyperparameter discrepancy encountered!')
    else:
        raise NameError('Path to default output folder either not provided or wrong')
        #raise NameError('Error: current directory does not contain necessary files for the function to operate: please set the currect working directory to default output directoery of Synthesizer module')
    try: TO
    except: raise NameError('Error: Could not find default output directory: Please run Synthesizer module first')
    else:
        if os.path.exists(TO):
            pass
        else:
            raise NameError('Error: Could not locate default output directory: Set directory either removed or does not exist')
        try: restore_path
        except: raise NameError('Error: no pre-trained model path found: please run MD_kfold function first')
        else:
            if os.path.exists(restore_path):
                pass
            else:
                raise NameError('Error: no pre-trained model for the selected option is found: please run MD_kfold function first')
        try: hpts
        except: raise NameError('Error: Cross validated hyper parameters not found: please run MD_kfold function first')
    if len(mdl)==0:
        raise NameError('Error: no pre-trained model path/hyper parameters found: please run MD_kfold function first')
    elif mdl=='GCM':
        CMQuery=ConvMol(query)
        CMmodel100 = dc.models.GraphConvModel(n_tasks=1,mode='classification',number_atom_features=hpts['number_atom_features'],graph_conv_layers=hpts['graph_conv_layers'],batch_size=hpts['batch_size'],dropout=hpts['dropout'],dense_layer_size=hpts['dense_layer_size'],model_dir=restore_path)
        CMmodel100.restore()
        probs = CMmodel100.predict(CMQuery)
        labels,preds = get_labels_2(probs)
        query['status']=labels
        query['Probability']=preds
        query=query.rename(columns={'smiles':'SMILES','status':'Status'})
        return query
    elif mdl=='AFP':
        MGQuery=MolGraphConv(query)
        AtFPmodel100 = dc.models.AttentiveFPModel(n_tasks=1,mode='classification',num_layers=hpts['num_layers'],num_timesteps=hpts['num_timesteps'],graph_feat_size=hpts['graph_feat_size'],dropout=hpts['dropout'],model_dir=restore_path)
        AtFPmodel100.restore()
        probs = AtFPmodel100.predict(MGQuery)
        labels,preds = get_labels(probs)
        query['status']=labels
        query['Probability']=preds
        query=query.rename(columns={'smiles':'SMILES','status':'Status'})
        return query
    elif mdl=='GCN':
        MGQuery=MolGraphConv(query)
        GCNmodel100 = dc.models.GCNModel(n_tasks=1,mode='classification',batch_size=hpts['batch_size'],graph_conv_layers=hpts['graph_conv_layers'],predictor_hidden_feats=hpts['predictor_hidden_feats'],learning_rate=hpts['learning_rate'],predictor_droput=hpts['predictor_droput'],model_dir=restore_path)
        GCNmodel100.restore()
        probs = GCNmodel100.predict(MGQuery)
        labels,preds = get_labels(probs)
        query['status']=labels
        query['Probability']=preds
        query=query.rename(columns={'smiles':'SMILES','status':'Status'})
        return query
    elif mdl=='GAT':
        MGQuery=MolGraphConv(query)
        GATmodel100 = dc.models.GATModel(n_tasks=1,mode='classification',alpha=hpts['alpha'],n_attention_heads=hpts['n_attention_heads'],dropout=hpts['dropout'],model_dir=restore_path)
        GATmodel100.restore()
        probs = GATmodel100.predict(MGQuery)
        labels,preds = get_labels(probs)
        query['status']=labels
        query['Probability']=preds
        query=query.rename(columns={'smiles':'SMILES','status':'Status'})
        return query
    else:
        raise NameError('Error: No model selected: please provide a valid model name to proceed with hyper parameter tuning')


def cutoff_optimize(path='',method='G-mean',mdl=''):
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

