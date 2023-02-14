# Interactome Prediction using Gcoupler
 <br>
<div align="center">
<img src="Images/Overview.png"></div>
<br>

### Introduction

Gcoupler, which leverages an integrative approach combining _de novo_ ligand design, statistical methods, and Graph Neural Networks for rational prediction of high-affinity ligands. Gcoupler offers an efficient and comparatively faster route to explore endogenous allosteric sites of GPCRs, including the GPCR-Gα interface. <br/><br/>


The only strong dependencies for this resource are [**RDKit**](https://www.rdkit.org/) and [**DeepChem**](https://github.com/deepchem/deepchem) which can be installed in a local [Conda](https://conda.io/) environment.


Click [here]() to download a copy of pre-complied conda environment yml file. The environment can be setup by 
```
$ conda env create -f Gcoupler.yml
$ conda activate Gcoupler
```

**Third-party dependencies**
1. [LigBuilder V3.0]()
2. [OpenBabel (2.4.1)]()

The installation procedure takes less than 5 minutes.

**Installation of LigBuilderV3**
```
$ tar –xvzf LigBuilderV3.0.tar.gz
$ cd LigBuilderV3.0/
$ chmod 775 configure
$ ./configure
```
**Installation of OpenBabel**
```
$ tar –xvzf openbabel-2.4.1.tar.gz
$ mkdir build
$ cd build
$ cmake ..
$ make -j2
$ sudo make install
```

## How to use Gcoupler?


### Installation using pip 
```
$ pip install Gcoupler
```

## Pipeline
Gcoupler supports 3 distinct modules:<br/>
1. Synthesizer
2. Authenticator
3. Generator

### Synthesizer

To generate target cavity specific sythetic compounds 
```
>>> import Synthesizer as sz
```
Set paths for the installed third-party softwares and default output folder to collect the intermediate result files and plots
```
>>> sz.Set_paths(LigBuilder_path='path to LigBuilderV3.0/',libiomp5_file='path to libiomp5.h file containing folder/',Output_dir='path to folder/')
```
To submit the query protein file of interest in PDB format
```
>>> sz.input_structure('path to pdbfile.pdb')
Cavity will output 16 cavity file(s)
```
Output shows total number of cavities predicted (in this case 16), which can be visualised by it's integer identifier
```
>>> cavity=4 #To select cavity number 4 as target cavity
>>> sz.cavity_view(cavity)
```

User can either directly choose a cavity number for the ligand synthesis
```
>>> cavity=4 #To select cavity number 4 as target cavity
>>> sz.compund_synthesis(cavity)
```
Or user can opt for cavity detection by submitting residue of interest in a TSV (Tab separated) file
```
>>> sz.cavity_detect(Residue_list.tsv)
```
Residue_list.tsv
```
$head -5 Residue_list.tsv
E       305
T       306
I       310
Y       316
V       466
```

##### Optional 

User can specify the number of compounds to synthesize (Default: 500)
```
>>> cavity=4 #To select cavity number 4 as target cavity
>>> lcount=800 #To synthesize 800 compounds 
>>> sz.compund_synthesis(cavity, ligand_count=lcount)
```
OR
```
>>> lcount=800 #To synthesize 800 compounds 
>>> sz.cavity_detect(Residue_list.tsv, ligand_count=lcount)
```

The output folder will contain the following files at the end of the successful execution of Synthesizer module
 1. Progress.st &emsp; #Status file containing Gcoupler progress
 2. Synth.csv &emsp; #CSV file containing SMILES of the synthetic compounds 
 3. PDBQT file &emsp; #Docking ready synthetic compounds



### Authenticator

To classify the synthetic compounds into binary classes based on molecular interaction
```
import Authenticator as au
```
To calculate interaction of individual synthetic compounds with the target cavity (which they are synthesized from)
```
>>> au.synthetic_interaction()
```
Additional arguments:
1. method: Statistical test to use for interaction cutoff estimation

| Parameter Name | Description |
| -------- | -------- |
| KS-test | Kolmogorov-Smirnov test (Default) |
| ES-test | Epps-Singleton test |
| AD-test | Anderson-Darling test |

3. p_val: Significance cutoff for the statistical test (Default: 0.05)
4. plot: Plot to visualize the distribution of HABs and LABs at each qualified cutoff

| Parameter Name | Description |
| -------- | -------- |
| Density | Density distribution plot (Default) |
| ECDF | Empirical cumulative distribution function plot |

To classify synthetic compounds into binary classes of HAB & LAB
```
>>> cutoff = -9  #user decided interaction cutoff for synthetic compound binary classification (Default: -7)
>>> au.synthetic_classify(cf=cutoff)
```

In case user want to opt for decoys as negative class against HABs
```
>>> au.synthetic_decoys()
```
Additional arguments:

| Arguments | Description |
| -------- | -------- |
| cf | User specified interaction cutoff for HABs, for which decoys will be generated(Default: -7) |
| decoy_csv | A csv file containg two coulmns. **SMILES** column containing the compund SMILES, and **Annotation** column containing it's annotation HAB(output from previous function) or Decoy |

The output folder will contain the following files at the end of the successful execution of Authenticator module
 1. Progress.st &emsp; #Status file containing Gcoupler progress
 2. Synth.csv &emsp; #CSV file containing SMILES of the synthetic compounds 
 3. PDBQT file &emsp; #Docking ready synthetic compounds
 4. Synth_BE.csv &emsp; #CSV file containing SMILES and interaction data of the synthetic compounds
 5. Labeled_cmps.csv &emsp; #CSV file containing SMILES and group data (HAB/LAB) of the synthetic compounds
 6. pdf file &emsp; #Distribution plots at each qualified balanced cutoff, containing information about the statistical test performed and respective p-value




### Generator

To classify the synthetic
