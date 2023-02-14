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
```
Output shows total number of cavities predicted which can be visualised by it's integer identifier
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

