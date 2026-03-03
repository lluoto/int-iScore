# Integrated AI Predictions Interface Analysis Scores Toolkit (int-iScore)

<div align="center">

```
**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​*************************************************
*    ██╗ ███╗  ██╗ ████████╗ ██╗ ███████╗ ███████╗  ██████╗   ██████╗ ███████╗  *
*    ██║ ████╗ ██║ ╚══██╔══╝ ██║ ██╔════╝ ██╔════╝ ██╔═══██╗ ██╔══██╗ ██╔════╝  *
*    ██║ ██╔██╗██║    ██║    ██║ ███████╗ █████╗   ██║   ██║ ██████╔╝ █████╗    *
*    ██║ ██║╚████║    ██║    ██║ ╚════██║ ██╔══╝   ██║   ██║ ██╔══██╗ ██╔══╝    *
*    ██║ ██║ ╚███║    ██║    ██║ ███████║ ███████╗ ╚██████╔╝ ██║  ██║ ███████╗  *
*    ╚═╝ ╚═╝  ╚══╝    ╚═╝    ╚═╝ ╚══════╝ ╚══════╝  ╚═════╝  ╚═╝  ╚═╝ ╚══════╝  *
**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​*************************************************
```





**A comprehensive toolkit for multi-modal cellular imaging analysis**
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](docs/)

</div>

---

**Developed by:** Yulin Luo ([luoyl2022@shanghaitech.edu.cn](mailto:luoyl2022@shanghaitech.edu.cn))  
**Institution:** ShanghaiTech University



## Introduction

   int-iScore is an open-source platform for comprehensive analysis of the interface of AI-prediction protein complex models. Our toolkit streamlines the workflow from raw models data extraction to advanced quantitative analysis, enabling researchers to evaluate physic-chemial/statistic potential feature of protein-protein interactions(PPI).


---

#### Dependencies

* Numpy
* Modeller(https://salilab.org/modeller/registration.html)
* Pandas
* MPI4py
* Freesasa (https://github.com/freesasa/freesasa-python)
* Intercaat (https://github.com/eved1018/Intercaat)
* Biopython (https://github.com/biopython/biopython)
* PDB2PQR (https://github.com/Electrostatics/pdb2pqr)
* Alphafold3-score (https://github.com/Mingchenchen/AF3Score)
* CPscore (https://github.com/wallnerlab/ProQDock)

* DockQ (https://github.com/bjornwallner/DockQ) for integrated score evaluating.



#### Softwares

* UCSF Chimera (https://www.cgl.ucsf.edu/chimera/)
* CCP4 (http://www.ccp4.ac.uk/)
* APBS (https://github.com/Electrostatics/apbs)
* Frustratometer2 (http://frustratometer.qb.fcen.uba.ar/new_job)



---





### Overview

*Table 1: Physiochemistry properties/Empirical principles-based scoring functions selected for evaluation*

| Feature | Range | Description | Feature_Type | 
| ------- | ----- | ----------- | ------------ |
Ec | [-1,1] | Electrostatic balance at the interface | Interface      
Sc | [-1,1] | Geometric / Steric fit at the interface | Interface       
iPLDDT | [0, 1] | PLDDT score within interface 5 Å | Interface        
iPTM | [0, 1] | PTM score within interface 5 Å | Interface       
BSA | [0, 1] | Size of the interface|Interface      
Frustration | [0, 1] | Frustration of the interface| Interface      
CPscore | [0, 1]| Inter-residue contacts preference | Interface 
SOAP | [0, 1] | Statistic Potential score of the interface | Interface 
DOPE | [0, 1] | Statistic Potential score of the overall Structure | All_Atom       
int_score | [0, 1] | overall evaluation on complex models | hybrid 

   Extra Ramachandran check:  
   The Ramachandran plot of steric angles (φ,ψ) relationship has been widely employed for protein structure determination, validation, model building, and for a great variety of applications and analyses.\
   We check it as complementary criteria for distinguishing improper models from prediction decoys.

   
   DockQ Statistics on CAPRI data:\
        0    <  DockQ <  0.23 - Incorrect\
        0.23 <= DockQ <  0.49 - Acceptable quality\
        0.49 <= DockQ <  0.80 - Medium quality\
                DockQ >= 0.80 - High quality
---
#### Install from local source

**Requirements:** Python 3.10+

```bash


# Clone the repository
git clone https://github.com/lluoto/int-iScore.git
cd int-iScore

# Create conda environment with Python 3.10+
conda create -n int_iScore -f env.yml
conda activate int_iScore
# mpi4py should extra process
conda install mpi4py -y

# Unzip CPScore and frustrate calculation module 
tar -jxvf cp_svm.tar.bz2
tar -jxvf frustratometer.tar.bz2
```
# Modify intercaat_config.ini for utilize the qhull for better intercaat calculate performance



