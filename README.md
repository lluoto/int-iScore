# Integrated AI Predictions Interface Analysis Scores Toolkit (int-iScore)

<div align="center">

```
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
*
*
*
*    ██╗  █████╗    ██╗ ██████╗       █████╗     *
*    ██║ ██╔══██╗ ██╔══██╗   ██║ ██╔══██╗     ██╔══██╗
*    ██║ ██████╔╝ ███████║   ██║ ██████╔╝     ███████║
*    ██║ ██╔═══╝  ██╔══██║   ██║ ██╔═══╝      ██╔══██║
*    ██║ ██║      ██║  ██║   ██║ ██║          ██║  ██║
*    ╚═╝ ╚═╝      ╚═╝  ╚═╝   ╚═╝ ╚═╝          ╚═╝  ╚═╝ 
*           
*        
*        
*       
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
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

iAI-PIA is an open-source platform for comprehensive analysis of the interface of AI-prediction protein complex models. Our toolkit streamlines the workflow from raw models data extraction to advanced quantitative analysis, enabling researchers to evaluate physic-chemial/statistic potential feature of protein-protein interaction(PPI).


---

#### Dependencies

* numpy
* Modeller(https://salilab.org/modeller/registration.html)
* pandas
* mpi4py
* freesasa (https://github.com/freesasa/freesasa-python)
* intercaat (https://github.com/eved1018/Intercaat)
* biopython (https://github.com/biopython/biopython)
* UCSF Chimera (https://www.cgl.ucsf.edu/chimera/)
* CCP4 (http://www.ccp4.ac.uk/)
* PDB2PQR (https://github.com/Electrostatics/pdb2pqr)
* APBS (https://github.com/Electrostatics/apbs)
* frustratometer2 (http://frustratometer.qb.fcen.uba.ar/new_job)
* Alphafold3-score (https://github.com/Mingchenchen/AF3Score)

* DockQ (https://github.com/bjornwallner/DockQ) for integrated score evaluating.



---

---

#### Softwares

* UCSF Chimera (https://www.cgl.ucsf.edu/chimera/)
* CCP4 (http://www.ccp4.ac.uk/)
* APBS (https://github.com/Electrostatics/apbs)
* frustratometer2 (http://frustratometer.qb.fcen.uba.ar/new_job)



---





### Overview

---

*Table 1: Physiochemistry properties/Empirical principles-based scoring functions selected for evaluation*

=================================================================================================
Feature   Range  Direction  Description                                            Feature_Type 
-------------------------------------------------------------------------------------------------
EC:       [-1,1] Positive   Electrostatic balance at the interface                 Interface      
Sc:       [-1,1] Positive   Geometric / Steric fit at the interface                Interface       
iPLDDT:   [0, 1] Positive   PLDDT score within interface 5 Å                       Interface        
iPTM:     [0, 1] Positive   PTM score within interface 5 Å                         Interface       
BSA:      [0, 1] Ambiguous  Size of the interface                                  Interface      
Frustra:  [0, 1] Positive   Frustration of the interface                           Interface      
CPscore:  [0, 1] Positive   inter-residue contacts preference                      Interface 
SOAP:     [0, 1] Positive   statistic potential score of the interface             Interface 
DOPE:     [0, 1] Positive   statistic potential score of the overall Structure      All_Atom       
=================================================================================================
   Extra Ramachandran check:  
   The Ramachandran plot of torsion angles (φ,ψ) has been widely employed for protein structure determination, validation, model building, and for a great variety of applications and analyses, since it was proposed in 1963 by late G.N. Ramachandran and coworkers ( Ramachandran et al., 1963; Ramakrishnan and Ramachandran, 1965; Ramachandran and Sasisekharan, 1968). We check it as complementary criteria for distinguishing improper models from prediction decoys.
    
	
   DockQ Statistics on CAPRI data:  
    0    <  DockQ <  0.23 - Incorrect
    0.23 <= DockQ <  0.49 - Acceptable quality
    0.49 <= DockQ <  0.80 - Medium quality
            DockQ >= 0.80 - High quality
---
#### Install from local source

**Requirements:** Python 3.7-3.8

```bash
# Create conda environment with Python 3.7 or 3.8
conda create -n interface_score python=3.7
conda activate interface_score

# Clone the repository
git clone https://github.com/SunLab-SH/iPA.git
cd iPA

# Install dependencies and the package
pip install -r requirements.txt
```




