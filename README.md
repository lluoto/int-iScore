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
* json
* Modeller
* ast
* pandas
* csv
* mpi4py
* math
* freesasa
* subprocess
* DockQ
* intercaat
* biopython
* chimera
* CCP4(http://www.ccp4.ac.uk/)
* PDB2PQR
* APBS
* frustratometer2
* Alphafold3-score



---





### Overview



*Table 1: iPA Analysis Pipeline Overview*

=================================================================================================================================
Feature   Range  Direction  Description                                            Feature_Type   Stability   Native_values (DB3)
---------------------------------------------------------------------------------------------------------------------------------
EC:       [-1,1] Positive   Electrostatic balance at the interface                 Interface      
Sc:       [-1,1] Positive   Geometric / Steric fit at the interface                Interface       
iPLDDT:   [0, 1] Positive   PLDDT score within interface 5 Å                       Interface        
iPTM:     [0, 1] Positive   PTM score within interface 5 Å                         Interface       
BSA:      [0, 1] Ambiguous  Size of the interface                                  Interface      
Frustra:  [0, 1] Positive   Frustration of the interface                           Interface      
CPscore:  [0, 1] Positive   inter-residue contacts preference                      Interface 
SOAP:     [0, 1] Positive   statistic potential score of the interface             Interface 
DOPE:     [0, 1] Positive   statistic potential score of the overall Structure      All_Atom       
=================================================================================================================================
   Extra Ramachandran check:  
    
	
   DockQ Statistics on CAPRI data:  
    0    <  DockQ <  0.23 - Incorrect
    0.23 <= DockQ <  0.49 - Acceptable quality
    0.49 <= DockQ <  0.80 - Medium quality
            DockQ >= 0.80 - High quality
#### Install from local source

**Requirements:** Python 3.7-3.8

```bash
# Create conda environment with Python 3.7 or 3.8
conda create -n interface_score python=3.7
conda activate interface_score

# Install Git LFS (if not already installed)
# On Ubuntu/Debian:
sudo apt install git-lfs
# On macOS:
brew install git-lfs
# On Windows: download from https://git-lfs.github.io/

# Initialize Git LFS
git lfs install

# Clone the repository
git clone https://github.com/SunLab-SH/iPA.git
cd iPA

# Pull large files (including .pth model files)
git lfs pull

# Install dependencies and the package
pip install -r requirements.txt
pip install -e .
```



<!-- ## Citation

If you use iPA in your research, please cite:

```bibtex
@software{iPA2024,
  title={Integrated Processing and Analysis Toolkit for Multi-Scale Cellular Imaging},
  author={Li, Angdi and others},
  year={2024},
  url={https://github.com/SunLab-SH/iPA}
}
``` -->
