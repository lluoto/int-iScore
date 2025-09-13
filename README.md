# Integrated Processing and Analysis Toolkit (iPA)

<div align="center">

```
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
*    ██╗ ██████╗   █████╗     *
*    ██║ ██╔══██╗ ██╔══██╗    *
*    ██║ ██████╔╝ ███████║    *
*    ██║ ██╔═══╝  ██╔══██║    *
*    ██║ ██║      ██║  ██║    *
*    ╚═╝ ╚═╝      ╚═╝  ╚═╝    *
​**​*​**​*​**​*​**​*​**​*​**​*​**​*​**​**​**​**​
```





**A comprehensive toolkit for multi-modal cellular imaging analysis**
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](docs/)

</div>

---

**Developed by:** Angdi Li ([liad@shanghaitech.edu.cn](mailto:liad@shanghaitech.edu.cn))  
**Institution:** ShanghaiTech University



## Introduction

iPA is an open-source Python platform for comprehensive analysis of cellular and subcellular imaging data. Our toolkit streamlines the workflow from raw image processing to advanced quantitative analysis, enabling researchers to extract biologically meaningful features from diverse microscopy modalities.



For detailed documentation, tutorials, and examples, please visit:  
[Documentation on ReadTheDocs](https://ipa.readthedocs.io/en/latest/)

---

#### Dependencies

* numpy
* scipy
* pandas
* matplotlib
* scikit-image
* tifffile
* Pillow
* opencv-python
* mrcfile
* n2v
* plotly
* seaborn
* torch
* torchvision
* tensorflow
* tqdm
* h5py
* scikit-learn


---

#### Example data
* Example data are available from: [google doc link](https://drive.google.com/drive/folders/12bhaITv_xdNvs-pBwr6SSJviRlOkQ9cW?usp=drive_link)   
* Download the `data.zip` file and extract it to the project root directory
* After extraction, you should have the following structure:
```
iPA/
├── data/          # <- Extracted data folder should be here
├── ipa/
├── examples/
├── README.md
└── ...
```




### Overview


![iPA Analysis Pipeline](./workflow_images/figure_1_v25.jpg)

*Figure 1: iPA Analysis Pipeline Overview*




#### Install from local source

**Requirements:** Python 3.7-3.8

```bash
# Create conda environment with Python 3.7 or 3.8
conda create -n ipa_env python=3.7
conda activate ipa_env

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
