import os
import json
import string
import csv
from Bio import PDB
import intercaat_functions as icaat
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from sklearn.cluster import DBSCAN,AffinityPropagation
from mpi4py import MPI
import subprocess
from modeller import *
from modeller import soap_pp
import math
import itertools
import numpy as np
import csv
import ast
import subprocess
import time
from tqdm import tqdm
from pathlib import Path
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import warnings
warnings.filterwarnings('ignore')

ss=np.load('3070_0082_representative_hi.npy')
print(ss)