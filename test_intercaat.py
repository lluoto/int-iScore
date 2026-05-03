#!/usr/bin/env python3
import os, sys
sys.path.insert(0, "/media/cuixi/data0/lluoto/int-iScore/src")
os.chdir("/media/cuixi/data0/lluoto/int-iScore/src/int_iscore")
print("CWD:", os.getcwd())
print("Config exists:", os.path.exists("intercaat_config.ini"))

# Test reading config directly
from configparser import ConfigParser
c = ConfigParser()
c.read("intercaat_config.ini")
print("Sections:", c.sections())

# Now import intercaat_functions and test
import intercaat_functions as icaat
print("intercaat_functions file:", icaat.__file__)

# Test the run_voro function
points = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]
try:
    result = icaat.run_voro(points)
    print("run_voro result:", result)
except Exception as e:
    print("run_voro error:", e)
    import traceback
    traceback.print_exc()
