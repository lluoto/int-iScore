#!/usr/bin/env python3
import sys, os
sys.path.insert(0, "/media/cuixi/data0/lluoto/int-iScore/src")
os.chdir("/media/cuixi/data0/lluoto/int-iScore/src/int_iscore")
print("CWD:", os.getcwd())
print("Config exists:", os.path.exists("intercaat_config.ini"))
from configparser import ConfigParser
config = ConfigParser()
config.read("intercaat_config.ini")
print("Sections:", config.sections())
print("qvoronoi_path:", config.get("qvoronoi_path", "run_python_version"))
