# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:34:52 2021

@author: m088378
"""
import argparse, os, sys, glob
from datetime import datetime
import tempfile
## Image by OpenCV
import cv2, base64
## Number crunching
import numpy as np
## Nice Graohing
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

quant_exe = "tsv"

def LoadAllQuantFiles (baseDir):
	targets = list(glob.glob( os.path.join( baseDir , "*"+quant_exe) ))
	tmpDFs = []
	for dH in targets:
		print("  Loading: "+os.path.basename(dH))
		tmpDFs.append( pd.read_csv(dH) )
	allFhs = pd.concat(tmpDFs, ignore_index=True)
	allFhs = allFhs[allFhs.columns.drop(list(allFhs.filter(regex='Cytoplasm:')))]







if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compile Single Cell metrics for assessment')
	parser.add_argument('-d', '--rootdir', help='Directory Containing QuPath Quant Files', nargs='?', type=str, dest="inDir", metavar="DIR",required=True)
	parser.add_argument('--multicomponents', action='store_true', default=False)
	parser.add_argument('--save_data', action='store_true', default=False)
	args = parser.parse_args()

	allDataDF LoadAllQuantFiles(args.inDir)
	baseDir = r'Y:\Studies\AnalysisPipelineProjects\GM_CSF\Quantifications'
