# -*- coding: utf-8 -*-
"""
@author: Raymond Moore
@institution: Mayo Clinic Foundation
"""
import argparse, os, sys, re, copy, glob  #, json
## Image by OpenCV
import cv2, base64
## Number crunching
import numpy as np
from scipy import stats
## Parse OME
import tifffile as tf
import xml.etree.ElementTree as ET
## Nice Graohing
import seaborn as sns
import pandas as pd
import matplotlib


from pprint import pprint

### GLOBAL CONSTANTS ###
omeExt = '*.ome.tiff'

########### OME.TIFF File Functions ###########
def findAllOMEFiles (baseDir, typ):
	mapFiles = {}
	nSamples = 0
	nFiles = 0
	if(typ):
		print("Assume Single OME.TIFF per Sample")
		
	else:
		print("Assume Multiple FOVs within Sample Dir")
		targets = list(filter(os.path.isdir, glob.glob( os.path.join( baseDir , "*") )))
		for dH in targets:
			mapFiles[os.path.basename(dH)] = list( glob.glob( os.path.join(dH,omeExt)) )
		nSamples = len( mapFiles.keys() )
		nFiles = sum([len(mapFiles[x]) for x in mapFiles if isinstance(mapFiles[x], list)])
	
	print("  Found "+str(nSamples)+" Samples and "+str(nFiles)+" OME.TIFF files")
	return mapFiles

def parseMarkerNamesFromOMETIFF(openOME):
	markers = []
	metadata = ET.fromstring(openOME.pages[0].description)
	for ele in metadata.findall('.//*'):
		chMeta = dict(ele.attrib)
		if 'Fluor' in chMeta:
			#pprint(chMeta)
			markers.append(chMeta['Name'])
	return markers
########### OME.TIFF Function ###########



########### Image Calculations Functions ###########
def getImageStats(opnImg):
	#opnImg = cv2.imread(imgFH)
	statDict = {
			'ImageMin':int(np.amin(opnImg)),
			'ImageMax':int(np.amax(opnImg)),
			'ImageMean':float(np.mean(opnImg)),
			'ImageMedian':float(np.median(opnImg)),
			'ImageStdDv':float(np.std(opnImg))
			}
	flatFull = opnImg.ravel()
	statDict["ImageP10"] = float(np.percentile(flatFull,10,0))
	statDict["ImageP90"] = float(np.percentile(flatFull,90,0))

	if statDict["ImageP10"] == 0:
		statDict["ImageP10"] = 0.1
		statDict["SNRp"] = float(statDict["ImageP90"] / 10)
	else:
		statDict["SNRp"] = float(statDict["ImageP90"] / statDict["ImageP10"])
	
	imgSub = flatFull[np.where(flatFull > 1)]
	statDict["NonBlankMean"] = float(np.mean(imgSub))
	statDict["NonBlankPercentage"] = float(len(imgSub)) / float(len(flatFull))
	return statDict
########### Image Calculations Functions ###########


def findAllImageMetrics(fileDict):
	premetrics = []
	for sample, omes in fileDict.items():
		print("  Looking into "+sample)
		totalOME = len(omes)-1
		for oI, omFh in enumerate(omes):
			fov = os.path.splitext( os.path.basename(omFh))[0].strip(".ome")
			## Open OME.TIFF
			openOME = tf.TiffFile(omFh);
			nDim = len(openOME.pages)
			if oI == 0:
				print("    "+fov+" with "+str(nDim)+" markers")
				print("    ...")
			if oI == totalOME:
				print("    "+fov+" with "+str(nDim)+" markers")
			makerList = parseMarkerNamesFromOMETIFF(openOME)
			## ASSERT METADATA INTEGRITY
			if nDim != len(makerList):
				sys.exit("ERROR: Metadata integritry issue!")
			for idx, mk in enumerate(makerList):
				tmpDF = getImageStats( openOME.pages[idx].asarray() )
				tmpDF['Sample'] = sample
				tmpDF['FOV'] = fov
				tmpDF['Marker'] = mk
				premetrics.append(tmpDF)
	metrics = pd.DataFrame(premetrics)
	return metrics



def generateWholeBatchBoxplots_FOVs(dataDict, typ):
	frst = list(dataDict.keys())[0]
	frstFov = list(dataDict[frst].keys())[0]
	makerList = list(dataDict[frst][frstFov].keys())
	mk = makerList[1]
	allMkVals = []
	for smp in dataDict.keys():
		allMkVals.append( list([ r[mk]['SNRp'] for k,r in dataDict[frst].items()]) )




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compile Percentile SNR metrics for assessment')
	parser.add_argument('-d', '--rootdir', help='Directory Containing OME.TIFFS', nargs='?', type=str, dest="inDir", metavar="DIR",required=True)
	parser.add_argument('--stitched', action='store_true', default=False)
	parser.add_argument('--save_data', action='store_true', default=False)
	args = parser.parse_args()
	
	allFileDict = findAllOMEFiles(args.inDir, args.stitched)
	allFileDict = findAllOMEFiles("B:\SHARED\MC1776 melanoma\MxIF\OME_TIFF_Images", False)
	
	allDataStatsDF = findAllImageMetrics(allFileDict)
	
	omFh = r'B:\SHARED\MC1776 melanoma\MxIF\OME_TIFF_Images\SU2C_9\SU2C_9_region_004.ome.tiff'
	# Add option to write Panda Dataframe out to csv
	if(args.save_data):
		allDataStatsDF.to_csv("snrp_data.csv")
		os.getcwd()
		
	if(typ):
		print("Assume Single OME.TIFF per Sample")
		
	else:
		print("Assume Multiple FOVs within Sample Dir")
		
	generateWholeBatchBoxplots_FOVs(allDataStats)
	dataDict = allDataStats
	
	