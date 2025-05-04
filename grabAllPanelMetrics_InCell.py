import os,json,sys,re,csv,glob,datetime
import gc
import numpy as np
import pandas as pd
from pprint import pprint
import copy
import xml.etree.ElementTree as ET

#import openslide
#print("Openslide :",openslide.__version__)

import cv2
print("OpenCV :",cv2.__version__)

from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import iqr

from skimage import img_as_float, filters

#import libtiff
#libtiff.libtiff_ctypes.suppress_warnings()
import tifffile
print("TiffFile :",tifffile.__version__)

import mahotas as mt
# from mahotas.features import surf

def markerNameUpdate(nom):
	if "HLA" in nom:
		nom = nom.replace('HLA_','HLA-')
	if "BMS" in nom:
		nom = nom.replace('_BMS','-BMS')
	if "Gores" in nom:
		nom = nom.replace('_Gores','-Gores')
	if "Caspase" in nom:
		nom = nom.replace('Caspase_','Caspase-')
	if "isotype" in nom:
		nom = nom.replace('isotype_','isotype-')
	if "T_bet" in nom:
		nom = nom.replace('T_bet','T-bet')
	if "PD_L" in nom:
		nom = nom.replace('PD_L','PD-L')
	if "Gran_B" in nom:
		nom = nom.replace('Gran_B','Gran-B')
	if "GATA_3" in nom:
		nom = nom.replace('GATA_3','GATA-3')
	return nom

def getPanelDesignMap(panelDirs):
	for dd in panelDirs:
		nom = markerNameUpdate( os.path.basename(dd) )
		dbase = nom.split("_")
		if dbase[4] == "fitc" and len(dbase) > 7:
			panel[dbase[0]] = [dbase[3],dbase[5],dbase[7]]
		elif dbase[4] == "cy3" and len(dbase) > 5:
			panel[dbase[0]] = ['-',dbase[3],dbase[5]]
		elif dbase[4] == "cy3":
			panel[dbase[0]] = ['-',dbase[3],'-']
		elif dbase[4] == "cy5":
			panel[dbase[0]] = ['-','-',dbase[3]]
		elif dbase[4] == "fitc" and dbase[6] == "cy5":
			panel[dbase[0]] = ['-','-',dbase[3]]
		else :
			pprint(dbase)
	ptbl=dict([(k, [x.upper() for x in v]) for (k,v) in panel.items()])
	return ptbl

def getCytoByPosition(p):
	if p == 0:
		return 'Cy2'
	elif p == 1:
		return 'Cy3'
	elif p == 2:
		return 'Cy5'
	else:
		return '???'



if __name__ == '__main__':
	slidePath = sys.argv[1]	
	roiIdx = sys.argv[2]
	sNom = os.path.basename(slidePath)
	outdir=r'/research/bsi/projects/staff_analysis/m088378/MetericMaking/DATA'
	
	panel = {}
	panelDirs = glob.glob(os.path.join(slidePath,"*_dapi_*"))
	if len(panelDirs) < 1 :
		print("No Scans panel!")
		sys.exit(1)
		
	panel2 = getPanelDesignMap(panelDirs)
	#pprint(panel2)
	
	allRoundsList = []
	targets = list(filter(os.path.isfile, glob.glob( os.path.join( slidePath,'AFRemoved', "*region_{}*".format(roiIdx)) ))) ## region_006
	tgNames = [os.path.basename(e).replace('-','_').upper() for e in targets]
	
	
	allRoundsList = []
	for stp, pnl in panel2.items():
		print(stp)
		for flrIdx, ele in enumerate(pnl):
			if ele == "-":
				continue
			elif ele == "BLEACH":
				continue
			elif ele == "BKGND":
				continue
			else:
				idxsF = [tgNames.index(l) for l in tgNames if l.startswith(ele.replace('-','_'))]
				subListFiles = [targets[index] for index in idxsF]
				for fh in subListFiles:
					runDesign = {}
					runDesign['Region'] = roiIdx
					runDesign['Round'] = stp
					runDesign['Marker'] = ele
					runDesign['Fullpath'] = fh
					runDesign['Fluor'] = getCytoByPosition(flrIdx)
					## Read in image
					chImage = tifffile.imread(fh)
					runDesign['Width'] = chImage.shape[0]
					runDesign['Height'] = chImage.shape[1]
					
					runDesign['Min'] = int(np.amin(chImage))
					runDesign['Max'] = int(np.amax(chImage))
					runDesign['Mean'] = int(np.mean(chImage))
					runDesign['StdDev'] = int(np.std(chImage))
					
					flatFull = chImage.flatten()
					runDesign["Percentile5"] = float(np.percentile(flatFull,5,0))
					runDesign["Percentile10"] = float(np.percentile(flatFull,10,0))
					runDesign["Percentile90"] = float(np.percentile(flatFull,90,0))
					runDesign["Percentile95"] = float(np.percentile(flatFull,95,0))
					runDesign["IQR"] = np.round(iqr(flatFull),3)
					runDesign["Variance"] = np.round(np.var(flatFull),-1)
					## Describe Histogram
					runDesign["HistogramVariance"] = np.var(flatFull)
					runDesign["HistogramSkewness"] = skew(flatFull)
					runDesign["HistogramKurtosis"] = kurtosis(flatFull)
					## Calculate percentages
					imgSub = flatFull[np.where(flatFull > 0)]
					runDesign["NonBlankPercentage"] = float(len(imgSub)) / float(len(flatFull))
					imgSub2 = flatFull[np.where(flatFull >= 65534)]
					runDesign["SaturationPercentage"] = float(len(imgSub2)) / float(len(flatFull))
					## SNR based on Percentile
					if runDesign["Percentile10"] == 0:
						runDesign["SNRp"] = float(runDesign["Percentile90"] / 0.000001)
					else:
						runDesign["SNRp"] = float(runDesign["Percentile90"] / runDesign["Percentile10"])
					## SNR based on Otsu
					runDesign['OtsuThreshold'] = 0
					if runDesign['Max'] != 0:
						threshold = filters.threshold_otsu(chImage)
						runDesign['OtsuThreshold'] = int(threshold)
						runDesign['OtsuSignalMean'] = float( np.mean( flatFull[np.where(flatFull > threshold)] ) )
						runDesign['OtsuNoiseMean'] = float( np.mean( flatFull[np.where(flatFull < threshold)] ) )
						if runDesign["OtsuNoiseMean"] == 0:
							runDesign["SNRo"] = float(runDesign["OtsuSignalMean"] / 0.0000001)
						else:
							runDesign["SNRo"] = float(runDesign["OtsuSignalMean"] / runDesign["OtsuNoiseMean"])
					## SNR based on Z Score
					if runDesign["Variance"] == 0:
						runDesign["SNRz"] = float(runDesign["Mean"] / 0.000001)
					else:
						runDesign["SNRz"] = float(runDesign["Mean"] / runDesign["Variance"])
					"""
					An implementation of the Brenner autofocus metric
					Brenner, J. F. et al (1976). An automated microscope for cytologic research
					a preliminary evaluation. Journal of Histochemistry & Cytochemistry, 24(1),
					100–111. http://doi.org/10.1177/24.1.1254907
					"""
					rows = runDesign['Height']
					columns = runDesign['Width'] - 2
					repImg = np.zeros((rows, columns))
					try:					
						repImg[:] = ((chImage[:, 0:-2] - chImage[:, 2:]) ** 2)
						runDesign["BrennerScore"] = repImg.sum()
					except:
						runDesign["BrennerScore"] = -1

					#def variance_of_laplacian(image):
					"""
					compute the Laplacian of the image and then return the focus
					measure, which is simply the variance of the Laplacian
					"""
					runDesign["LaplacianVariance"] = cv2.Laplacian(chImage, cv2.CV_64F).var()
					"""
					These are texture features, based on the adjacency matrix (the adjacency matrix stores in position (i,j) the number of 
					times that a pixel takes the value i next to a pixel with the value j. Given different ways to define next to, you obtain 
					slightly different variations of the features. Standard practice is to average them out across the directions to get some rotational invariance.
					"""
					# calculate haralick texture features for 4 types of adjacency
					textures = mt.features.haralick(chImage)
					# Only the first 13 features are implemented. The last (14th) feature is normally considered to be unstable
					ht_means = textures.mean(axis=0)
					for idx, htD in enumerate(ht_means):
						runDesign["Haralick"+str(idx)] = htD
					"""
					This is a companion to the paper Determining the subcellular location of new proteins from microscope 
					images using local features by Coelho et al. (2013).
					"""
					# https://mahotas.readthedocs.io/en/latest/surf.html
					#spoints = surf.surf(chImage)
					#runDesign["SurfPoints"] = len(spoints)
					#print("Nr points: {}".format())
					
					#Calculate a threshold according to the Riddler-Calvard method.
					runDesign["Threshold_RiddlerCalvard"] = mt.rc(chImage, ignore_zeros=True)
					allRoundsList.append(runDesign)
					
					# Force a garbage collection
					gc.collect()

	stepwiseTable = pd.DataFrame(allRoundsList)
	stepwiseTable['Platform'] = 'INCell'
	stepwiseTable['Sample'] = sNom
	stepwiseTable.sort_values(['Sample','Round'], inplace=True)
	stepwiseTable.reset_index(drop=True, inplace=True)
	# pprint(stepwiseTable)
	today = datetime.date.today().strftime('%Y%m%d')
	stepwiseTable.to_csv('{}/Panel_Metrics_{}_{}_{}.csv'.format(outdir,sNom,roiIdx,today), index=False)


