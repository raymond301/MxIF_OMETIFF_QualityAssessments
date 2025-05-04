import os,json,sys,re,csv,glob,datetime
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


if __name__ == '__main__':
	slidePath = sys.argv[1]	
	roiIdx = sys.argv[2]
	sNom = os.path.basename(slidePath)
	outdir=r'/research/bsi/projects/staff_analysis/m088378/MetericMaking/DATA'
	
	afiFh = list(filter(os.path.isfile, glob.glob( os.path.join(slidePath,sNom+"_Final","CellDIVE_"+sNom+"_"+roiIdx+".afi") ) ))[0]
	metadata = ET.parse(afiFh)
	roiIdx = os.path.basename(afiFh).replace(".afi","").split("_")[-1]

	allRoundsList = []
	for img in metadata.getroot().findall('Image'):
		runDesign = {}
		runDesign['Region'] = roiIdx
		ch, path = img.findall('./')
		runDesign['Marker'] = ch.text
		runDesign['Round'] = str(path.text).replace(sNom+"_",'').split('_')[0]
		fullPath = os.path.join(slidePath,sNom+"_Final",path.text)
		runDesign['Fullpath'] = fullPath
		runDesign['Fluor'] = str(path.text).replace(sNom+"_",'').split('_')[2]
		## Read in image
		chImage = tifffile.imread(fullPath)
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
		
	stepwiseTable = pd.DataFrame(allRoundsList)
	stepwiseTable['Platform'] = 'CellDive'
	stepwiseTable['Sample'] = sNom
	stepwiseTable.sort_values(['Sample','Round'], inplace=True)
	stepwiseTable.reset_index(drop=True, inplace=True)
	today = datetime.date.today().strftime('%Y%m%d')
	stepwiseTable.to_csv('{}/Panel_Metrics_{}_{}_{}.csv'.format(outdir,sNom,roiIdx,today), index=False)


