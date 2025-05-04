import os,json,sys,re,csv,glob
import numpy as np
import pandas as pd
from pprint import pprint
import copy
import xml.etree.ElementTree as ET

import openslide
print("Openslide :",openslide.__version__)

import cv2
print("OpenCV :",cv2.__version__)

from scipy.stats import kurtosis
from scipy.stats import skew
from scipy.stats import iqr

from skimage import img_as_float, filters

import libtiff
libtiff.libtiff_ctypes.suppress_warnings()

base=r'/research/bsi/archive/PI/'

def getPIArchiveFolder(pi, di):
	return os.path.join(base,pi,'tertiary',di,'integrated','*BMSAim2*')

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
	#base=r'/research/bsi/archive/PI/Ilyas_Sumera_m052774/tertiary/300919.Human_Cholangiocarcinoma_MxIF/integrated'
	# Markovic_Svetomir_snm02/tertiary/s210155.CellSegmentation/integrated/


	potentialProjects = []
	potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Markovic_Svetomir_snm02', 's210155.CellSegmentation') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Ilyas_Sumera_m052774', '300919.Human_Cholangiocarcinoma_MxIF') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Yoon_Harry_m061620', 's210155.MxIF_MarkovicCollaboration') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Block_Matthew_msb04', 's210155.MxIF_MarkovicCollaboration') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Wang_Chen_m092469', 's302493.MxIF_Ovarian_Cancer') ) )))

	#pprint(potentialProjects)
	#sys.exit()

	for pProj in potentialProjects:
		dfList = []
		ppNom = os.path.basename(pProj)
		print(" Project: "+ppNom)

		targetSamples = list(filter(os.path.isdir, glob.glob( os.path.join(pProj,'*')) ))
		# Remove "SLIDE-" from CellDive
		#targetSamples = [x for x in targetSamples if not x.startswith("SLIDE")]
		print(" Samples: "+str(len(targetSamples)))

		for sample in targetSamples:
			sNom = os.path.basename(sample)
			if sNom.startswith("SLIDE"):
				afiFile = os.path.join(pProj,sNom,sNom+"_Final","CellDIVE_"+sNom+"_R0.afi") ### need to modify to catch other FOV regions
				if not os.path.isfile(afiFile):
					print("  > ERROR in AFI: "+afiFile)
					continue
				metadata = ET.parse(afiFile)

				allRoundsList = []
				for img in metadata.getroot().findall('Image'):
					runDesign = {}
					runDesign['Region'] = "R0"
					ch, path = img.findall('./')
					runDesign['Round'] = str(path.text).replace(sNom+"_",'').split('_')[0]
					fullPath = os.path.join(pProj,sNom,sNom+"_Final",path.text)
					runDesign['Fullpath'] = fullPath
					runDesign['Fluor'] = str(path.text).replace(sNom+"_",'').split('_')[2]
					#runDesign[roundName][cyto] = [ch.text, os.path.isfile(fullPath)]
					if os.path.isfile(fullPath):
						#print(fullPath)
						slide = openslide.OpenSlide(fullPath)
						runDesign['Width'] = slide.dimensions[0]
						runDesign['Height'] = slide.dimensions[1]
						runDesign['Pyramid'] = slide.level_count
						region = slide.read_region((0,0), 0, slide.dimensions)
						r, g, b, a = region.split()
						r_np = np.asarray(r)
						b_np = np.asarray(b)
						g_np = np.asarray(g)
						tmp1 = cv2.add(r_np, b_np, g_np)
						tmp2 = cv2.normalize(tmp1, None, 0, 255, cv2.NORM_MINMAX)
						runDesign['Min'] = int(np.amin(tmp2))
						runDesign['Max'] = int(np.amax(tmp2))
						runDesign['Mean'] = int(np.mean(tmp2))
						runDesign['StdDev'] = int(np.std(tmp2))

						flatFull = tmp2.flatten()
						runDesign["Percentile5"] = float(np.percentile(flatFull,5,0))
						runDesign["Percentile10"] = float(np.percentile(flatFull,10,0))
						runDesign["Percentile90"] = float(np.percentile(flatFull,90,0))
						runDesign["Percentile95"] = float(np.percentile(flatFull,95,0))
						runDesign["IQR"] = np.round(iqr(flatFull),3)
						## Describe Histogram
						runDesign["HistogramVariance"] = np.var(flatFull)
						runDesign["HistogramSkewness"] = skew(flatFull)
						runDesign["HistogramKurtosis"] = kurtosis(flatFull)
						## Calculate percentages
						imgSub = flatFull[np.where(flatFull > 0)]
						runDesign["NonBlankPercentage"] = float(len(imgSub)) / float(len(flatFull))
						
						### HANDLE SPECIAL EXCEPTION - WHEN IMAGE IS 100% BLANK!!!
						if runDesign['Max'] == 0:
							print("WARNING: IMAGE IS COMPLETLY BLANK!!")
							runDesign['OtsuThreshold'] = 0
							runDesign['OtsuSignalMean'] = 0.0
							runDesign['OtsuNoiseMean'] = 0.0
							runDesign["SNR"] = -1
						else:
							"""
							Signal to Noise Ratio based on otsu thresholding
							"""
							### Blur Detection: https://www.pyimagesearch.com/2015/09/07/blur-detection-with-opencv/
							threshold = filters.threshold_otsu(tmp2)
							runDesign['OtsuThreshold'] = int(threshold)
							runDesign['OtsuSignalMean'] = float( np.mean( flatFull[np.where(flatFull > threshold)] ) )
							runDesign['OtsuNoiseMean'] = float( np.mean( flatFull[np.where(flatFull < threshold)] ) )
							if runDesign["OtsuNoiseMean"] == 0:
								runDesign["SNR"] = float(runDesign["OtsuSignalMean"] / 0.0000001)
							else:
								runDesign["SNR"] = float(runDesign["OtsuSignalMean"] / runDesign["OtsuNoiseMean"])			
					
					allRoundsList.append(runDesign)

				stepwiseTable = pd.DataFrame(allRoundsList)
				stepwiseTable['Platform'] = 'CellDive'
				stepwiseTable['Sample'] = sNom
				stepwiseTable['Location'] = pProj
				dfList.append( copy.deepcopy(stepwiseTable) )
				#pprint(stepwiseTable)
				#sys.exit()

			else:
				panel = {}
				panelDirs = glob.glob(os.path.join(pProj,sample,"*_dapi_*"))
				#pprint(panelDirs)
				totReg = len(glob.glob( os.path.join(pProj,sample,'RegisteredImages','S002',"*mono_dapi_*.tif")))
				if totReg < 1:
					print("  "+os.path.basename(sample)+"  =>  "+str(totReg)+" FOVs")

				if len(panelDirs) < 1 :
					print("No Scans panel!")
					continue
				panel2 = getPanelDesignMap(panelDirs)
				
				allRoundsList = []
				targets = list(filter(os.path.isfile, glob.glob( os.path.join( pProj,sample,'AFRemoved', "*") ))) ## region_006
				tgNames = [os.path.basename(e).replace('-','_').upper() for e in targets]
				
				for stp, pnl in panel2.items():
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
								runDesign['Round'] = stp
								runDesign['Fullpath'] = fh
								runDesign['Fluor'] = getCytoByPosition(flrIdx)
								baseFh = os.path.basename(fh)
								runDesign['Region'] = "R"+baseFh.split('region_')[1].replace('.tif','')
							
								img16 = cv2.imread(fh, -1)
								runDesign['Width'] = img16.shape[1]
								runDesign['Height'] = img16.shape[0]
								runDesign['Pyramid'] = 2
								
								img8 = (img16/256).astype('uint8')
								runDesign['Min'] = int(np.amin(img8))
								runDesign['Max'] = int(np.amax(img8))
								runDesign['Mean'] = int(np.mean(img8))
								runDesign['StdDev'] = int(np.std(img8))
								
								flatFull = img8.flatten()
								runDesign["Percentile5"] = float(np.percentile(flatFull,5,0))
								runDesign["Percentile10"] = float(np.percentile(flatFull,10,0))
								runDesign["Percentile90"] = float(np.percentile(flatFull,90,0))
								runDesign["Percentile95"] = float(np.percentile(flatFull,95,0))
								runDesign["IQR"] = np.round(iqr(flatFull),3)
								## Describe Histogram
								runDesign["HistogramVariance"] = np.var(flatFull)
								runDesign["HistogramSkewness"] = skew(flatFull)
								runDesign["HistogramKurtosis"] = kurtosis(flatFull)
								## Calculate percentages
								imgSub = flatFull[np.where(flatFull > 0)]
								runDesign["NonBlankPercentage"] = float(len(imgSub)) / float(len(flatFull))
								
								### HANDLE SPECIAL EXCEPTION - WHEN IMAGE IS 100% BLANK!!!
								if runDesign['Max'] == 0:
									print("WARNING: IMAGE IS COMPLETLY BLANK!! "+fh)
									runDesign['OtsuThreshold'] = 0
									runDesign['OtsuSignalMean'] = 0.0
									runDesign['OtsuNoiseMean'] = 0.0
									runDesign["SNR"] = -1
								else:
									"""
									Signal to Noise Ratio based on otsu thresholding
									"""
									### Blur Detection: https://www.pyimagesearch.com/2015/09/07/blur-detection-with-opencv/
									threshold = filters.threshold_otsu(img8)
									runDesign['OtsuThreshold'] = int(threshold)
									runDesign['OtsuSignalMean'] = float( np.mean( flatFull[np.where(flatFull > threshold)] ) )
									runDesign['OtsuNoiseMean'] = float( np.mean( flatFull[np.where(flatFull < threshold)] ) )
									if runDesign["OtsuNoiseMean"] == 0:
										runDesign["SNR"] = float(runDesign["OtsuSignalMean"] / 0.0000001)
									else:
										runDesign["SNR"] = float(runDesign["OtsuSignalMean"] / runDesign["OtsuNoiseMean"])	
							
								allRoundsList.append(runDesign)

				stepwiseTable = pd.DataFrame(allRoundsList)
				stepwiseTable['Platform'] = 'InCell'
				stepwiseTable['Sample'] = sNom
				stepwiseTable['Location'] = pProj
				dfList.append( copy.deepcopy(stepwiseTable) )
				#pprint(stepwiseTable)
				#sys.exit()
		# sys.exit()

		#print(dfList)
		df = pd.concat(dfList, axis=0) # join='inner'
		df.sort_values(['Sample','Round'], inplace=True)
		df.reset_index(drop=True, inplace=True)
		df.to_csv('Panel_Metrics_{}.csv'.format(ppNom), index=False)
	



