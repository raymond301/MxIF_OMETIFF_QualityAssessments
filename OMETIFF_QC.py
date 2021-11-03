# -*- coding: utf-8 -*-
"""
@author: Raymond Moore
@institution: Mayo Clinic Foundation
"""
import argparse, os, sys, glob
from datetime import datetime
import tempfile
## Image by OpenCV
import cv2, base64
## Number crunching
import numpy as np
#from scipy import stats
from skimage import data, img_as_float, io, filters
## Parse OME
import tifffile as tf
import xml.etree.ElementTree as ET
## Nice Graohing
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


from pprint import pprint
import warnings
warnings.filterwarnings("ignore")

### GLOBAL CONSTANTS ###
omeExt = '*.ome.tiff'
verbose = False

########### OME.TIFF File Functions ###########
def findAllOMEFiles (baseDir, typ):
	mapFiles = {}
	nSamples = 0
	nFiles = 0
	if(typ):
		print("Assume Single OME.TIFF per Sample")
		## Future development for when no longer rely on FOV selection
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
#opnImg = openOME.pages[idx].asarray()
def getImageStats(opnImg):
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

	## Calculate percentages
	imgSub = flatFull[np.where(flatFull > 1)]
	statDict["NonBlankMean"] = float(np.mean(imgSub))
	statDict["NonBlankPercentage"] = float(len(imgSub)) / float(len(flatFull))

	### HANDLE SPECIAL EXCEPTION - WHEN IMAGE IS 100% BLANK!!!
	if statDict['ImageMax'] == 0:
		if verbose:
			print("WARNING: IMAGE IS COMPLETLY BLANK!!")
		statDict['OtsuThreshold'] = 0
		statDict['OtsuSignalMean'] = 0.0
		statDict['OtsuNoiseMean'] = 0.0
		statDict["SNRo"] = 0.0
		statDict["BrennerScore"] = 0.0
	else:
		"""
		Signal to Noise Ratio based on otsu thresholding
		"""
		### Blur Detection: https://www.pyimagesearch.com/2015/09/07/blur-detection-with-opencv/
		threshold = filters.threshold_otsu(opnImg)
		statDict['OtsuThreshold'] = int(threshold)
		#pprint("Len of AF Img Flat:"+str(len(afImgFlat)))
		statDict['OtsuSignalMean'] = float( np.mean( flatFull[np.where(flatFull > threshold)] ) )
		statDict['OtsuNoiseMean'] = float( np.mean( flatFull[np.where(flatFull < threshold)] ) )

		if statDict["OtsuNoiseMean"] == 0:
			statDict["SNRo"] = float(statDict["OtsuSignalMean"] / 10)
		else:
			statDict["SNRo"] = float(statDict["OtsuSignalMean"] / statDict["OtsuNoiseMean"])

		"""
	    An implementation of the Brenner autofocus metric
	    Brenner, J. F. et al (1976). An automated microscope for cytologic research
	    a preliminary evaluation. Journal of Histochemistry & Cytochemistry, 24(1),
	    100–111. http://doi.org/10.1177/24.1.1254907
	    """
		rows = opnImg.shape[0]
		columns = opnImg.shape[1] - 2
		repImg = np.zeros((rows, columns))
		repImg[:] = ((opnImg[:, 0:-2] - opnImg[:, 2:]) ** 2)
		statDict["BrennerScore"] = repImg.sum()

	return statDict
########### Image Calculations Functions ###########

#fileDict = t2
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
			if verbose:
				print("    "+fov+" with "+str(nDim)+" markers")
			else:
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



def generateWholeBatchBoxplots_FOVs(pDf):
	nMarks = pDf['Marker'].nunique()
	grouped = pDf.loc[:,['Marker', 'SNRp']].groupby(['Marker']).median().sort_values(by='SNRp', ascending=False)
	#sns.set_theme(style="ticks")
	sns.set(style="ticks")
	# Initialize the figure with a logarithmic x axis
	f, ax = plt.subplots(figsize=(25, nMarks))
	ax.set_xscale("log")
	# Plot the orbital period with horizontal boxes
	sns.boxplot(x="SNRp", y="Marker", data=pDf, order=grouped.index, whis=[0, 100], width=.5, palette="vlag")
	# Add in points to show each observation
	sns.stripplot(x="SNRp", y="Marker", data=pDf, order=grouped.index, size=2, color=".3", linewidth=0)
	# Tweak the visual presentation
	ax.xaxis.grid(True)
	ax.set(ylabel="")
	sns.despine(trim=True, left=True)

	with tempfile.TemporaryFile(suffix=".png") as tmpfile:
		plt.savefig(tmpfile, format="png") # File position is at the end of the file.
		tmpfile.seek(0) # Rewind the file. (0: the beginning of the file)
		mainPlot64 = base64.b64encode(tmpfile.read()).decode('utf-8')

	template = """
	</br><div style="width:99%;">
	<h3>SNR Distributions Across Dataset</h3>
	<img style="max-width: 100%;" src="data:image/png;base64, {}" alt="SNRp Boxplot Graph" />
	</div>
	"""
	return template.format(mainPlot64)


def generateWholeBatchSamples_FOVs(pDf):
	nSamples = pDf['Sample'].nunique()
	sns.set(style="ticks", palette="pastel")
	#sns.set_theme(style="ticks", palette="pastel")
	f, ax = plt.subplots(figsize=(int(nSamples/1.2)+3, 6))
	ax.set_yscale("log")
	# Draw a nested boxplot to show bills by day and time
	sns.boxplot(x="Sample", y="SNRp",width=4, hue="Sample", data=pDf)
	ax.get_legend().remove()
	sns.despine(offset=10, trim=True)
	ax.tick_params(axis='x', rotation=20)

	with tempfile.TemporaryFile(suffix=".png") as tmpfile:
		plt.savefig(tmpfile, format="png") # File position is at the end of the file.
		tmpfile.seek(0) # Rewind the file. (0: the beginning of the file)
		mainPlot64 = base64.b64encode(tmpfile.read()).decode('utf-8')

	template = """
	</br><div  style="width:90%;">
	<h3>SNR Distributions Across Samples</h3>
	<img style="max-width: 100%;" src="data:image/png;base64, {}" alt="By Sample Boxplot Graph" />
	</div>
	"""
	return template.format(mainPlot64)



html_str_css = """
<style>
table.minimalistBlack {
		border: 3px solid #000000;
		text-align: center;
		border-collapse: collapse;
}
table.minimalistBlack td, table.minimalistBlack th {
		border: 1px solid #000000;
		padding: 5px 8px;
}
table.minimalistBlack tbody td {
		font-size: 15px;
}
table.minimalistBlack thead {
		background: #CFCFCF;
		background: -moz-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%);
		background: -webkit-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%);
		background: linear-gradient(to bottom, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%);
		border-bottom: 3px solid #000000;
}
table.minimalistBlack thead th {
		font-size: 19px;
		font-weight: bold;
		color: #000000;
		text-align: left;
}
table.minimalistBlack tfoot td {
		font-size: 14px;
}
</style>
"""

def getSummaryHTMLTable(pDf):
	template = """
	<table class="minimalistBlack">
	<tr>
		<th>Number of Samples</th><td>{}</td>
	</tr> <tr>
		<th>Number of Images</th><td>{}</td>
	</tr> <tr>
		<th>Number of Markers</th><td>{}</td>
	</tr> <tr>
		<th>Timestamp</th><td>{}</td>
	</tr>
	</table>
	"""
	nSamples = pDf['Sample'].nunique()
	nMarks = pDf['Marker'].nunique()
	nRow = len(pDf.index)
	tm = str(datetime.now().strftime('%Y-%m-%d T%H:%M:%S'))
	return template.format(nSamples, nRow, nMarks, tm)


def getCutoffImages(dfSub):
	t = dfSub.loc[(dfSub.Quality == 'Poor'),]
	t['DropOuts'] = t['Sample']+":"+t['FOV']
	return ",  ".join(t['DropOuts'])

def getQualityPieChart(dfSub):
	#create pie chart
	plt.clf()
	pData = dfSub.groupby("Quality")["Quality"].count()
	pData.plot.pie(autopct="%.1f%%", figsize=(5,4))
	fig = plt.gcf()
	with tempfile.TemporaryFile(suffix=".png") as tmpfile:
		fig.savefig(tmpfile, format="png")
		#plt.savefig(tmpfile, format="png") # File position is at the end of the file.
		tmpfile.seek(0) # Rewind the file. (0: the beginning of the file)
		mainPlot64 = base64.b64encode(tmpfile.read()).decode('utf-8')
	template = """
	<img style="max-width: 250px;" src="data:image/png;base64, {}" alt="Quality PieChart" />
	"""
	plt.close()
	return template.format(mainPlot64)

def getQualityScatterPlot(dfSub, metric):
	f, ax = plt.subplots(figsize=(14, 5))
	ax.set_xscale("log")
	graph = sns.scatterplot(y="NonBlankPercentage", x=metric, palette="deep", hue="Quality", data=dfSub, ax=ax)
	graph.axvline(0.916, color='r') ### 2.5 ~ 0.916 Log
	graph.axhline(0.03, color='r')
	plt.legend(bbox_to_anchor=(0.93, 1), loc=2, borderaxespad=0.)
	with tempfile.TemporaryFile(suffix=".png") as tmpfile:
		plt.savefig(tmpfile, format="png") # File position is at the end of the file.
		tmpfile.seek(0) # Rewind the file. (0: the beginning of the file)
		mainPlot64 = base64.b64encode(tmpfile.read()).decode('utf-8')
	template = """
	<img style="max-width:700px;" src="data:image/png;base64, {}" alt="Quality Scatterplot" />
	"""
	plt.close()
	return template.format(mainPlot64)


def generateByMarkerSNRTable(pDf):
	htmlList = ["</br><div><h3>SNR Thresholds Per Marker</h3>","<table class=\"minimalistBlack\">",
			 "<tr><th>Marker</th><th>Poor Quality</th><th>Breakdown</th><th>Percentile Scatterplot</th><th>Otsu Scatterplot</th><th>Brenner Score</th></tr>"]
	template = """
	<tr>
	<td>{}</td>
	<td>{}</td>
	<td>{}</td>
	<td>{}</td>
	<td>{}</td>
	<td>{}</td>
	</tr>
	"""
	for mk in pDf['Marker'].unique():
		#mk = pDf['Marker'].unique()[2]
		print("  Assessing "+mk)
		dfSub = pDf[pDf['Marker']==mk]
		lostFovs = getCutoffImages(dfSub)
		pieChart = getQualityPieChart(dfSub)
		scatterPlot = getQualityScatterPlot(dfSub,"SNRp")
		scatterPlot2 = getQualityScatterPlot(dfSub,"SNRo")
		scatterPlot3 = getQualityScatterPlot(dfSub,"BrennerScore")
		htmlList.append(  template.format(mk, lostFovs, pieChart, scatterPlot, scatterPlot2, scatterPlot3)  )
	htmlList.append( "</table></div>" )
	return '\n'.join(htmlList)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Compile Percentile SNR metrics for assessment')
	parser.add_argument('-i', '--inputdir', help='Input OME.TIFFs directory.', nargs='?', type=str, dest="inDir",
						metavar="DIR", required=True)
	parser.add_argument('-o', '--outputdir', help='Output directory.', nargs='?', type=str, dest="outDir",
						metavar="DIR")
	#parser.add_argument('--stitched', action='store_true', default=False)
	parser.add_argument('--save_data', action='store_true', default=False)
	args = parser.parse_args()

	allFileDict = findAllOMEFiles(args.inDir, False)
	# TESTING
	#allFileDict = findAllOMEFiles("Q:\OME_TIFFs", False)
	#tmp1 = 'IH_37' #next(iter(allFileDict))
	#testing1 = { tmp1 : [allFileDict[tmp1][9]] }

	allDataStatsDF = findAllImageMetrics(allFileDict)
	#allDataStatsDF = findAllImageMetrics(testing1)


	if allDataStatsDF.empty:
		raise Exception('No Statistics Calculated!')

	#Generate Quality Cutoffs => Make Parameters in future
	allDataStatsDF['Quality'] = "Okay"
	allDataStatsDF.loc[(allDataStatsDF.ImageP10 < 2),'Quality']='Lower Bound Issue'
	allDataStatsDF.loc[(allDataStatsDF.SNRp <= 1.3),'Quality']='Questionable'
	allDataStatsDF.loc[(allDataStatsDF.SNRp <= 0.8),'Quality']='Poor'
	allDataStatsDF.loc[(allDataStatsDF.SNRo <= 1.3),'Quality']='Questionable'
	allDataStatsDF.loc[(allDataStatsDF.SNRo <= 0.8),'Quality']='Poor'
	allDataStatsDF.loc[(allDataStatsDF.NonBlankPercentage <= 0.06),'Quality']='Questionable'
	allDataStatsDF.loc[(allDataStatsDF.NonBlankPercentage < 0.02),'Quality']='Poor'

	# Add option to write Panda Dataframe out to csv
	if(args.save_data):
		outPath = "snr_data.csv"
		if args.outDir is not None:
			print( "outDir has been set (value is {})".format(args.outDir) )
			outPath = os.path.normpath(os.path.join(args.outDir, "snr_data.csv"))
		allDataStatsDF.to_csv(outPath)

	outPath = "qc_report.html"
	if args.outDir is not None:
		outPath = os.path.normpath(os.path.join(args.outDir, "qc_report.html"))
	Html_Out_file= open(outPath,"w")

	Html_Out_file.write(html_str_css)
	h1 = getSummaryHTMLTable(allDataStatsDF)
	Html_Out_file.write(h1)
	i2 = generateWholeBatchBoxplots_FOVs(allDataStatsDF)
	Html_Out_file.write(i2)
	i3 = generateWholeBatchSamples_FOVs(allDataStatsDF)
	Html_Out_file.write(i3)
	tbl = generateByMarkerSNRTable(allDataStatsDF)
	Html_Out_file.write(tbl)
	Html_Out_file.close()


