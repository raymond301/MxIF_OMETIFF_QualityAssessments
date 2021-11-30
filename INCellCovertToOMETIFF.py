# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:23:12 2020
@author: Raymond Moore
"""
import os, sys, glob, random, argparse
import tifffile
import xarray as xr
import numpy as np
from pprint import pprint
from PIL import Image, ImageFont, ImageDraw

Image.MAX_IMAGE_PIXELS = None
import cv2
from subprocess import check_output
# pip install imagecodecs -user

#####################
###   FUNCTIONS   ###
#####################

def get_marker_name(filehandle):
	if 'RegisteredImages' in filehandle and 'dapi' in filehandle:
		return 'DAPI'
	elif 'AFRemoved' in filehandle:
		return os.path.basename(filehandle).split("_AFRemoved")[0]
	else:
		return '???'


def get_channel_name(paneljson, mark):
	if mark == 'DAPI':
		return mark
	else:
		phld = 'null'
		for k, val in paneljson.items():
			if mark in val:
				idx = val.index(mark)
				if idx == 0:
					phld = 'Cy2'
				elif idx == 1:
					phld = 'Cy3'
				elif idx == 2:
					phld = 'Cy5'
				else:
					phld = 'Unlabled'
		return phld


def get_channel_xml(i, mark, cy):
	outStr = f"""<Channel ID="Channel:{i}" Name="{mark}" Fluor="{cy}" SamplesPerPixel="1" ContrastMethod="Fluorescence" />"""
	return outStr


def write_ometiff(inDir, panel, sample, fov, outpath):
	### Get Primary DAPI IMAGE
	topDapi = glob.glob(os.path.join(inDir, 'RegisteredImages', 'S002', "*mono_dapi_*{}*.tif".format(fov)))
	if len(topDapi) < 1:
		raise Exception('Unable to find DAPI TIFF file!')
	if not os.path.exists(topDapi[0]) or not os.path.getsize(topDapi[0]) > 0:
		raise Exception('Missing/Empty Original DAPI TIFF file!')
	### Get All other AFRemoval Images
	afRemovs = sorted(glob.glob(os.path.join(inDir, 'AFRemoved', "*AFRemoved_*{}.tif".format(fov))))
	allTiffs = topDapi + afRemovs
	dat = []
	xmlChan = []
	cDim = []
	u = 0
	outChannelTxt = open(outpath.replace('.ome.tiff', '.channels.txt'), "w")
	for fh in allTiffs:
		#print( "\t\tLoading file: {}".format(os.path.basename(fh)) )
		## Check to ensure import TIFF is single image type.
		tmpImg = tifffile.TiffFile(fh)
		## Get MetaData
		mark = get_marker_name(fh)
		channel = get_channel_name(panel, mark)
		cDim.append(mark)
		print("    {}: {}  {} Image: {} {}".format(u, mark, channel, str(tmpImg.pages[0].shape),
												   os.path.basename(fh)))
		outChannelTxt.write('\t'.join([str(u), mark, channel, os.path.basename(fh)]) + "\n")
		xmlChan.append(get_channel_xml(u, mark, channel))
		# dat.append([tmpImg.pages[0].asarray(),tmpImg.pages[1].asarray()])
		# pprint(tmpImg.pages[0])
		dat.append(xr.DataArray(tmpImg.pages[0].asarray(), name=mark, dims=("y", "x")))
		u += 1
	outChannelTxt.close()

	imarr = xr.concat(dat, dim="c")
	imarr.assign_coords({"c": cDim, "x": range(2040), "y": range(2040)})

	# imarr = imarr.transpose("c", "y", "x")
	Nc, Ny, Nx = imarr.shape
	pprint(["Nc", Nc, "Ny", Ny, "Nx", Nx])
	channels_xml = '\n'.join(xmlChan)
	outname = os.path.splitext(os.path.basename(outpath))[0]
	xml = f"""<?xml version="1.0" encoding="UTF-8"?>
	<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
			xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
			xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
		<Image ID="Image:0" Name="{outname}">
			<Pixels BigEndian="false"
					DimensionOrder="XYCZT"
					ID="Pixels:0"
					Interleaved="false"
					SizeC="{Nc}"
					SizeT="1"
					SizeX="{Nx}"
					SizeY="{Ny}"
					SizeZ="1"
					PhysicalSizeX="0.37"
					PhysicalSizeY="0.37"
					Type="float">
				<TiffData />
				{channels_xml}
			</Pixels>
		</Image>
	</OME>
	"""

	print("    Number of Channels = " + str(len(xmlChan)))
	tifffile.imwrite(outpath, data=imarr.values, description=xml, contiguous=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Generate OME TIFF files from INCell directory.')
	parser.add_argument('-i', '--inputdir', help='Input INCell Raw directory.', nargs='?', type=str, dest="inDir",
						metavar="DIR", required=True)
	parser.add_argument('-o', '--outputdir', help='Output OME.TIFFs directory.', nargs='?', type=str, dest="outDir",
						metavar="DIR", required=True)
	args = parser.parse_args()

	sample = os.path.basename(os.path.dirname(args.inDir))
	print("OUTPUT:"+args.outDir)
	if not os.path.isdir(os.path.join(args.outDir, sample)):
		print("Output Dir doesn't exist. Creating now...")
		os.mkdir(os.path.join(args.outDir, sample))

	panel = {}
	panelDirs = glob.glob(os.path.join(args.inDir, "*_dapi_dapi_*"))
	if len(panelDirs) < 1:
		print("No panel found!")
		sys.exit()
	for dd in panelDirs:
		nom = os.path.basename(dd)
		if "HLA" in nom:
			nom = nom.replace('HLA_', 'HLA-')
		if "Caspase" in nom:
			nom = nom.replace('Caspase_', 'Caspase-')
		if "isotype" in nom:
			nom = nom.replace('isotype_', 'isotype-')
		dBase = nom.split("_")

		if dBase[4] == "fitc":
			panel[dBase[0]] = [dBase[3], dBase[5], dBase[7]]
		elif dBase[4] == "cy3" and len(dBase) > 5:
			panel[dBase[0]] = ['-', dBase[3], dBase[5]]
		elif dBase[4] == "cy3":
			panel[dBase[0]] = ['-', dBase[3], '-']
		elif dBase[4] == "cy5":
			panel[dBase[0]] = ['-', '-', dBase[3]]
		else:
			pprint(dBase)

	panel2 = dict([(k, [x.upper() for x in v]) for (k, v) in panel.items()])
	pprint(panel2)

	totReg = len(glob.glob(os.path.join(args.inDir, 'RegisteredImages', 'S002', "*mono_dapi_*.tif")))
	print("> Sample {} contains {} regions".format(sample, totReg))
	if not os.path.isdir(os.path.join(args.outDir, sample)):
		os.mkdir(os.path.join(args.outDir, sample))

	for nn in range(totReg):
		fov = 'region_' + str(nn + 1).zfill(3)
		outpath = os.path.join(args.outDir, sample, "{}_{}.ome.tiff".format(sample, fov))

		print("Loading: {} - {} to {}".format(sample, fov, outpath))
		write_ometiff(args.inDir, panel2, sample, fov, outpath)
