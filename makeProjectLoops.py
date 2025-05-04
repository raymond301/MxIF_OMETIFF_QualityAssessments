import os,sys,re,glob
import subprocess

base=r'/research/bsi/archive/PI/'


def getPIArchiveFolder(pi, di):
	return os.path.join(base,pi,'tertiary',di,'integrated','*')

if __name__ == '__main__':
	potentialProjects = []
	potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Markovic_Svetomir_snm02', 's210155.CellSegmentation') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Ilyas_Sumera_m052774', '300919.Human_Cholangiocarcinoma_MxIF') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Yoon_Harry_m061620', 's210155.MxIF_MarkovicCollaboration') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Block_Matthew_msb04', 's210155.MxIF_MarkovicCollaboration') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Wang_Chen_m092469', 's302493.MxIF_Ovarian_Cancer') ) )))
	#potentialProjects.extend(list(filter(os.path.isdir, glob.glob( getPIArchiveFolder('Khazaie_Khashayarsha_m123285', '300679.Collaboration_MxIF') ) )))
	
	for pProj in potentialProjects:
		dfList = []
		ppNom = os.path.basename(pProj)
		print(" Project: "+ppNom)
		
		targetSamples = list(filter(os.path.isdir, glob.glob( os.path.join(pProj,'*')) ))
		print(" Samples: "+str(len(targetSamples))+"    (Y/n)")

		prjtFlags=input()
		if prjtFlags == "n":
			print("        Skip.")
			continue
		else:
			for sample in targetSamples:
				sNom = os.path.basename(sample)
				#print("    "+sNom)
				if sNom.startswith("SLIDE"):
					afiFiles = list(filter(os.path.isfile, glob.glob( os.path.join(sample,sNom+"_Final","CellDIVE_"+sNom+"_R*.afi") ) ))
					print("  > CellDive:   {} => {} AFIs   (Y/n)".format(sNom, len(afiFiles)))
					lstFlags=input()
					if lstFlags == "n":
						print("        Skip.")
						continue
					else:
						for afiFh in afiFiles:
							roiIdx = os.path.basename(afiFh).replace(".afi","").split("_")[-1]
							# print("run_slurm_celldive.sh {} {}".format(sample, roiIdx))
							subprocess.Popen("sbatch run_slurm_celldive.sh {} {}".format(sample, roiIdx), shell=True)
							#sys.exit()
				else:
					totReg = list(filter(os.path.isfile, glob.glob( os.path.join(sample,'RegisteredImages','S002',"*mono_dapi_*.tif"))))
					totReg = [os.path.basename(e).replace(".tif","").split("_")[-1] for e in totReg]
					print("  & InCell:  {} => {} FOVs   (Y/n)".format(sNom, len(totReg)))
					lstFlags=input()
					if lstFlags == "n":
						print("        Skip.")
						continue
					else:
						for afRegFh in totReg:
							roiIdx = os.path.basename(afRegFh).replace(".tif","").split("_")[-1]
							# print("run_slurm_incell.sh {} {}".format(sample, roiIdx))
							subprocess.Popen("sbatch run_slurm_incell.sh {} {}".format(sample, roiIdx), shell=True)
							#sys.exit()

					
					
				
				
				
