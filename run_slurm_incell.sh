#!/bin/bash
#SBATCH --job-name=INCellMetericExtraction
#SBATCH --time=12:00:00
#SBATCH --mem=300G
#SBATCH --partition=cpu-med
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
#SBATCH --mail-user=moore.raymond@mayo.edu

INDIR="$1"
SLIDEINDX="$2"
WKFL=/research/bsi/projects/staff_analysis/m088378/MetericMaking

tmsp=$(date +"%D %T")
echo "Start: $SLIDEINDX [ $tmsp ]\n  $INDIR "

python $WKFL/grabAllPanelMetrics_InCell.py $INDIR $SLIDEINDX





