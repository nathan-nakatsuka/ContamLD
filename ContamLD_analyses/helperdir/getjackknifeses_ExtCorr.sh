#!/bin/bash

while read sampleID panel; do
cd $2/directories/${sampleID}_jackknife

var=$(<${sampleID}_${panel}_$4panel_scores_reads_ExtCorr_meanscore.txt)
$1/helperdir/dowtjack -i ${sampleID}_${panel}_$4panel_scores_reads_ExtCorr.jin -m ${var} -o ${sampleID}_${panel}_$4panel_scores_reads_ExtCorr.jout

done < $2/$3
