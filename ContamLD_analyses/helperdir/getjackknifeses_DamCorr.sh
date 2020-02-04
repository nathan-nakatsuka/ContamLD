#!/bin/bash

while read sampleID panel; do
cd $2/directories/${sampleID}_jackknife

var=$(<${sampleID}_${panel}_$4panel_scores_reads_DamCorr_meanscore2.txt)
$1/helperdir/dowtjack -i ${sampleID}_${panel}_$4panel_scores_reads_DamCorr.jin -m ${var} -o ${sampleID}_${panel}_$4panel_scores_reads_DamCorr.jout

var2=$(<${sampleID}_${panel}_$4panel_kscores_reads_DamCorr_meanscore2.txt)
$1/helperdir/dowtjack -i ${sampleID}_${panel}_$4panel_kscores_reads_DamCorr.jin -m ${var2} -o ${sampleID}_${panel}_$4panel_kscores_reads_DamCorr.jout

done < $2/$3
