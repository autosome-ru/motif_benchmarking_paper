#!/usr/bin/env bash
MATRIX=$1
CENTERS_FOLDER=$2
echo -n $MATRIX $'\t'
java -cp ape-3.0.2.jar ru.autosome.macroape.ScanCollection all_hoco_motifs/${MATRIX}.pwm ${CENTERS_FOLDER} --all | grep -vPe '^#' | ruby -e 'puts readlines.sort_by{|l| l.split("\t")[1].to_f }.last.split("\t").first(2).join("\t")'