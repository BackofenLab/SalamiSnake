#!/bin/bash

INPUT_FOLDER=$1"/"
OUTPUT_FOLDER=$2"/"
INPUT_FORMAT=$3

numfiles=$((-1))

# get all files in input folder 
declare -a FILELIST
for f in $INPUT_FOLDER*.$INPUT_FORMAT; do 
	numfiles=$(($numfiles+1))
    FILELIST=("${FILELIST[@]}" "$f")
done

# pairwise intersections
for i in `seq 0 $numfiles`; do
    for j in `seq $(($i+1)) $numfiles`; do
    	a=${FILELIST[$i]}
    	a=${a//.$INPUT_FORMAT/}
    	a=${a//$INPUT_FOLDER/}
    	b=${FILELIST[$j]}
    	b=${b//.$INPUT_FORMAT/}
    	b=${b//$INPUT_FOLDER/}
    	outfile=$OUTPUT_FOLDER$a"_vs_"$b".bed"
		bedtools intersect -a ${FILELIST[$i]} -b ${FILELIST[$j]} -s -u -f 0.1 > $outfile
    done
done

#check if there are only two files 
if [ $(($numfiles)) \> 1 ]; then
    # between all intersections
    outfile_name=$OUTPUT_FOLDER"robust_between_all"
    a=${FILELIST[0]}
    for i in `seq 1 $numfiles`; do
    	outfile=$outfile_name"_"$i".bed"
    	b=${FILELIST[$i]}
    	echo $a
    	echo $b
    	bedtools intersect -a $a -b $b -s -u -f 0.1 > $outfile
    	a=$outfile
    done

    mv $a $OUTPUT_FOLDER"robust_between_all.bed"
    rm $OUTPUT_FOLDER*"robust_between_all_"*".bed"
else
    mv $outfile $OUTPUT_FOLDER"robust_between_all.bed"
fi