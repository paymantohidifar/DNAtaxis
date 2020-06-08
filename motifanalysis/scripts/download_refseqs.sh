#!/bin/bash

# set up variables for directories
rootdir=/home/payman/Research/dnataxis/motifanalysis
refdir=$rootdir/reference_data

# downaload Refseqs zip filed using 'assembly_accession' numbers of selected species
file=$refdir/speciesAcc.txt
while IFS= read -r line
do
  IFS=',' #setting comma as delimiter
  read -a strarr <<<"$line"
  taxid=${strarr[0]}
  acc=${strarr[1]}
  echo "Downloading Refseq for Tax ID: '$taxid'"
  datasets download assembly $acc --filename $refdir/$taxid.zip
done < "$file"
