#!/bin/bash

# set up variables for directories
rootdir=/home/payman/Research/dnataxis/motifanalysis
metadir=$rootdir/meta
refdir=$rootdir/reference_data

# get number of species
count=$(cat "$metadir/speciesNames.txt" | wc -l)

# retrieve description of Refseq for each species and write them into a file
for ((i=1; i<=count; i++)); do
    name=$(sed -n "$i"p "$metadir/speciesNames.txt")
    datasets assembly_descriptors tax_name "$name" -r | \
    jq-linux64 '.datasets[] | {scientific_name: .org.sci_name, strain: .org.strain,
    tax_id: .org.tax_id, Assembly_level: .assembly_level, Assembly_accession: .assembly_accession,
    sequence_length: .seq_length,  Submission_date: .submission_date}' >> $refdir/refSeqInfo.txt
    echo "--------------------------" >> $refdir/refSeqInfo.txt
done