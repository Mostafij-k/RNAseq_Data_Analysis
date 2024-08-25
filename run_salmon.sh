#!/bin/bash

# Read input
GEO=$(cat geo_accessions.txt)

for i in ${GEO}
    do
       SRR=$(grep ${i} SraRunTable.txt | cut -d ',' -f 1)
       #echo ${SRR}
       SRR=$(echo ${SRR} | sed 's/ /.fastq.gz /g')
       #echo ${SRR}
       SRR=${SRR}.fastq.gz
       #echo ${SRR}
       salmon quant -i gencode_v38_index -l A -o ${i} -r ${SRR}
    done 

