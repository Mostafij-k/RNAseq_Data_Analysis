#!/bin/bash

#get input data (uniport ids)
VAR=$(tail -n +2 SraRunTable.txt | cut -d "," -f 1)

for i in ${VAR}
   do
       if [ -f ${i}.fastq.gz ]
          then 
            echo "${i} is already downloded"
       else
            echo '(o) Downloading SRA entry: ${i}'
            fastq-dump --gzip --defline-qual '+' ${i}
            echo '(o) Done Download SRA entry'
       fi
  done
