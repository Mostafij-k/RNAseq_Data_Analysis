#!/bin/bash

# Loop for downloading data
SRA=$(tail -n +2 SraRunTable.txt | cut -d "," -f 1)

for i in ${SRA}
   do
      prefetch ${i}
      if [ -f ${i}.fastq.gz]
          then 
                  echo "${i} is already Finished"
       else
                  echo '(o) Convering  SRA entry: ${i}'
                  fastq-dump --gzip --defline-qual '+' ${i}/${i}.sra
                  echo '(o) Done convert SRA entry'
       fi
  done
