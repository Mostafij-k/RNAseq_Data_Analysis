***RNA-seq Data Analysis Pipeline***

This repository contains a complete step-by-step workflow for RNA sequencing (RNA-seq) data analysis, starting from raw data retrieval to downstream differential expression and visualization.

📌**Overview**
  
    ☀ Downloading raw sequencing data from GEO/SRA
    
    ☀ Processing raw FASTQ files using Salmon
    
    ☀ Generating transcript-level quantification
    
    ☀ Converting transcript-level data to gene-level counts
    
    ☀ Performing differential gene expression analysis using DESeq2
    
    ☀ Visualizing results using heatmaps, volcano plots, and other figures
  
🔬 **Workflow Summary**

    ☀ Data Acquisition
    
      > Data was downloded from GEO databases
      
      > SRA Toolkit (fasterq-dump) and SRA Run Selector were used to obtain FASTQ files.
  
    ☀ Quality and Quantification
  
      > Fustqc use to check the quality of data
      > Transcript-level quantification was performed using Salmon
  
    ☀ Quality and Quantification
  
      > R package tximport was used to gene level count
  
    ☀ Differential Expression Analysis
    
      > Gene-level count matrix was analyzed using DESeq2 in R
  
    ☀ Data Visualization
    
      > Heatmap, Volcano Plot were draw to visualize the DEG

📊 **Tools & Technologies**

    ✔ SRA Toolkit
    ✔ Salmon
    ✔ R Programming
    ✔ Bioconductor packages:
    ✔ DESeq2
    ✔ tximport
    ✔ ggplot2
    ✔ Linux shell scripting



    
