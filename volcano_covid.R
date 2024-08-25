getwd()
dplyr::select()  # it is useful to call select function

## Call all necessary packages
install.packages('dplyr')
library(dplyr)
library(BiocManager)  # Call BiocManager to download bioconductor packages
library(DESeq2)
library(tximport)
library(ggplot2)
library(readr)
library(dplyr)
library(magrittr)  #Load this packages for piping
library(biomaRt)
library(plotly)
library(tidyverse)    

## Filter metadata for further analysis from SraRunTable.txt
sample_table = read_csv("SraRunTable.txt") %>%
  select(`Sample Name`, severity, source_name) %>%
  slice(seq(1,68, by=2))
## Need produce files for tximport()
## That's why we want to create 'GSM4432378/quant.sf' figure

sample_files=paste0(pull(sample_table,
                         `Sample Name`), '/quant.sf')

names(sample_files) = pull(sample_table, `Sample Name`)

## Need gene_map csv file for tx2gene contain enst and ensg id 

gene_map= read_csv("gene_map.csv", col_names = c('enstid', 'ensgid'))

## Summarize gene level from transcript level count produce by salmon

count_data=tximport(files = sample_files,
                    type = 'salmon',
                    tx2gene =gene_map ,
                    ignoreTxVersion = TRUE)     # Here ignoretxversion (ignore ENST2676767677.2) .2

## Produce more friendly data for DESeq2     

sample_table = as.data.frame(sample_table)  # Sample_table change into data frame from tibble

colnames(sample_table)[1]= "Sample"         # It will change column header

sample_table

conditions= c('treatment', 'condition')

conditions=rep(conditions, each=17)  # As we have 12 sample

conditions= factor(conditions)   # It produce factors 

sample_table$conditions = conditions  # Add experimental conditions in sample_table as column

deseq_dataset= DESeqDataSetFromTximport(txi=count_data, 
                                        colData =sample_table ,
                                        design =~conditions )     # Creates deseq data set from tximport

counts(deseq_dataset)   # Extract to show



deseq_dataset =estimateSizeFactors(deseq_dataset)  #Step 1

deseq_dataset = estimateDispersions(deseq_dataset)  #Step 2

plotDispEsts(deseq_dataset)

deseq_dataset = nbinomWaldTest(deseq_dataset)  # Step 3


result_table=results(deseq_dataset)

summary(result_table)
View(as.data.frame(result_table))
# Result_table is a data.frame!

result_df=as.data.frame(result_table)  # Produce data frame

View(result_df)

# Remove NA from data frame

complete.cases(result_df)  # Show logical vector IF contain NA give FALSE

sum(complete.cases(result_df))

filter_df1=result_df[complete.cases(result_df),]
View(filter_df1)

#Volcano plot
library(ggplot2)
library(plotly)

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_point() #simple way of volcano plot(){Basic}



# add a column of NAs
filter_df1$Category <- "Not Sig"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 

filter_df1$Category[filter_df1$log2FoldChange > 2 & filter_df1$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"

filter_df1$Category[filter_df1$log2FoldChange < -2 & filter_df1$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

p <- ggplot(data=filter_df1, aes(x=log2FoldChange, y=-log10(pvalue), col=Category)) +
  geom_point(size = 1.5, alpha =1 , na.rm = T, shape = 21)


# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-2,2), col="black",linetype="dashed") +
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  labs(x="log2 fold change",
       y="-Log10 pvalue",
       title="GSE152418(COVID-19)")+
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed")+
  xlim(-6,6)+
  ylim(0,20)
  #scale_y_continuous(trans="log1p", breaks=c(0,30,3000))   # here you can use 2600


ggplotly(p2)

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "Not Sig")
p3 <- p2 + scale_colour_manual(values = mycolors)
ggplotly(p3)



# Or
p2 <- p + geom_vline(xintercept=c(-2,2), col="red",linetype="dashed") +
  theme_bw() +         # theme_bw will create border (fixed size)
  theme(legend.position = 'none')+
  theme(legend.position = "right") +
  ggtitle(label = 'Differential Expression- COVID-19')+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)
  #scale_y_continuous(trans="log1p", breaks=c(0,30,2600))


ggplotly(p2)






















