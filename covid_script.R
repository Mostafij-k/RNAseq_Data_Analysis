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
dim(filter_df1)


filter_df2=filter_df1[filter_df1$padj<.05,]
dim(filter_df2)

filter_df3=filter_df2[abs(filter_df2$log2FoldChange>1),]



#MA plot
plotMA(result_table)   #contain count value from RNAseq or microarry



plotMA( result_table, ylim = NULL,
        colNonSig = "gray32", colSig = "red3", colLine = "#ff000080")





#Volcano plot

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point() #simple way of volcano plot(){Basic}

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_vline(xintercept=1)+
  geom_vline(xintercept=-1)+
  geom_hline(yintercept = -log10(.05))               #Here add fold change cutoff (geom_vline=vertical line and hline=horizontal)


#Differentially express genes (significantly)

filter_df1$padj <.05 & abs (filter_df1$log2FoldChange) >1
sum(filter_df1$padj <.05 & abs (filter_df1$log2FoldChange) >1)
filter_df1$test=filter_df1$padj <.05 & abs (filter_df1$log2FoldChange) >1 # here add new column test


ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=test))+
  geom_vline(xintercept=1)+
  geom_vline(xintercept=-1)+
  geom_hline(yintercept = -log10(.05)) 


#change colour by manual

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=test))+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1)+
  geom_vline(xintercept=-1)+
  geom_hline(yintercept = -log10(.05)) 



#other change

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=test),size=1,alpha=.3)+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1)+
  geom_vline(xintercept=-1)+
  geom_hline(yintercept = -log10(.05)) #change giom point size


ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=test))+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1,colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept = -log10(.05), colour='green',linetype=3)  #here linetype(1)=solid line and linetype(2)=split line, (3)=Dot 



#add theme

ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(aes(colour=test))+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1,colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept = -log10(.05), colour='green',linetype=3) +
  xlim(-6,6)+
  ylim(0,10)+
  theme_bw()+
  theme(legend.position = 'top')



library(BiocManager)
install('plotly')
library(ggplot2)
library(plotly)
library(tibble)

#now add row in column
filter_df1=rownames_to_column(filter_df1, var='ensgene') #this tibble function create column from row

#now add gene id in volcano plot by using ggplotly()
g=ggplot(filter_df1, aes(x=log2FoldChange, y=-log10(padj),
                         name=ensgene))+
  geom_point(aes(colour=test),size=1, alpha=.3)+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1,colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept = -log10(.05), colour='green',linetype=3) +
  theme_bw()+
  theme(legend.position = 'none')     

ggplotly(g)


install('biomaRt')
library(biomaRt)
??biomaRt
listMarts()
ensembl99=useEnsembl(biomart = 'ensembl', version=99)
View=listDatasets(ensembl99)
ensembl99=useDataset('hsapiens_gene_ensembl', 
                      mart=ensembl99)

listAttributes(ensembl99)
View(listAttributes(ensembl99))
getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version', 'ensembl_transcript_id',
                     'ensembl_transcript_id_version', 'external_gene_name'),
      filters = c('ensembl_gene_id'),
      values = filter_df1$ensgene[1:6],
      mart=ensembl99)
#annotation
annotation=getBM(attributes = c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                 filters = c('ensembl_gene_id'),
                 values=filter_df1$ensgene,
                 mart=ensembl99)

#join data frame 
annotated_df=left_join(filter_df1,annotation, by=c('ensgene'='ensembl_gene_id')) #by left join dplyr fuction.
View(annotated_df)

#now volcano plot for annotated_df 

g=ggplot(annotated_df, aes(x=log2FoldChange, y=-log10(padj),
                           name=external_gene_name))+
  geom_point(aes(colour=test),size=1, alpha=.3)+
  scale_colour_manual(values=c('black','red'))+
  geom_vline(xintercept=1,colour='green', linetype=3)+
  geom_vline(xintercept=-1, colour='green',linetype=3)+
  geom_hline(yintercept = -log10(.05), colour='green',linetype=3) +
  theme_bw()+
  theme(legend.position = 'none')

ggplotly(g)

# Add gene name , here use pvalue, you can use padj
p <- ggplot(data=annotated_df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
ggplotly(p)

# Add more simple "theme"
p <- ggplot(data=annotated_df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()


#Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
ggplotly(p2)



# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
annotated_df$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP" 
annotated_df$diffexpressed[annotated_df$log2FoldChange > 1 & annotated_df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
annotated_df$diffexpressed[annotated_df$log2FoldChange < -1 & annotated_df$pvalue < 0.05] <- "DOWN"



# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=annotated_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

ggplotly(p2)



## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

ggplotly(p3)


# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
annotated_df$delabel <- NA
annotated_df$delabel[annotated_df$diffexpressed != "NO"] <- annotated_df$external_gene_name[annotated_df$diffexpressed != "NO"]

ggplot(data=annotated_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=annotated_df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")









#heatmap
?heatmap
anno_df2=annotated_df[annotated_df$pvalue<.05,]
anno_df3=anno_df2[abs(anno_df2$log2FoldChange)>1,]
degs=anno_df3$ensgene
#vst_nhbe=varianceStabilizingTransformation(dds_nhbe)
#vst_nhbe_mat=assay(vst_nhbe)

data_for_hm=vst_nhbe_mat[degs,]
rownames(data_for_hm)=anno_df3$external_gene_name
heatmap(data_for_hm)

#To improve heatmap use pheatmap packages
install('pheatmap')
library(pheatmap)
pheatmap(data_for_hm)

?pheatmap()

pheatmap(data_for_hm, fontsize = 3, scale='row')

##Install RcolorBrewer
install('RColorBrewer')
library(RColorBrewer)
display.brewer.all()
brewer.pal(9,'Greys')  #To know hex value

pheatmap(data_for_hm, fontsize = 4, scale='row')

colorRampPalette(brewer.pal(9,'Greys'))(100)
pheatmap(data_for_hm, fontsize = 4, scale='row', color=greys)  # problem

pheatmap(data_for_hm, fontsize = 4, scale='row', color=greys)

pairs= colorRampPalette(brewer.pal(12,'paired'))(100)
pheatmap(da_for_hm, fontsize_row=4, scale='row', color=pairs)
last_scheme = colorRampPalette(brewer.pal(7,'PuRd'))(100)
pheatmap(data_for_hm,fontsize_row=4,scale='row', color=last_scheme)



#gene ontology
library(BiocManager)
install('clusterProfiler')
library(clusterProfiler)

install('org.Hs.eg.db')
library(org.Hs.eg.db)


#go enrichment
??enrichGO

ent_gene =getBM(attributes=c('entrezgene_id'),
                filters=c('ensembl_gene_id'),
                values=anno_df3$ensgene,
                mart= ensembl99)
ent_gene=ent_gene$entrezgene_id
ent_gene=as.character(ent_gene)


ent_uni= getBM(attributes=c('entrezgene_id'),
               filters=c('ensembl_gene_id'),
               values=annotated_df$ensgene,
               mart=ensembl99)
ent_uni=as.character(ent_uni$entrezgene_id)
ent_uni=as.character(ent_uni)

ego=enrichGO(gene=ent_gene,
             OrgDb = org.Hs.eg.db,
             ont="BP",
             universe =ent_uni )


library(enrichplot)
??clusterProfiler
barplot(ego)
barplot(ego, showCategory=20)
dotplot(ego)
dotplot(ego,showCategory=20)

??enrichKEGG
ekg = enrichKEGG(gene =ent_gene, universe=ent_uni)
View(summarry(ekg))

# Save
write_tsv(anno_df3, 'filtered_nhbe_results.txt')
















