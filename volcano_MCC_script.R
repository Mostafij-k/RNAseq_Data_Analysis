library(affy)
library(oligo)
library(Biobase)
library(GEOquery)
library(splitstackshap)
library(tidyr)
library(dplyr)
library(arrayQualityMetrics)
library(limma)
library(plotly)

#load series and platform data from GEO

gset <- getGEO("GSE39612", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXXXXX111111111111111111111111111111111100000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
hist(gset)

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("control","Case"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, .05)
data <- topTable(fit2, adjust="fdr", sort.by="B", number=60000)

#data2 <- subset(data, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))

data2=as.data.frame(data)

View(result_df)

# Remove NA from data frame

#complete.cases(result_df)  # Show logical vector IF contain NA give FALSE

#sum(complete.cases(result_df))

#filter_df1=result_df[complete.cases(result_df),]
View(filter_df1)

#Volcano plot
library(ggplot2)
library(plotly)

ggplot(data2, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point() #simple way of volcano plot(){Basic}



# add a column of NAs
data2$Category <- "Not Sig"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 

data2$Category[data2$logFC > 2 & data2$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"

data2$Category[data2$logFC < -2 & data2$adj.P.Val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

p <- ggplot(data=data2, aes(x=logFC, y=-log10(adj.P.Val), col=Category)) +
  geom_point(size = 1.5, alpha =1 , na.rm = T, shape = 21)

ggplotly(p)

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-2,2), col="red",linetype="dashed") +
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  theme(legend.position ='right')+
  labs(x="log2 fold change",
       y="-Log10 adj p value",
       title="GSE39612(MCC)")+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)+
  ylim(0,80)
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









