library(ggplot2)
library(plotly)

d=read.delim(file = 'MC_tab.txt', header = TRUE, sep='\t')

d1=as.data.frame(d)

View(d1)

# Remove NA from data frame

#complete.cases(d1)  # Show logical vector IF contain NA give FALSE

#sum(complete.cases(d1))

dg=d1[complete.cases(d1),]
View(dg)

#Volcano plot
library(ggplot2)
library(plotly)

ggplot(dg, aes(x=logFC, y=-log10(adj.P.Val)))+
  geom_point() #simple way of volcano plot(){Basic}



# add a column of NAs
dg$Category <- "Not Sig"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 

dg$Category[dg$logFC > 2 & dg$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"

dg$Category[dg$logFC < -2 & dg$adj.P.Val < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"

p <- ggplot(data=dg, aes(x=logFC, y=-log10(adj.P.Val), col=Category)) +
  geom_point(size = 1.5, alpha =1 , na.rm = T, shape = 21)



# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-2,2), col="red",linetype="dashed") +
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  theme(legend.position = "right") +
  labs(x="log2 fold change",
       y="-Log10 adj p value",
       title="GSE15605(MC)")+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)+
  ylim(0,10)   
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
  theme_bw(base_size = 25) +         # theme_bw will create border (fixed size)
  theme(legend.position = "right") +
  labs(x="log2 fold change",
       y="-Log10 adj p value",
       title="GSE15605(MC)")+
  geom_hline(yintercept=-log10(0.05), col="red", linetype="dashed")+
  xlim(-6,6)+
  scale_y_continuous(trans="log1p", breaks = c(0,10,2600))   # 
ggplotly(p2)


p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "Not Sig")
p3 <- p2 + scale_colour_manual(values = mycolors)
ggplotly(p3)












