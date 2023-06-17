#First set your working directory:
setwd("C:/Mi unidad/Analisis/DEG")

#Then call all the libraries from the RTCGA package. They allows us to download and integrate the variety and volume of TCGA data.
library("TCGAbiolinks")
library("RTCGA")
library("RTCGA.clinical")
llibrary("SummarizedExperiment")

#If you don't have them already installed you must do it using the BiocManager. (https://rtcga.github.io/RTCGA)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RTCGA")

#Then call the packages for the DEG analysis. I'm using edgeR which also needs limma
library(limma)
library(edgeR)

#Last but not least, the libraries for the manipulation of the data
library("tidyverse") #contain tidyr, ggplot2 and dplyr, among others

#Download the RNAseq
query.exp <- GDCquery(
  project = "TCGA-SKCM", # I use the SKCM project, you may use BRCA or breast cancer or for MMRF bone marrow, etc.
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts") # For raw counts

GDCdownload(
  query = query.exp,
  files.per.chunk = 100)

skcm.exp <- GDCprepare(
  query = query.exp,
  save = TRUE,
  save.filename = "skcmExp.rda")

rnaseq <- assay(skcm.exp)
rnaseq <- rnaseq[which(!duplicated(rownames(rnaseq))),]   
write.table(rnaseq, "raw counts.txt") # To save the data set

# Manipulate an clean the data set
# - Mutate the ids 
id<-as.data.frame(colnames(rnaseq))
colnames(id)<-c("bcr_patient_barcode")
id<-id%>%mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
colnames(rnaseq)<-id$bcr_patient_barcode

# - Depending on the package you may need ENTREZ, GENE SYMBOL or ENSEMBL 
library(org.Hs.eg.db) # To acces to all the gene notations
hs <- org.Hs.eg.db 
my.symbols <- c(rnaseq$ENSEMBL)
IDS<-select(hs, 
            keys = my.symbols,
            columns = c("ENSEMBL", "SYMBOL","ENTREZID"),
            keytype = "ENSEMBL")
IDS<-na.omit(IDS)
rnaseq<-merge(rnaseq,IDS,by="ENSEMBL")
write.table(rnaseq,"para edgeR.txt")

# - Clinical information
clinica<-as.data.frame(SKCM.clinical)
clinica$patient.bcr_patient_barcode<-toupper(clinica$patient.bcr_patient_barcode) 
colnames(id)<-c("patient.bcr_patient_barcode")
x<-merge(clinica,id,by="patient.bcr_patient_barcode")
x<-x[,c(1,9,11,13,16,28,30,31,438,439,448,449,450,770,775,778,779,780,785,917,935,939,940,941,944,945,846,949,956,964,1060,1867)]
write.table(x,"clinical information.txt")

# - Categorize according to Vav1,2 y 3 expression
phenotype<-as.data.frame(rnaseq[rnaseq$SYMBOL=="VAV1",])
phenotype[2,]<-rnaseq[rnaseq$SYMBOL=="VAV2",]
phenotype[3,]<-rnaseq[rnaseq$SYMBOL=="VAV3",]
rownames(phenotype)<-phenotype$SYMBOL
phenotype[,1]<-NULL # to remove the first column, the ENSEMBL
phenotype[,474:475]<-NULL # to remove the first column, the SYMBOL and ENTREZID
phenotipe<-cpm(phenotipe)

# To perform the DEG analysis you should order the patients according to the expression of your protein:
secuencias <- secuencias[order(secuencias$VAV3),] # In my case I study VAV3.
ICVAV3<-secuencias[-(93:276),] # This is because I wanted to work with Q1 and Q3. My dataset had 368 patients, then became 188.
write.table(ICVAV3, "ICVAV3.txt") # I'm quite strict about saving every little step.

# Then you have to make the condition vector, I have two conditions for my protein: High and Low, you may have more conditions. In this case 1 is for Low expression and 2 is for High.
condition<-c(rep(1,92),rep(2,92))
description<-data.frame(barcode,condition) # This may be helpfull for future analysis, you'll have your patients categorized.

#The DEG itself
tRNASeq_VAV3<-data.frame(t(ICVAV3)) # You'll need to transpose the data
y <- DGEList(counts=tRNASeq_VAV3,group=condition) # Here you'll be creating the DGElist object.
keep <- filterByExpr(y) # In this step, the edgeR algorithm will be keeping the rows that are 'worthwhile' keeping (doesn't have too many zeros).
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~condition) # Your design comes from the condition previously set. Here you can add another interaction using ~ condition A + condition B
y0 <- estimateDisp(y,design) # In this step I change the name because in some analysis the estimateDisp function isn't needed.
et <- exactTest(y0) # Here the pseudocounts matrix is made.
is.de <- decideTestsDGE(et) # We examin the total nmuber if DEG.
summary(is.de)
df_et<-as.data.frame(topTags(et,n=Inf))
dim(df_et) # Just for control
head(df_et)
df_et$Expression = ifelse(df_et$FDR < 0.01 & abs(df_et$logFC) >= 0.6, 
                               ifelse(df_et> 0.6 ,'Up','Down'),
                               'Stable') # I add a column with the expression for furhter use.
df_et$Genes = rownames(df_et) # Prepare for saving the data.
write.table(df_et, "Exact test (ICVAV3).txt") 

#Performing the VULCANO PLOT to visualize de Genes and its expression
plot <- ggplot(data = df_et, 
               aes(x = logFC, 
                   y = -log10(FDR), 
                   colour=Expression
               )) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-6,6)) +
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 2.001,lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Differential expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",
        legend.title = element_blank())
plot

# To perform the ORA using GO, KEGG (both with edgeR) or REACTOME (with ReactomPA) you first need to change the Genes ID to Entrez.
## So, call the libraries
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db) # In my case the organism is Human so I use Hs, to use otrher species check the annotation.fit <- glmFit(y0, design)

## Creating the new data
hs <- org.Hs.eg.db
Exact_test_ICVAV3<-subset(df_et,df_et$Expression!="Stable") # I'll be only using the up and down DEGs.
my.symbols <- c(Exact_test_ICVAV3$Genes) # The string with the genes names to map with ENTREZ
IDS<-AnnotationDbi::select(org.Hs.eg.db, keys=my.symbols, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
names_ids<-IDS$ENTREZID
b<-na.omit(IDS) # I've never had na but, just to make sure.
c<-merge(tRNASeq_SNCG,b,by="SYMBOL") # Here I create my new data as the package requires.
c$SYMBOL<-NULL  
row.names(c)<-c$ENTREZID
c$ENTREZID<-NULL

## Finally, the enrichment analysis
x <- enrichPathway(gene=names_ids,pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory=30, title="DotPlot VAV3 DEG enrichment analysis (IC)",font.size=6) # You can also create a CNet Plot and Bar Plot

#To know the KEGG pathway enriched in this contrast
contrast <- makeContrasts(2 - 1, levels=design) # I make the constrast, in this case High/Low
fit<-glmQLFit(y0,design,robust=TRUE)
qlf <- glmQLFTest(fit, contrast=contrast)
keg <- kegga(qlf, species="Hs") # Pay atention to the species set, I'm working whith Homo sapiens
KEGG_ICVAV3<-data.frame(topKEGG(keg)) # I'll create a data frame of the most significative KEGG pathways
write.table(KEGG_ICVAV3, "KEGG (ICVAV3).txt") # Allways save my work.
(barcode,condition)

# To know the GO terms that are enriched, the algorithm uses the same "fit" already created.
lrt <- glmLRT(fit)
go <- goana(lrt)
goICVAV3_bp<-topGO(go, ont="BP", sort="Up", n=100, truncate=30) # for biological process
goICVAV3_cc<-topGO(go, ont="CC", sort="Up", n=100, truncate=30) # for cellular component
goICVAV3_mf<-topGO(go, ont="MF", sort="Up", n=100, truncate=30) # for molecular functions
write.table(goVAV3_bp, "go_bp (ICVAV3).txt")
write.table(goVAV3_cc, "go_cc (ICVAV3).txt")
write.table(goVAV3_mf, "go_mf (ICVAV3).txt")

#That's all. I hope you'll find it usefull.
