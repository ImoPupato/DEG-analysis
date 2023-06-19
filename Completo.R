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

# - Acces to the clinical information
clinica<-as.data.frame(SKCM.clinical)
clinica$bcr_patient_barcode<-toupper(clinica$patient.bcr_patient_barcode) 
id<-as.data.frame(clinica$bcr_patient_barcode)
x<-x[,c(1,9,11,13,16,28,30,31,438,439,448,449,450,770,775,778,779,780,785,917,935,939,940,941,944,945,846,949,956,964,1060,1867)]
write.table(x,"clinical information.txt")

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
id2<-as.data.frame(colnames(rnaseq))
colnames(id2)<-c("bcr_patient_barcode")
id2<-id%>%mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
colnames(rnaseq)<-id2$bcr_patient_barcode

# - Select according to the clinical informataion available
rnaseq<-select(rnaseq, by=c(x$id))
CPM<-rnaseq[,-1] # eliminate the ENSEMBL column and save the CPM for later
id3<-as.data.frame(colnmaes(rnaseq))

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

# - Survival analysis and cutpoints
VAVs_SKCM.surv<-survivalTCGA(SKCM.clinical)
colnames(VAVs_SKCM.surv)<-c("times","ID","status")
CPM<-cpm(CPM)
CPM<-cbind(rnaseq$SYMBOL,CPM)
phenotype<-t(rbind(CPM[CPM$SYMBOL=="VAV1",],CPM[CPM$SYMBOL=="VAV2",],CPM[CPM$SYMBOL=="VAV3",]))
phenotype<-cbind(id3,phenotype)
colnames(phenotype)<-c("ID","VAV1","VAV2,"VAV3)
VAVs_SKCM.surv_rnaseq<- VAVs_SKCM.surv %>% left_join(CPM,by = "ID")
VAVs_SKCM.surv_rnaseq.cut<-surv_cutpoint(
  VAVs_SKCM.surv_rnaseq,
  time = "times",
  event = "status",
  variables = c("VAV1","VAV2", "VAV3")
)
summary(VAVs_SKCM.surv_rnaseq.cut)

library(survival)
fit <- survfit(Surv(times, status) ~ Vav*,
               data = VAVs_SKCM.surv_rnaseq.cat)
ggsurvplot(
   fit,                     
   risk.table = TRUE,       
   pval = TRUE,             
   conf.int = FALSE,                  
   xlim = c(0,5500),      
   break.time.by = 1000,    
   ggtheme = theme_bw(), 
   risk.table.y.text.col = T, 
  font.legend=8, 
  legend.title = "Expression",
  risk.table.y.text = FALSE,                          
)


# Categorize according to Vav1,2 y 3 expression
phenotype$PhenoV1<-ifelse(phenotype$VAV1<7.717641,"Low","High")
phenotype$PhenoV2<-ifelse(phenotype$VAV2<57.475259,"Low","High")
phenotype$PhenoV3<-ifelse(phenotype$VAV3<9.791567,"Low","High")

#The DEG itself (when you see a *, means you decide what Vav, 1-3)
y <- DGEList(counts=CPM[,-1],group=c(phenotype$PhenoV*) # Here you'll be creating the DGElist object.
keep <- filterByExpr(y) # In this step, the edgeR algorithm will be keeping the rows that are 'worthwhile' keeping (doesn't have too many zeros).
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~c(phenotype$PhenoV*) # Your design comes from the condition previously set. Here you can add another interaction using ~ condition A + condition B
y0 <- estimateDisp(y,design) # In this step I change the name because in some analysis the estimateDisp function isn't needed.
et <- exactTest(y0) # Here the pseudocounts matrix is made.
is.de <- decideTestsDGE(et) # We examin the total nmuber if DEG.
summary(is.de)
df_etV*<-as.data.frame(topTags(et,n=Inf))
dim(df_etV*) # Just for control
head(df_etV*)
df_etV*$Expression = ifelse(df_etV*$FDR < 0.01 & abs(df_etV*$logFC) >= 1, 
                               ifelse(df_etV*> 1 ,'Up','Down'),
                               'Stable') # I add a column with the expression for furhter use.
df_etV*$Genes = rownames(df_etV*) # Prepare for saving the data.
write.table(df_etV*, "Exact test (Vav*).txt") 

#Performing the VULCANO PLOT to visualize de Genes and its expression
plot <- ggplot(data = df_etV*, 
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

#- Venn Diagram (Shown for comon, then repeat for Up & Down)
library("ggVennDiagram")
gene_list <- list(VAV1 = c(subset(df_etV1$Gene,df_etV1$Expression!="Stable")),
                  VAV2 = c(subset(df_etV2$Gene,df_etV2$Expression!="Stable")),
                  VAV3 = c(subset(df_etV3$Gene,df_etV3$Expression!="Stable")))
p1 <- ggVennDiagram(gene_list, label_alpha = 0,
                    category.names = c("VAV1","VAV2","VAV3"))
p1
comungenes<-intersect(VAV1,intersect(VAV2,VAV3))#repeat fo Up & Down
                       
# To perform the ORA 
library(ReactomePA)
library(enrichplot)

## Creating the new data
hs <- org.Hs.eg.db
Exact_test_Vav*<-subset(df_et,df_et$Expression!="Stable") # I'll be only using the up and down DEGs.
my.symbols <- c(Exact_test_Vav*$Genes) # The string with the genes names to map with ENTREZ
IDS<-AnnotationDbi::select(org.Hs.eg.db, keys=my.symbols, columns=c("SYMBOL", "ENTREZID"), keytype="SYMBOL")
names_ids<-IDS$ENTREZID
b<-na.omit(IDS) # I've never had na but, just to make sure.
                     
x <- enrichPathway(gene=names_ids,pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory=10, title="DotPlot Vav*",font.size=8) # You can also create a CNet Plot and Bar Plot

#That's all. I hope you'll find it usefull.
