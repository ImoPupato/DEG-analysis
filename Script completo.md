Setting your working directory:  
```R
setwd("C:/Mi unidad/Analisis/DEG")
```
# Libraries. 
**From the RTCGA package**  
_This libraries allows us to download and integrate the variety and volume of TCGA._
```R
library("TCGAbiolinks")
library("RTCGA")
library("RTCGA.clinical")
llibrary("SummarizedExperiment")
```
 _If you don't have them already installed you must do it using the BiocManager. For more information you can visit this [site](https://rtcga.github.io/RTCGA)_
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RTCGA")
```
_Then call the packages for the DEG analysis. I'm using edgeR_
```R
library("edgeR")
```
_Last but not least, the libraries for the manipulation of the data and other analysis_
```R
library("tidyverse") # contains tidyr, ggplot2 and dplyr, among others
library("survival") # for the Survival analysis
library("survminer") # for the cutpoints values
library("org.Hs.eg.db") # To acces to all the gene notations
library("ReactomePA") # for the ORA 
library("enrichplot") # for the enrichment analysis
library("ggVennDiagram") # for the Venns diagramas
```
# Access to the clinical information.
```R
clinica<-as.data.frame(SKCM.clinical)
clinica$bcr_patient_barcode<-toupper(clinica$patient.bcr_patient_barcode) # to change the format of the id
id<-as.data.frame(clinica$bcr_patient_barcode) # extracting the id to an object that will be used later
write.table(clinica,"clinical information.txt") # saving the data
```
# Query and download the RNAseq
```R
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
```
_Reading the skcm.rda file_
```R
rnaseq <- assay(skcm.exp)
rnaseq <- rnaseq[which(!duplicated(rownames(rnaseq))),]   # to eliminate the duplicated rows
write.table(rnaseq, "raw counts.txt") # To save the data set
```
# Manipulate an clean the data set
## - Mutate the ids 
```R
auxiliar_id<-as.data.frame(colnames(rnaseq))
colnames(auxiliar_id)<-c("bcr_patient_barcode")
auxiliar_id<-auxiliar_id%>%mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
colnames(rnaseq)<-auxiliar_id$bcr_patient_barcode
```
## - Select according to the clinical informataion available
```
rnaseq<-select(rnaseq, by=c(x$id)) # remember the id object generated above
CPM<-rnaseq[,-1] # eliminate the ENSEMBL column and save the CPM for later
auxiliar_id<-as.data.frame(colnmaes(rnaseq))
```
## - Add the ENTREZ, GENE SYMBOL or ENSEMBL gene ID  
_Depending on the package used, it may need differet gene notations_  
```
IDS<-select(org.Hs.eg.db, 
            keys = c(rnaseq$ENSEMBL),
            columns = c("ENSEMBL", "SYMBOL","ENTREZID"),
            keytype = "ENSEMBL")
IDS<-na.omit(IDS)
rnaseq<-merge(rnaseq,IDS,by="ENSEMBL")
```
# Survival analysis 
## - Prepare the data
```R
VAVs_SKCM.surv<-survivalTCGA(SKCM.clinical) # acces to survival status
colnames(VAVs_SKCM.surv)<-c("times","ID","status")
CPM<-cpm(CPM) # here you can use FPKM, CPM, TMM.
CPM<-cbind(rnaseq$SYMBOL,CPM)
phenotype<-t(rbind(CPM[CPM$SYMBOL=="VAV1",],CPM[CPM$SYMBOL=="VAV2",],CPM[CPM$SYMBOL=="VAV3",]))
phenotype<-cbind(auxiliar_id,phenotype)
colnames(phenotype)<-c("ID","VAV1","VAV2,"VAV3)
VAVs_SKCM.surv_rnaseq<- VAVs_SKCM.surv %>% left_join(CPM,by = "ID")
```
## - Select the cutpoints
```R
VAVs_SKCM.surv_rnaseq.cut<-surv_cutpoint(
  VAVs_SKCM.surv_rnaseq,
  time = "times",
  event = "status",
  variables = c("VAV1","VAV2", "VAV3")
)
summary(VAVs_SKCM.surv_rnaseq.cut)
```
## - Fit and graph the KM
```R
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
```
# Categorize according to Vav1,2 y 3 expression
```R
phenotype$PhenoV1<-ifelse(phenotype$VAV1<7.717641,"Low","High")
phenotype$PhenoV2<-ifelse(phenotype$VAV2<57.475259,"Low","High")
phenotype$PhenoV3<-ifelse(phenotype$VAV3<9.791567,"Low","High")
```
# The DEG itself  
_when you see a *, means you decide what Vav, 1-3_
```R
y <- DGEList(counts=rnaseq[,-1],group=c(phenotype$PhenoV*) # Here you'll be creating the DGElist object. You should check y[["samples"]]$group to make sure the groups hd been asigned. If not, you can use y[["samples"]]$group<-auxiliar2$phenotype$PhenoV* 
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
```
# VULCANO PLOT to visualize the Genes and its expression
```R
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
```
# Venn Diagram  
_Shown for common, then repeat for Up & Down_
```R
gene_list <- list(VAV1 = c(subset(df_etV1$Gene,df_etV1$Expression!="Stable")),
                  VAV2 = c(subset(df_etV2$Gene,df_etV2$Expression!="Stable")),
                  VAV3 = c(subset(df_etV3$Gene,df_etV3$Expression!="Stable")))
p1 <- ggVennDiagram(gene_list, label_alpha = 0,
                    category.names = c("VAV1","VAV2","VAV3"))
p1
comungenes<-intersect(VAV1,intersect(VAV2,VAV3))#repeat fo Up & Down
```                       
# ORA 
## Creating the new data
```R
Exact_test_Vav*<-subset(df_et,df_et$Expression!="Stable") # I'll be only using the up and down DEGs.
IDS<-AnnotationDbi::select(org.Hs.eg.db,
                           keys= c(Exact_test_Vav*$Genes), # The string with the genes names to map with ENTREZ
                           columns=c("SYMBOL", "ENTREZID"),
                           keytype="SYMBOL")
names_ids<-IDS$ENTREZID
b<-na.omit(IDS) # I've never had na but, just to make sure.
                     
x <- enrichPathway(gene=names_ids,pvalueCutoff=0.05, readable=T)
dotplot(x, showCategory=10, title="DotPlot Vav*",font.size=8) # You can also create a CNet Plot and Bar Plot
```
# That's all. I hope you'll find it usefull.
