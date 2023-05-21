args <- commandArgs(trailingOnly = TRUE)
gmt<- args[1]
expr<- args[2]
wd<- args[3]
suppressMessages(library(Seurat))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(corrplot))

#Read gene expression profile and hypoxic gene sets data
setwd(wd)
geneSets <- getGmt(gmt)
if(grepl("RData$",expr,ignore.case = TRUE)){
  data_matrix=get(load(expr))
} else if(grepl("txt$|csv$",expr,ignore.case = TRUE)){
  data_matrix=read.table(expr,header = T,sep = "\t",row.names = 1)
}
scRNA_obj<-CreateSeuratObject(counts=data_matrix)
scRNA_matrix = GetAssayData(scRNA_obj, assay="RNA", slot = "data")
save(scRNA_matrix,file="scRNA_matrix.RData")
#Calculating the hypoxia activity score for each gene set in each cell
hypoxia_score<- gsva(as.matrix(scRNA_matrix), geneSets, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)
save(hypoxia_score, file = 'hypoxia_score.RData')

#The correlation between cell activity scores in pairwise gene sets
res <- cor(t(hypoxia_score), method = 'pearson')
pdf('GSVA_corr.pdf')
corrplot(res, method = "squar", shade.col = NA, tl.col ="black", tl.srt = 60, order = "AOE",addCoef.col="grey")
dev.off()