args <- commandArgs(trailingOnly = TRUE)
wd<- args[1]
setwd(wd)

load("score_cluster.RData")
load("scRNA_matrix.RData")
#Count the frequency of each cell classified as 1 and 0
score_cluster_count <- apply(score_cluster, 1, function(x) {
  return(c(length(which(x == 1)), length(which(x == 0))))
})
rownames(score_cluster_count) <- c('hypoxia', 'normoxia')
max_num<-max(score_cluster_count)

#Identify high confidence cells by classification frequency
Hypoxia <- score_cluster_count['hypoxia',]
Hypoxia <- names(Hypoxia[which(Hypoxia == max_num)])
Normoxia <- score_cluster_count['normoxia',]
Normoxia <- names(Normoxia[which(Normoxia == max_num)])
cell_number<-min(length(Hypoxia),length(Normoxia))
high_confi_cellID <- list(Hypoxia = Hypoxia, Normoxia = Normoxia)
save(high_confi_cellID, file = 'high_confi_cellID.RData')

train_cells <- unname(unlist(high_confi_cellID ))
all_cells<- colnames(scRNA_matrix)
unclassified_cells <- all_cells[-which(all_cells%in% train_cells)]

#The top500 differential expression of genes as input features of the classifier
expr<-scRNA_matrix[which(rowSums(as.matrix(scRNA_matrix) > 0)>dim(scRNA_matrix)[2]*0.01),]
expr_train<-expr[,train_cells]
gene_pvalue=apply(expr_train,1,function(x) wilcox.test(x[Hypoxia],x[Normoxia])$p.value)
gene_diff_sort <- sort(gene_pvalue)
gene_feature <- names(gene_diff_sort)[1:500]
expr_train<-expr_train[gene_feature,]
expr_train<-as.matrix(expr_train)
expr_unclassified<-scRNA_matrix[gene_feature,unclassified_cells]
expr_unclassified<-as.matrix(expr_unclassified)
#Expression profile of high confidence cells(training set)
write.csv(t(expr_train), file = 'expr_highconfi.csv',quote=F)
#Expression profile of other cells (cells to be classified)
write.csv(t(expr_unclassified), file = 'expr_others.csv', quote = F)

#Training set label
label<- matrix(0, ncol = length(train_cells), nrow = 1)
colnames(label) <- train_cells
rownames(label) <- 'group'
label[1, high_confi_cellID$Hypoxia] <- 1
label[1, high_confi_cellID$Normoxic] <- 0
write.csv(t(label), file = 'label_highconfi.csv', row.names = train_cells, quote = F)


