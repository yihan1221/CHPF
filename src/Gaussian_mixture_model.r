suppressMessages(library(mclust))
suppressMessages(library(pheatmap))
args <- commandArgs(trailingOnly = TRUE)
wd<- args[1]

setwd(wd)
load("hypoxia_score.RData")
geneset_name <- rownames(hypoxia_score)

#GMM clustering 
Cluster <- list()
for (gsname in geneset_name) {
	fit_GMM <- Mclust(hypoxia_score[gsname,], G = 2)
	Cluster <- c(Cluster, list(fit_GMM["classification"]))
}
names(Cluster) <- geneset_name

# Compare the mean activity scores in cluster1 and cluster2, 1 for high score(hypoxia) and 0 for low score (normoxia)
score_cluster <- list()
for (gsname in geneset_name) {
	if (mean(hypoxia_score[gsname, which(unlist(Cluster[gsname]) == '1')]) > mean(hypoxia_score[gsname, which(unlist(Cluster[gsname]) == '2')])) {
		tmp_cluster <- 2 - unname(unlist(Cluster[gsname]))
	} else {
		tmp_cluster <- unname(unlist(Cluster[gsname])) - 1
	}
	score_cluster <- c(score_cluster, list(tmp_cluster))
}
score_cluster <- as.data.frame(score_cluster)
colnames(score_cluster) <- geneset_name
rownames(score_cluster) <- colnames(hypoxia_score)
save(score_cluster, file = 'score_cluster.RData')

pdf("GMM_classification.pdf")
pheatmap(t(score_cluster),clustering_method = "ward.D", show_colnames = F,show_rownames =T,col=c("#FDE9DC","#C6403D"))
dev.off()

