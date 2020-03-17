require(dplyr)
require(Seurat)
require(Matrix)
set.seed(12819)
ngenes <- 200
mt_perc=40
projname="Caudate"
npc <- 15
coarse_lvl <- "knn_70_res_1.3"
fine_lvl <- "knn_40_res_1.5"
use <- "Fine_clusters"
knn = 40
file <- "Caudate/filtered_gene_bc_matrices/GRCh38/"
#projname="Flush"
#ngenes <- 750
#npc <- 17
#coarse_lvl <- "knn_80_res_0.7"
#fine_lvl <- "knn_40_res_1.3"
#use <- "Coarse_clusters"
#knn=80
#file <- "Flush/filtered_gene_bc_matrices/GRCh38/"

# Read the data
mydata <- Read10X(data.dir =file)
myseur <- CreateSeuratObject(counts = mydata, project = projname, min.cells=10, min.features=100)

# Mitochondrial filter & QC filters
myseur[["percent.mt"]] <- PercentageFeatureSet(myseur, pattern = "^MT-")
png(paste(projname, "_MTqc.png", sep=""), width=6, height=6, units="in", res=50) 
plot(myseur@meta.data$nFeature_RNA, myseur@meta.data$percent.mt)
abline(v=ngenes, col="red")
abline(h=mt_perc, col="red")
dev.off()

myseur <- subset(myseur, subset = nFeature_RNA > ngenes  & percent.mt < mt_perc)

# PCA
myseur <- NormalizeData(myseur)
myseur <- ScaleData(myseur);
myseur <- FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)
myseur <- RunPCA(myseur, features = VariableFeatures(object = myseur))
ElbowPlot(myseur)

npcs = npc
# Clustering
myseur <- FindNeighbors(myseur, dims = 1:npcs)
res <- seq(from=0.3, to=2, by=0.2)
nkNN <- seq(from=30, to=90, by=10)
for(res_param in res) {
for(nkNN_param in nkNN){
	myseur <- FindNeighbors(myseur, dims = 1:npcs, k.param=nkNN_param)
	myseur <- FindClusters(myseur, resolution = res_param, k.param=nkNN_param)
	name <- paste("knn_",nkNN_param,"_res_", res_param, sep="");
	myseur@meta.data[[name]] <- myseur@meta.data$seurat_clusters;
}}
# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <-myseur@meta.data[, grepl("^knn_", names(myseur@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)
clust_similr1 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="nmi")}, by_rows=FALSE)
clust_similr2 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="adjusted.rand")}, by_rows=FALSE)
# Find robust exemplar clustering(s)
require("apcluster")
require("gplots")
set.seed(18371)
res1 <- apcluster(-1*as.matrix(clust_dists), p=-2.5)
myseur$Coarse_clusters <- myseur@meta.data[,coarse_lvl]
myseur$Fine_clusters <- myseur@meta.data[,fine_lvl]
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==fine_lvl, colnames(clust_table)==fine_lvl] <- "F"
lab[colnames(clust_table)==coarse_lvl, colnames(clust_table)==coarse_lvl] <- "C"

png(paste(projname, "_clustOptim.png", sep=""), width=6, height=6, units="in", res=50) 
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()
myseur$Use_clusters <- myseur@meta.data[,use]


nkNN <- knn
myseur <- RunTSNE(myseur, dims=1:npcs)
myseur <- RunUMAP(myseur, dims=1:npcs, parallel=FALSE, n.neighbour=nkNN)
png(paste(projname, "_umap.png", sep=""), width=6, height=6, units="in", res=300) 
DimPlot(myseur, reduction = "umap", group.by="Use_clusters")
dev.off()
png(paste(projname, "_tsne.png", sep=""), width=6, height=6, units="in", res=300) 
DimPlot(myseur, reduction = "tsne", group.by="Use_clusters")
dev.off()


set.seed(8410)
source("~/Tallulah/AutoAnnotation/scripts/AutoAnno_LiverMap1.0/Setup_autoannotation.R")
myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);
norm <- myseur@assays$RNA@data
clus_lab <- myseur@meta.data$Use_clusters
res <- Use_markers_for_anno(norm, clus_lab)

myseur@meta.data$marker_anno <- res$cell_assign

saveRDS(myseur, paste(projname, "_Seurat.rds", sep=""))

png(paste(projname, "_markAnno_umap.png", sep=""), width=7.5, height=6, units="in", res=300) 
DimPlot(myseur, reduction = "umap", group.by="marker_anno")
dev.off()

png(paste(projname, "_scmapAnno2_umap.png", sep=""), width=7.5, height=6, units="in", res=300) 
DimPlot(myseur, reduction = "umap", group.by="scmap_cell_anno")
dev.off()

png(paste(projname, "_scmapAnno2_tsne.png", sep=""), width=7.5, height=6, units="in", res=300) 
DimPlot(myseur, reduction = "tsne", group.by="scmap_cell_anno")
dev.off()

cluster_entropy <- function(x) {
        freqs <- table(x)/length(x);
        Shannon <- -sum(freqs * log2(freqs))
	# Equitability
	Shannon <- Shannon/log2(length(freqs))
	if (length(freqs)==1){return(0)}
        return(Shannon);
}
clustering_entropy <- function(cell_labs, anno_labs) {
	cluster_entropies <- apply(unique(cell_labs), function(c) {cluster_entropy(anno_labs[cell_labs==c])})
	nclusters <- length(unique(cell_labs));
	return(c(nclusters, mean(cluster_entropies)))
}

# Find clusters based on minimal entropy of autoannotations

clust_table <-myseur@meta.data[, grepl("^knn_", names(myseur@meta.data))]

res <- apply(clust_table, 2, clustering_entropy, myseur@meta.data$scmap_cell_anno)


