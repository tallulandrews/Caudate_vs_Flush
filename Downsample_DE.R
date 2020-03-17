flush <- readRDS("Flush_Seurat.rds")
caudate <- readRDS("Caudate_Seurat.rds")
common_genes <- rownames(flush)[rownames(flush) %in% rownames(caudate)]
flush <- flush[match(common_genes, rownames(flush)),]
caudate <- caudate[match(common_genes, rownames(caudate)),]

require("Seurat")
require("Matrix")
set.seed(10192)
caudate_downsample <- SampleUMI(caudate@assays$RNA@counts, max.umi=1500)
flush_downsample <- SampleUMI(flush@assays$RNA@counts, max.umi=1500)


# Use autoannotations directly for DE
label_caudate <- caudate@meta.data$scmap_cell_anno
label_flush <- flush@meta.data$scmap_cell_anno


simplify_labels <- function(labs) {
	labs<- as.character(labs)
	labs[grep("Hep", labs)] <- "Hepatocyte"
	return(labs);
}

label_caudate <- simplify_labels(label_caudate);
label_flush <- simplify_labels(label_flush);

# DE - only if at least 30 cells in each group.
tab_caud <- table(label_caudate);
tab_caud <- tab_caud[tab_caud > 30];

tab_flush <- table(label_flush);
tab_flush <- tab_flush[tab_flush > 30];

to_de <- names(tab_flush)[names(tab_flush) %in% names(tab_caud)]

mat <- cbind(caudate_downsample, flush_downsample);
mat <- t(t(mat)/colSums(mat)*1500)
labs <- c(label_caudate, label_flush);
conds <- rep(c("caudate", "flush"), times=c(ncol(caudate_downsample), ncol(flush_downsample)));

my_wilcox.test <- function(mat, labels, condition, group) {
	set <- mat[,labels==group];
	cond <- condition[labels==group]
	ps <- apply(set, 1, function(x){wilcox.test(x[cond=="flush"], x[cond=="caudate"])$p.value})	 
	mean_c <- Matrix::rowMeans(set[,cond=="caudate"])
	detect_c <- Matrix::rowMeans(set[,cond=="caudate"]>0)
	mean_f <- Matrix::rowMeans(set[,cond=="flush"])
	detect_f <- Matrix::rowMeans(set[,cond=="flush"]>0)
	return(data.frame(mean_flush=mean_f, mean_caudate=mean_c, detect_flush=detect_f, detect_caudate=detect_c, pval=ps, qval=p.adjust(ps, method="fdr")))

}
all_out <- list()
for(ctype in to_de) {
	out <- my_wilcox.test(mat, labs, conds, ctype);
	all_out[[ctype]] <- out;
}

saveRDS(all_out, "Downsampled_Wilcox_Flush_vs_Caudate.rds");
