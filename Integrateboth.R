require("Seurat")

set.seed(3921)

if (!file.exists("CR_harmony_integrated.rds")) {
seurfiles <- c("Flush_Seurat.rds", "Caudate_Seurat.rds")

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x[[1]]}))
obj_list <- list()
common_genes <- c();
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(seurfiles[i]);
	obj@meta.data$cell_barcode <- colnames(obj);
	if (length(common_genes) == 0) {
		common_genes <- rownames(obj);
	} else {
		common_genes <- intersect(common_genes, rownames(obj));
	}
	obj@meta.data$donor <- rep(n, ncol(obj));
	obj_list[[n]] <- obj
}

common_genes <- sort(common_genes);
common_genes <- common_genes[!grepl("^MT-", common_genes)]
for (i in 1:length(obj_list)) {
	obj_list[[i]] <- obj_list[[i]][match(common_genes, rownames(obj_list[[i]])),]
}

# Harmony Workflow

merged_obj <- merge(obj_list[[1]], obj_list[[2]], project="Caudate_vs_Flush")


merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE) 
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableGenes(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)

png("CF_merged_not_integrated.png", width=6, height =6, units="in", res=200)
DimPlot(merged_obj, reduction="tsne", group.by="donor", pt.size=0.1)
dev.off();
png("CF_merged_not_integrated_umap.png", width=6, height =6, units="in", res=200)
DimPlot(merged_obj, reduction="umap", group.by="donor", pt.size=0.1)
dev.off();

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, "donor", plot_convergence = TRUE)
DimPlot(object = merged_obj, reduction = "harmony", pt.size = .1, group.by = "donor")

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)

png("CF_harmony_integrated_umap.png", width=6, height =6, units="in", res=150)
DimPlot(merged_obj, reduction = "umap", group.by = "donor", pt.size = .1)
dev.off();
png("CF_harmony_integrated_tsne.png", width=6, height =6, units="in", res=150)
DimPlot(merged_obj, reduction = "tsne", group.by = "donor", pt.size = .1)
dev.off();
png("CF_harmony_integrated_type.png", width=7.5, height =6, units="in", res=150)
DimPlot(merged_obj, reduction = "umap", group.by = "scmap_cell_anno", pt.size = .1)
dev.off();

# hand tuned to get separate inflammatory and non-inflammatory macrophages.
merged_obj <- FindNeighbors(merged_obj, dims = 1:20, k.param=30, reduction="harmony")
merged_obj <- FindClusters(merged_obj, resolution=1, k.param=30)

png("CF_hires_clusters_on_merged.png", width=5, height=5, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", pt.size = .1)
dev.off()
merged_obj$Tuned_clusters <- merged_obj$seurat_clusters
saveRDS(merged_obj, "CR_harmony_integrated.rds");
} else {
merged_obj <- readRDS("CR_harmony_integrated.rds")
} 


if (!file.exists("Merged_CR_dataobj.rds")) {
### Plot Frequencies of different cell types ###
merge_categories<- function(anno) {
	lvls <- levels(anno)
	anno <- as.character(anno)
	anno[grepl("Bcell", anno)] <- "Bcell"
	anno[grepl("gdTcell", anno)] <- "gdTcell"
	anno[grepl("LSEC", anno)] <- "LSEC"
	anno[grepl("Hep", anno)] <- "Hep"
	anno[grepl("Non-inf", anno)] <- "Non-Mac"
	anno[grepl("Infla", anno)] <- "Inf-Mac"
	anno[grepl("NK-", anno)] <- "NK-like"
	anno[grepl("Eryth", anno)] <- "Eryth"
	return(factor(anno, levels=c(lvls,"Bcell", "gdTcell", "LSEC", "Hep", "Inf-Mac", "Non-Mac","NK-like", "Eryth")))
}



caudate_freq <- table(merge_categories(obj_list$Caudate@meta.data$scmap_cell_anno))

flush_freq <- table(merge_categories(obj_list$Flush@meta.data$scmap_cell_anno))

type_freq_table <- cbind(caudate_freq/sum(caudate_freq)*100, flush_freq/sum(flush_freq)*100)
type_freq_table <- type_freq_table[apply(type_freq_table,1,max)>1,]
colnames(type_freq_table) <- c("Caudate", "Flush")
type_freq_table <- type_freq_table[order(apply(type_freq_table,1,max)),]

# Significance:
# Start from biggest difference to smallest
# if Sig remove all those cells from both datasets before continuing
# -> avoid sig because of proportionality b/c of other group that is sig

diff <- abs(type_freq_table[,1]-type_freq_table[,2])
p <- list()
for (type in names(diff)[order(diff, decreasing=T)]) {
	a = caudate_freq[type]
	b = sum(caudate_freq[names(caudate_freq) != type])
	d = flush_freq[type]
	e = sum(flush_freq[names(flush_freq) != type])
	
	out <- fisher.test(cbind(c(a,b),c(d,e)))
	p[[type]] <- out$p.value
	if (p[[type]] < 0.05/length(diff)){
		caudate_freq <- caudate_freq[names(caudate_freq) != type]
		flush_freq <- flush_freq[names(flush_freq) != type]
	}

}


png("CF_Cell_type_freqs.png", width=7, height=6, units="in", res=300)
par(mar=c(3,6,1,1))
loc <- barplot(t(type_freq_table), beside=T, las=1, horiz=T, xlab="% of cells", col=c("firebrick", "grey85"))
mtext("% of cells", 1,2)
sig <- names(p)[p< 0.01]
colnames(loc) <- rownames(type_freq_table)
xes <- colMeans(loc)
yes <- apply(type_freq_table, 1, max)
p_vals <- signif(unlist(p)[match(names(xes), names(p))], digits=1)
labs <- as.character(p_vals)
labs[p_vals > 0.01] <- ""
text(yes, xes, labs, pos=4)
legend("bottomright", c("Flush", "Caudate"), bty="n", fill=c("grey85", "firebrick"))
dev.off()


#merge clusters
merged_obj$orig.cluster <- paste(merged_obj$orig.ident,merged_obj$Use_clusters)
tab <- table(merged_obj$orig.cluster, merged_obj$scmap_cell_anno)
name_assign <- apply(tab, 1, function(x){colnames(tab)[which(x==max(x))]})
tab2 <- table(merged_obj$orig.cluster, merged_obj$seurat_clusters)
require("gplots")
png("CF_HarmonyC_vs_DatasetC.png", width=6, height=6, units="in", res=100)
heatmap.2(t(t(tab2)/colSums(tab2)), trace="none", scale="none", col=colorRampPalette(c("white", "black"))(20))
dev.off()


# compare dataset clusters with merged clusters.
heatmap.2(t(t(tab2)/colSums(tab2)), trace="none", scale="none", col=colorRampPalette(c("white", "black"))(20))
c_v_c <- AverageExpression(merged_obj)
heatmap.2(cor(c_v_c$RNA, method="spearman"))


# Manually merge similar clusters
merged_obj$Manual_clusters <- merged_obj$Tuned_clusters
merged_obj$Manual_clusters[merged_obj$Tuned_clusters %in% c(5,6,8)] <- 8
merged_obj$Manual_clusters[merged_obj$Tuned_clusters %in% c(9,0)] <- 9
merged_obj$Manual_clusters[merged_obj$Tuned_clusters %in% c(4,12)] <- 12

# Label merged clusters
tab <- table(merged_obj$Manual_clusters, merged_obj$scmap_cell_anno)
name_assign <- apply(tab, 1, function(x){if(sum(x)==0){return(NA)}; colnames(tab)[which(x==max(x))]})

i = 2;
while(sum(duplicated(name_assign)) > 0) {
	fix <- name_assign[duplicated(name_assign)]
	if (i > 2) {
		fix <- substr(fix, 1, nchar(fix)-1)
	} 
	fix <- paste(fix, i, sep="")
	name_assign[duplicated(name_assign)] <- fix
	i <- i+1;
}

merged_obj$Manual_clusters <- name_assign[merged_obj$Manual_clusters]

DimPlot(merged_obj, reduction = "umap", pt.size = .1, group.by="Manual_clusters")
merged_obj$seurat_clusters <- merged_obj$Manual_clusters
saveRDS(merged_obj, "Merged_CR_dataobj.rds")
} else {
merged_obj <- readRDS("Merged_CR_dataobj.rds")
}

## Differential Expression ##
if (!file.exists("CF_DE_output.rds")) {
de_obj <- merged_obj;
de_obj$seurat_clusters <- factor(merged_obj@meta.data$Manual_clusters)
de_obj@meta.data <- de_obj@meta.data[,c("orig.ident", "seurat_clusters")]

out <- list()
set.seed(38917)
for (group in 1:length(levels(de_obj$seurat_clusters))) {
	type_makers <- FindMarkers(de_obj, ident.1=group, logfc.threshold=-Inf)
	de_genes <- FindMarkers(de_obj, ident.1 = "Flush", ident.2 = "Caudate", group.by="orig.ident", subset.ident = group, logfc.threshold=-Inf)
	# Add mean expression in each group to double check DE
	out[[group]] <- list(type=type_makers, de=de_genes);
	file_name=paste("CF", levels(de_obj$seurat_clusters)[group],"quickWilcoxDE.txt", sep="_")
	write.table(type_makers, file=file_name)
}
names(out) <- levels(de_obj$seurat_clusters)
saveRDS(out, "CF_DE_output.rds")
} else {
out <- readRDS("CF_DE_output.rds")
}

# Turn this into a matrix of effect sizes across cell-types
# This will facilitate removing effects that are consistent across cell-types 
# which are more likely to be due to batch.
#
l2fc_tab <- vector()
diff_frac_tab <- vector()
fdr_tab <- vector()
for(thing in names(out)) {
	l2fc <- out[[thing]]$de[,2]
	names(l2fc) <- rownames(out[[thing]]$de)
	diff_frac <- out[[thing]]$de[,3]-out[[thing]]$de[,4]
	names(diff_frac) <- rownames(out[[thing]]$de)
	fdr <- out[[thing]]$de[,5]
	names(fdr) <- rownames(out[[thing]]$de)
	if (length(l2fc_tab) == 0) {
		l2fc_tab <- matrix(l2fc, ncol=1)
		rownames(l2fc_tab) <- names(l2fc)
		colnames(l2fc_tab) <- thing
	} else {
		tmp_rnames <- unique(sort(c(rownames(l2fc_tab), names(l2fc))));
		tmp_cnames <- colnames(l2fc_tab)
		l2fc_tab <- cbind(l2fc_tab[match(tmp_rnames, rownames(l2fc_tab)),], l2fc[match(tmp_rnames, names(l2fc))])
		colnames(l2fc_tab) <- c(tmp_cnames, thing) #probs need to fix that.
		l2fc_tab[is.na(l2fc_tab)] <- 0;
		rownames(l2fc_tab) <- tmp_rnames
	}
	if (length(diff_frac_tab) == 0) {
		diff_frac_tab <- matrix(diff_frac, ncol=1)
		rownames(diff_frac_tab) <- names(diff_frac)
		colnames(diff_frac_tab) <- thing
	} else {
		tmp_rnames <- unique(sort(c(rownames(diff_frac_tab), names(diff_frac))));
		tmp_cnames <- colnames(diff_frac_tab)
		diff_frac_tab <- cbind(diff_frac_tab[match(tmp_rnames, rownames(diff_frac_tab)),], diff_frac[match(tmp_rnames, names(diff_frac))])
		colnames(diff_frac_tab) <- c(tmp_cnames, thing) #probs need to fix that.
		diff_frac_tab[is.na(diff_frac_tab)] <- 0;
		rownames(diff_frac_tab) <- tmp_rnames
	}
	if (length(fdr_tab) == 0) {
		fdr_tab <- matrix(fdr, ncol=1)
		rownames(fdr_tab) <- names(fdr)
		colnames(fdr_tab) <- thing
	} else {
		tmp_rnames <- unique(sort(c(rownames(fdr_tab), names(fdr))));
		tmp_cnames <- colnames(fdr_tab)
		fdr_tab <- cbind(fdr_tab[match(tmp_rnames, rownames(fdr_tab)),], fdr[match(tmp_rnames, names(fdr))])
		colnames(fdr_tab) <- c(tmp_cnames, thing) #probs need to fix that.
		fdr_tab[is.na(fdr_tab)] <- 0;
		rownames(fdr_tab) <- tmp_rnames
	}
}

# I suspect bias in diff_frac due to different numbers of cells & thus reads/cell across the two.
#disagree <- sign(diff_frac_tab) * sign(l2fc_tab)
#diff_frac_tab[disagree < 0] <- 0
#diff_l2fc_tab[disagree < 0] <- 0

sig_l2fc_tab <- l2fc_tab
sig_l2fc_tab[fdr_tab > 0.05] <- 0

consistent_up <- apply(sig_l2fc_tab, 1, function(x) {sum(x > 0)})
consistent_dw <- apply(sig_l2fc_tab, 1, function(x) {sum(x < 0)})
consistent_up <- consistent_up[consistent_up>0]
consistent_dw <- consistent_dw[consistent_dw>0]

png("Flush_vs_Caudate_DE.png", width=5, height=5, units="in", res=100)
#Create Data
par(mfrow=c(2,1))
 
#Make the plot
par(mar=c(0,5,3,3))
hist(consistent_up , main="" , xlim=c(0,16), ylab="Genes Increases", xlab="", ylim=c(0,2000) , xaxt="n", las=1 , col="tomato3", breaks=14)
par(mar=c(5,5,0,3))
hist(consistent_dw , main="" , xlim=c(0,16), ylab="Genes Decreases", xlab="Number of Celltypes", ylim=c(2000,0) , las=1 , col="slateblue1"  , breaks=14)
dev.off()

###
require(GSEABase)
require(GSVA)
library(parallel)
pathways <- getGmt("~/Tallulah/GeneSets/Human_GO_AllPathways_no_GO_iea_February_01_2020_symbol.gmt")
test <- gsva(l2fc_tab, pathways, parallel.sz=1)


score <- apply(l2fc_tab,1,function(x){sum(x>0)})-apply(l2fc_tab,1,function(x){sum(x<0)})
heatmap.2(l2fc_tab[abs(score) < 3,], trace="none", symbreaks=TRUE, col=colorRampPalette(c("magenta","black","yellow"))(20) )

score <- apply(l2fc_tab,1,function(x){sum(x>0)})-apply(l2fc_tab,1,function(x){sum(x<0)})
heatmap.2(l2fc_tab[abs(score) < 5 & apply(abs(l2fc_tab),1,max) > 1,], trace="none", symbreaks=TRUE, col=colorRampPalette(c("magenta","black","yellow"))(20) )


tmp <- prcomp(diff_frac_tab)

# GSVA?

test <- gsva(l2fc_tab, pathways, parallel.sz=1)

apply(test,2, function(x){rownames(test)[x >= x[order(x, decreasing=T)[10]]]})
#This works but I don't know how to interpret it...



# Filter for cell-type specific effects
# heatmap, pathway analysis.

# Seurat's default colour scheme
#require(scales)
#color_palette <- hue_pal(n_groups);
