tab <- read.table("Flush_vs_Caudate_signifDE_table.txt", header=T)

# dissociation heatmap
dissociation <- read.table("Brink_Dissociation_marks.txt", header=T)

#dissociation score
types <- unique(tab$Ctype)
is.diss <- tab$Gene %in% dissociation$Hgene
is.Mt <- grepl("MT-", tab$Gene);
diss.score <- unlist(lapply(types, function(a) {mean(tab$mean_flush[is.diss & tab$Ctype == a]) - 
				mean(tab$mean_caudate[is.diss & tab$Ctype == a])}))
names(diss.score) <- as.character(types)

diss.n <- unlist(lapply(types, function(a) {
		length(tab$mean_flush[is.diss & tab$Ctype == a & tab$dir == 1])}))
names(diss.n) <- as.character(types)



Mt.score <- unlist(lapply(types, function(a) {mean(tab$mean_flush[is.Mt & tab$Ctype == a]) - 
				mean(tab$mean_caudate[is.Mt & tab$Ctype == a])}))
names(Mt.score) <- as.character(types)
Mt.n <- unlist(lapply(types, function(a) {
		length(tab$mean_flush[is.Mt & tab$Ctype == a & tab$dir == -1])}))
names(Mt.n) <- as.character(types)



png("Flush_vs_Caudate_MT_Diss_scores.png", width=4, height=6, units="in", res=300)
par(mfrow=c(2,1))
par(mar=c(3,6,1,1))
pos <- barplot(Mt.score, horiz=T, las=1, col="magenta", xlim=c(-23,0))
mtext("MT expression in Flush vs Caudate", 1, 2)
text(Mt.score, pos, Mt.n, pos=2)
pos2 = barplot(diss.score, horiz=T, las=1, col="goldenrod1", xlim=c(0,1))
mtext("Dissociation in Flush vs Caudate", 1, 2)
text(diss.score, pos2, diss.n, pos=4)
dev.off()

# Remove MT and ribosomal
tab <- tab[!grepl("MT-", tab$Gene),]
tab <- tab[!grepl("^RP", tab$Gene),]

tab <- tab[order(tab$qval),]

# each one
require("gprofiler2")

set.seed(101)

up_res <- list();
down_res <- list();
for (type in unique(tab$Ctype)) {
	genes <- tab[tab$dir > 0 & tab$Ctype == type,"Gene"]
	res <- gprofiler2::gost(genes, ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result

	# filter results
	res <- res[res$intersection_size > res$query_size/10 & res$intersection_size > 1 & res$term_size <= res$query_size*10,]
	res <- res[order(res$p_value),]

	up_res[[type]] <- head(res, 25)

	genes <- tab[tab$dir < 0 & tab$Ctype == type,"Gene"]
	res <- gprofiler2::gost(genes, ordered_query=T, correction_method="fdr", sources=c("GP:BP", "KEGG", "REAC"))$result

	# filter results
	res <- res[res$intersection_size > res$query_size/10 & res$intersection_size > 1 & res$term_size <= res$query_size*10,]
	res <- res[order(res$p_value),]

	down_res[[type]] <- head(res, 25)

}

create_table <- function(gprofiler2_output) {
	all_tab <- c();
	for (type in names(gprofiler2_output)) {
		this_tab <- gprofiler2_output[[type]][,c("term_name", "intersection_size", "p_value")]
		this_tab$type <- rep(type, nrow(this_tab));
		all_tab <- rbind(this_tab, all_tab);
	}
	return(all_tab)
}

up_tab <- create_table(up_res)
up_tab$condition <- rep("Flush", nrow(up_tab));
dw_tab <- create_table(down_res)
dw_tab$condition <- rep("Caudate", nrow(dw_tab));

all <- rbind(up_tab, dw_tab)
write.table(all, "top25_enrichments_Flush_vs_Caudate.csv", sep=",", row.names=F, col.names=T)

toplot <- read.table("top25_enrichments_Flush_vs_Caudate_trimmed.csv", sep=",", header=T)
toplot <- toplot[toplot$p_value >0 & !is.na(toplot$p_value),]
require("ggplot2")
toplot$score <- -log10(toplot$p_value)

# Both plotted
my_plot <- toplot
my_plot <- my_plot[my_plot$term_name != "Metabolism",]
my_plot[my_plot$condition == "Caudate","score"] <- -1*my_plot[my_plot$condition == "Caudate","score"]
my_plot_mat <- matrix(0, nrow=length(unique(my_plot$term_name)), ncol=length(unique(my_plot$type)))
rownames(my_plot_mat) <- unique(my_plot$term_name)
colnames(my_plot_mat) <- unique(my_plot$type)
for (i in 1:nrow(my_plot)) {
	my_plot_mat[as.character(my_plot[i,1]), as.character(my_plot[i,4])] <- my_plot[i,7]
}

png("top_enrichments_Flush_vs_Caudate.png", width=8, height=15, units="in", res=300)
require(gplots)
heatmap.2(my_plot_mat, dendrogram="none", trace="none", symbreaks=T, 
		col=colorRampPalette(c("magenta", "black", "yellow"))(20),
		cexRow=1, cexCol=1, mar=c(15,15), key.title="", key.xlab="-log(p)")
dev.off()