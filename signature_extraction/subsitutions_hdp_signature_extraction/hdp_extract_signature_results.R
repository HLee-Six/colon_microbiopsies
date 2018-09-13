.libPaths("/nfs/users/nfs_h/hl11/R/R-3.4.0")
library("hdp")
library("RColorBrewer")

hdpouts <- list.files(pattern="RData")
for (hdpout in hdpouts) {
  load(hdpout)
}

outs <- ls(pattern="hdp_")
test_chlist <- list()
for (out in outs) {
  test_chlist[[out]] <- eval(parse(text = out))
}


test_pr_multi <- hdp_multi_chain(test_chlist)
test_pr_ec <- hdp_extract_components(test_pr_multi, cos.merge = 0.9)


###
tnt <- read.csv("sbs_category_counts.txt", sep="\t", header = T, stringsAsFactors =F, row.names = 1)
patients <- sapply(strsplit(rownames(tnt), ":"), "[[", 1)
sigs <- read.csv("PCAWG_sigProfiler_SBS_signatures.csv", header = T, row.names = 1, stringsAsFactors = F)
gdsigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS25", "SBS28", "SBS30", "SBS37", "SBS40", "SBS41", "SBS43", "SBS45", "SBS49") 
sigs <- as.matrix(sigs[,gdsigs])

###

toutmeans <- (comp_dp_distn(test_pr_ec)$mean)
mymeans <- (toutmeans[(nrow(toutmeans)-nrow(tnt) + 1):nrow(toutmeans),])
rownames(mymeans) <- rownames(tnt)

# rename the signatures - because of a and b sigs in the file the labelling has become confused. 
tochange <- colnames(mymeans)[grep("P", colnames(mymeans))]
colnames(mymeans)[grep("P", colnames(mymeans))] <- colnames(sigs)[as.numeric(gsub("P", "", tochange))]

write.table(mymeans, "HDP_conditioned_on_pcawg_sigs_assignments_cosmerge_pt9.txt", sep="\t", col.names = T, row.names = T, quote=F)

mm <- as.data.frame(mymeans)
mm$patient <- sapply(strsplit(rownames(mymeans), ":"), "[[", 1)


#########
mut_count <- tnt

pdf("HDP_conditioned_on_PCAWG_plots_cosmerge_pt9.pdf")
par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
p1 <- lapply(chains(test_pr_multi), plot_lik, bty="L", start=10000) # would really need a longer burn in.
p2 <- lapply(chains(test_pr_multi), plot_numcluster, bty="L")
p3 <- lapply(chains(test_pr_multi), plot_data_assigned, bty="L")

par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(test_pr_ec, bty="L")

bases <- c("A", "C", "G", "T")
trinuc_context <- paste0(rep(rep(bases, times=6), each=4),
                         rep(c("C", "T"), each=48),
                         rep(bases, times=24))
group_factor <- as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each=16))

for (i in 1:ncol(mymeans)) {
  plot_comp_distn(test_pr_ec, comp=i-1, cat_names=trinuc_context,
                  grouping=group_factor, col=c("blue", "black", "red", "grey", "limegreen", "pink"),
                  col_nonsig="lightcyan", show_group_labels=TRUE, plot_title=colnames(mymeans)[i])
}

plot_dp_comp_exposure(test_pr_ec, 1+ncol(sigs)+1+(1:nrow(mut_count)), incl_nonsig = FALSE,
                      col=rep(RColorBrewer::brewer.pal(9, "Set1"), 4), incl_numdata_plot = F)
dev.off()


sigcategs <- t(comp_categ_distn(test_pr_ec)$mean)
write.table(sigcategs, "Signature_category_counts.txt", sep="\t", col.names=T, row.names = F, quote=F)