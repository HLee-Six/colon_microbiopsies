
## The purpose of this script is to time the coalescent events in the observed data and write out a dataframe of event_counts that
## can be compared with simulated data. 

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

library(ape)
library(ggtree)
source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")

tree_data_dir <- "/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/tree_by_each_sig/"

biopsy_list <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data_combined.txt", h=T, stringsAsFactors = F)

biopsy_list$patient <- unlist(strsplit(biopsy_list$biopsy, split = "_"))[c(TRUE, FALSE)]
biopsy_list$location <- "Left"
biopsy_list$location[grep("trans", biopsy_list$biopsy)] <- "Trans"
biopsy_list$location[grep("right", biopsy_list$biopsy)] <- "Right"
biopsy_list$location[grep("ileum", biopsy_list$biopsy)] <- "Ileum"


## SBS1 mutation rates as stated in the manuscript
mut_rate_right <- 16.8
mut_rate_trans <- 16.1
mut_rate_left <- 12.8
mut_rate_ileum <- 12.7 

event_counts <- data.frame()
for(i in 1:nrow(biopsy_list)) {
  
  ## O438_trans just keeps failing. Distance between crypts is 45. Almost certainly no events but better skip.
  if(i==43) {
    next
  }
  tree_tmp <- read.tree(paste(tree_data_dir, biopsy_list$patient[i], "/", biopsy_list$patient[i], "_sbs_SBS1.tree", sep=""))
  
  ## Extract from the trees the crypts that belong to the biopsy we're testing. 
  ## If there are any coalescent events linking crypts from different "crypt-groups" within a biopsy
  ## then these will be removed so the comparison with the simulated data is valid. 
  ## Similarly, if there are crypts with all-missing spatial relationship data then these will be excluded as well. 
  simulation_crypts <- unlist(strsplit(biopsy_list$crypt_list[i], split=","))
  pruned.tree<-drop.tip(tree_tmp,tree_tmp$tip.label[-match(simulation_crypts, tree_tmp$tip.label)])
  
  
  ## Convert mutation numbers into years
  if(biopsy_list$location[i]=="Right") {
    pruned.tree$edge.length <- pruned.tree$edge.length/mut_rate_right
  } else if(biopsy_list$location[i]=="Trans") {
    pruned.tree$edge.length <- pruned.tree$edge.length/mut_rate_trans
  } else if(biopsy_list$location[i]=="Left") {
    pruned.tree$edge.length <- pruned.tree$edge.length/mut_rate_left
  } else if(biopsy_list$location[i]=="Ileum") {
    pruned.tree$edge.length <- pruned.tree$edge.length/mut_rate_ileum
    } else {
    stop(paste("Unrecognized biopsy location for biopsy", i))
  }
  
  ## For each node that isn't a leaf and not the root, count it. 
  event_counts <- rbind(event_counts, count_event_times(pruned.tree, timeIncl_limit = 0))
  #ggtree(pruned.tree) +theme_tree2()+xlim(0,max(fortify(pruned.tree)$x)*1.1) 
  
}

colnames(event_counts) <- c(paste("t", 1:8, sep=""))
event_counts$biopsy <- biopsy_list$biopsy[-43]

write.table(event_counts, file="/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/observed_event_counts_excl0yr.txt", row.names = F, quote = F)



