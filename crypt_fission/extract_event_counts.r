
## The purpose of this script is to count the number of coalescence events occurring in the simulated data. 

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  source("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
  #load("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data.RData")
  output_dir="/lustre/scratch114/projects/crohns/somatic_ibd_p1/ABC/normal_fission_rate/"
} else {
  source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
  #load("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data.RData")
  output_dir <- "/Users/so11/phd/so11_lustre/somatic_ibd_p1/ABC/normal_fission_rate/"
}


parameters <- commandArgs(TRUE)
## Crypt fission rate in the unit fissions per crypt per year.
crypt_fission_rate <- as.numeric(parameters[1])

## Read in the index for this job in the job-array. 
simulation_nr <- as.numeric(parameters[2])
biopsy_list <- read.table(parameters[3], h=T, stringsAsFactors = F)
tree_dir <- parameters[4]
grid_dir <- parameters[5]
event_times_dir <- parameters[6]

## Skip O438_trans as it keeps failing
biopsy_list <- biopsy_list[-43,]


event_counts_excl0yr <- data.frame()
event_counts_excl1yr <- data.frame()
event_counts_excl2yr <- data.frame()
event_counts_excl3yr <- data.frame()
event_counts_excl4yr <- data.frame()
tree_comb_df <- data.frame()
for(i in 1:nrow(biopsy_list)) {
  ## Need a patch to deal with the fact that for a while, PD34201_trans4 appeared twice in the list.
  ## To preserve the seeds I need to add 1 to the iterator.
  if(i>49) {
    ## A separate seed for each biopsy and each iteration.
    mySeed=simulation_nr + i + (simulation_nr -1)*nrow(biopsy_list)
  } else {
    ## A separate seed for each biopsy and each iteration.
    mySeed=simulation_nr + (i - 1) + (simulation_nr -1)*nrow(biopsy_list)
  }
  set.seed(mySeed)
  
  ## If the grid for this biopsy exists, read it. Otherwise, write that it is missing. 
  if(file.exists(paste(grid_dir,"sim", simulation_nr, "/grid_", biopsy_list$biopsy[i], "_sim", simulation_nr ,sep=""))) {
    grid <- read.table(paste(grid_dir, "sim", simulation_nr, "/grid_", biopsy_list$biopsy[i], "_sim", simulation_nr ,sep=""), h=T, stringsAsFactors = F)
    
  } else {
    
    system(paste("echo \"", biopsy_list$biopsy[i]," ", simulation_nr," ", crypt_fission_rate, "\" >> ", output_dir, "unfinished_grids.txt", sep=""))
    next
  }
  
  
  
  ## Select x and y coordinates for the crypt that will serve as reference for sampling. 
  ## randomly choose a position for the 'reference crypt', that is the first crypt in crypt_order.
  ## The position must be chosen so that we don't get crypt positions outside the edges of the grid.
  ## And the reference should not fall on the edge of the grid.
  ref_x_coord <- sample(c((1+biopsy_list$largest_distance[i]):(2*round(biopsy_list$largest_distance[i]))), 1)
  ref_y_coord <- sample(c((1+biopsy_list$largest_distance[i]):(2*round(biopsy_list$largest_distance[i]))), 1)
  ## Unexpected behavior of the sample() function. When length(x)==1 sampling takes place from 1:x
  while(ref_x_coord==1 | ref_y_coord==1) {
    ref_x_coord <- sample(c((1+biopsy_list$largest_distance[i]):(2*round(biopsy_list$largest_distance[i]))), 1)
    ref_y_coord <- sample(c((1+biopsy_list$largest_distance[i]):(2*round(biopsy_list$largest_distance[i]))), 1)
  }
  
  
  ## Sample grid in a way that is representative of the actual LCM sampling
  sampled_crypts <- take_representative_sample(grid,ref_coordinates =  c(ref_x_coord, ref_y_coord), crypt_distances_x=as.numeric(unlist(strsplit(biopsy_list$crypt_distances_x[i], split=","))),
                                               crypt_distances_y=as.numeric(unlist(strsplit(biopsy_list$crypt_distances_y[i], split = ","))), mySeed=mySeed)
  
  ## Make the phylogenic tree for the sampled crypts
  sampled_tree <- make_tree(grid = grid, cryptList = sampled_crypts) 
  tree_df <- fortify(sampled_tree)
  tree_df$biopsy <- biopsy_list$biopsy[i]
  tree_comb_df <- rbind(tree_comb_df, tree_df)
  
  ## Count the number of coalescence events in the tree
  ## Calculate the number of coalescence events that occur in each 10y time window. Write this to a results dataframe
  
  branching_events <- count_event_times(sampled_tree, timeIncl_limit = 0)
  event_counts_excl0yr <- rbind(event_counts_excl0yr, branching_events)

  branching_events <- count_event_times(sampled_tree, timeIncl_limit = 1)
  event_counts_excl1yr <- rbind(event_counts_excl1yr, branching_events)
  
  branching_events <- count_event_times(sampled_tree, timeIncl_limit = 2)
  event_counts_excl2yr <- rbind(event_counts_excl2yr, branching_events)
  
  branching_events <- count_event_times(sampled_tree, timeIncl_limit = 3)
  event_counts_excl3yr <- rbind(event_counts_excl3yr, branching_events)
  
  branching_events <- count_event_times(sampled_tree, timeIncl_limit = 4)
  event_counts_excl4yr <- rbind(event_counts_excl4yr, branching_events)
  

  mySeed <- mySeed+1
}

write.table(event_counts_excl0yr, file=paste(event_times_dir, "event_counts_", simulation_nr, "_excl0yr.txt", sep=""), quote = F, row.names = F)
write.table(event_counts_excl1yr, file=paste(event_times_dir, "event_counts_", simulation_nr, "_excl1yr.txt", sep=""), quote = F, row.names = F)
write.table(event_counts_excl2yr, file=paste(event_times_dir, "event_counts_", simulation_nr, "_excl2yr.txt", sep=""), quote = F, row.names = F)
write.table(event_counts_excl3yr, file=paste(event_times_dir, "event_counts_", simulation_nr, "_excl3yr.txt", sep=""), quote = F, row.names = F)
write.table(event_counts_excl4yr, file=paste(event_times_dir, "event_counts_", simulation_nr, "_excl4yr.txt", sep=""), quote = F, row.names = F)


write.table(tree_comb_df, file=paste(tree_dir, "trees_", simulation_nr,".txt", sep=""), quote = F, row.names = F)



