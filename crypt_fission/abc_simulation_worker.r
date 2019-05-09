
## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk


if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  source("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
  load("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data.RData")
  output_dir="/lustre/scratch119/humgen/teams/anderson/users/so11/somatic_ibd_p1/ABC/normal_fission_rate/"
} else {
  source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
  load("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/biopsy_data.RData")
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

second_batch=T
event_counts <- data.frame()
tree_comb_df <- data.frame()
for(i in 1:nrow(biopsy_list)) {

  print(paste("Running",i))
  
  ## When I'm running the second batch of biopsies, I need to make sure the seed is unique
  if(second_batch) {
    ## Need a patch to deal with the fact that for a while, PD34201_trans4 appeared twice in the list.
    ## To preserve the seeds I need to add 1 to the iterator.
    if(i>49) {
      ## A separate seed for each biopsy and each iteration.
      mySeed=simulation_nr + i+66 + (simulation_nr -1)*(nrow(biopsy_list)+66)
    } else {
      ## A separate seed for each biopsy and each iteration.
      mySeed=simulation_nr + (i - 1 + 66) + (simulation_nr -1)*(nrow(biopsy_list)+66)
    }
  } else {
    ## Need a patch to deal with the fact that for a while, PD34201_trans4 appeared twice in the list.
    ## To preserve the seeds I need to add 1 to the iterator.
    if(i>49) {
      ## A separate seed for each biopsy and each iteration.
      mySeed=simulation_nr + i + (simulation_nr -1)*nrow(biopsy_list)
    } else {
      ## A separate seed for each biopsy and each iteration.
      mySeed=simulation_nr + (i - 1 ) + (simulation_nr -1)*nrow(biopsy_list)
    }
  }

  set.seed(mySeed)
  
  if(file.exists(paste(grid_dir, "grid_", biopsy_list$biopsy[i], "_sim", simulation_nr, sep=""))) {
    next
  } else {
    ## Simulate a colon
    ## Had this really weird problem where the numbers are actually not integers but stored as say 2.99999999999 rather than 3. Use round
    grid <- simulate_colon(crypt_fission_rate=crypt_fission_rate, largest_distance = round(biopsy_list$largest_distance[i]), 
                           patient_age = biopsy_list$patient_age[i], mySeed = mySeed)
    
    ## Save the grid
    write.table(grid, file=paste(grid_dir, "grid_", biopsy_list$biopsy[i], "_sim", simulation_nr, sep=""), row.names = F, quote = F)
  }

  if(FALSE) {
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


  ## Sample the one-dimensional slice
  sampled_crypts <- take_representative_sample(grid,ref_coordinates =  c(ref_x_coord, ref_y_coord), crypt_distances_x=as.numeric(unlist(strsplit(biopsy_list$crypt_distances_x[i], split=","))),
                                               crypt_distances_y=as.numeric(unlist(strsplit(biopsy_list$crypt_distances_y[i], split = ","))), mySeed=mySeed)

  ## Make the phylogenic tree for the sampled crypts
  sampled_tree <- make_tree(grid = grid, cryptList = sampled_crypts)
  tree_df <- fortify(sampled_tree)
  tree_df$biopsy <- biopsy_list$biopsy[i]
  tree_comb_df <- rbind(tree_comb_df, tree_df)

  ## Count the number of coalescence events in the tree
  ## Calculate the number of coalescence events that occur in each 10y time window. Write this to a results dataframe
  branching_events <- count_event_times(sampled_tree)
  event_counts <- rbind(event_counts, branching_events)
  }

  mySeed <- mySeed+1
}

## Save the dataframes containing the results for all biopsies. 
colnames(event_counts) <- c(paste("t", 1:8, sep=""))
event_counts$biopsy <- biopsy_list$biopsy


## If files exist, append
## Not doing this ATM. Need to extract trees and event counts in a separate step.
if(FALSE) {
  if(file.exists(paste(event_times_dir, "event_counts_", simulation_nr, ".txt", sep=""))) {
  write.table(event_counts, file=paste(event_times_dir, "event_counts_", simulation_nr, ".txt", sep=""), quote = F, row.names = F, col.names = F, append=T)

} else {
  write.table(event_counts, file=paste(event_times_dir, "event_counts_", simulation_nr, ".txt", sep=""), quote = F, row.names = F)
}

if(file.exists(paste(tree_dir, "trees_", simulation_nr,".txt", sep=""))) {
  write.table(tree_comb_df, file=paste(tree_dir, "trees_", simulation_nr,".txt", sep=""), quote = F, row.names = F, col.names = F, append=T)
} else {
  write.table(tree_comb_df, file=paste(tree_dir, "trees_", simulation_nr,".txt", sep=""), quote = F, row.names = F)
}
}







