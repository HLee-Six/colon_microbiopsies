
## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk


if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  source("/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
  output_dir="/lustre/scratch119/humgen/teams/anderson/users/so11/somatic_ibd_p1/ABC/normal_fission_rate/"
} else {
  source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")
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

  mySeed <- mySeed+1
}








