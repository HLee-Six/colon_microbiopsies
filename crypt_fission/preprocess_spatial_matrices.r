
## The purpose of this script is to read in individual spatial relationship matrices and output a file containing a list
## of all necessary information for each matrix=biopsy.

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

source("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/abc_function_archive.r")

spatial_data_dir <- "/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/spatial_relationships/biopsy_matricies/"
output_dir <- "/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/"
biopsy_list <- read.table("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/spatial_relationships/biopsy_matricies/biopsy_list_combined.txt", h=F, stringsAsFactors = F)

colnames(biopsy_list) <- c("biopsy", "patient_age")

biopsy_list$largest_y_distance <- NA
biopsy_list$largest_x_distance <- NA
biopsy_list$crypt_distances_y <- NA
biopsy_list$crypt_distances_x <- NA
biopsy_list$crypt_list <- NA
for(i in 1:nrow(biopsy_list)) {
  
  assign(paste("P", biopsy_list$biopsy[i], sep=""), read.table(paste(spatial_data_dir, biopsy_list$biopsy[i], ".txt", sep="")))
  ## Convert df to matrix
  assign(paste("P", biopsy_list$biopsy[i], "m", sep=""), format_matrix(get(paste("P", biopsy_list$biopsy[i], sep=""))))

  ## Find largest distance between two crypts in the matrix
  biopsy_list$largest_y_distance[i] <- max(abs(get_crypt_distances(get(paste("P", biopsy_list$biopsy[i], "m", sep="")))[,2]))
  biopsy_list$largest_x_distance[i] <- max(abs(get_crypt_distances(get(paste("P", biopsy_list$biopsy[i], "m", sep="")))[,1]))
  
  ## Get a really strange conversion of 6 to 6.00000000000001 when i=10. Put the round here although risky - could mask other bugs
  biopsy_list$crypt_distances_x[i] <- paste(round(get_crypt_distances(get(paste("P", biopsy_list$biopsy[i], "m", sep="")))[,1]), collapse = ",")
  biopsy_list$crypt_distances_y[i] <- paste(round(get_crypt_distances(get(paste("P", biopsy_list$biopsy[i], "m", sep="")))[,2]), collapse = ",")
  
  biopsy_list$crypt_list[i] <- paste(get_crypt_distances(get(paste("P", biopsy_list$biopsy[i], "m", sep="")), get_order=T), collapse = ",")
}

biopsy_list$largest_distance <- ifelse(biopsy_list$largest_y_distance>=biopsy_list$largest_x_distance,biopsy_list$largest_y_distance, biopsy_list$largest_x_distance )

write.table(biopsy_list, file=paste(output_dir, "biopsy_data_combined.txt", sep=""), row.names = F, quote = F)


## Split the biopsy_list into 10 parts such that each part has approx the same sum(largest_distance)
## was playing around with this to parallelize more but all code after this point can actually be ignored.
total <- sum(biopsy_list$largest_distance)
stub <- round(total/10)
sub_sum <- 0
j <- 1
count <- 1
for( i in 1:nrow(biopsy_list)) {
  if(sub_sum>stub) {
    write.table(biopsy_list[c(j:i),], file=paste(output_dir, "biopsy_data_", count, ".txt", sep=""), row.names = F, quote = F)
    sub_sum <- 0
    count <- count +1
    j <- i+1
  }
  sub_sum <- sub_sum + biopsy_list$largest_distance[i]

  if(i==nrow(biopsy_list)) {
    write.table(biopsy_list[c(j:i),], file=paste(output_dir, "biopsy_data_", count, ".txt", sep=""), row.names = F, quote = F)
  }
}



