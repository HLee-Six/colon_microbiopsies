
## This file contains various functions that are used in the ABC simulation of crypt fission rate in the 
## normal colon.
## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk

if(file.exists("/nfs/users/nfs_s/so11/phd/doIexist")) {
  ## Script is running on the farm
  lib_location="/lustre/scratch114/projects/crohns/somatic_ibd_p1/ABC/R_packages_farm4/"
  library(crayon, lib.loc=lib_location)
  library(backports, lib.loc=lib_location)
  library(vctrs, lib.loc=lib_location)
  library(data.tree, lib.loc=lib_location)
  library(ape, lib.loc=lib_location)
  #library(phytools, lib.loc=lib_location)
  library(ggtree, lib.loc=lib_location)
  .libPaths( c(lib_location, .libPaths()) )
} else {
  ## Running the script locally - Testing
  library(data.tree)
  library(ape)
  library(phytools)
  library(ggtree)
}



## This function takes the x and y coordinates of a crypt and returns the x and y coordinates of
## one of its 8 neighbors on the grid by random.
## The x and y coordinates are integers in {2..(3*largest_distance-1)}
choose_neighbor <- function(x_coordinate, y_coordinate) {
  
  x_allNeighbors <- c(rep(x_coordinate-1, 3), x_coordinate, rep(x_coordinate+1, 3), x_coordinate)
  y_allNeighbors <- c(y_coordinate+1,y_coordinate, rep(y_coordinate-1, 3), y_coordinate, rep(y_coordinate+1,2))
  neighbor_nr <- sample(c(1:8), 1)
  return(as.numeric(c(x_allNeighbors[neighbor_nr], y_allNeighbors[neighbor_nr])))
}




## Input: crypt_fission_rate is a float and has the units fissions per crypt per year.
##        largest_distance is the largest distance seen for any two crypts for the biopsy being simulated
##        patient_age is the age of the patient from which the biopsy was taken.
## Output: A dataframe containing the whole history of the grid throughout the simulation. 
## Algorithm:
##            1. Define a grid with dimensions that are three times as large as the largest distance seen for the biopsy
##            2. Draw the time until next crypt fission from exponential distribution with expected value
##            crypt_fission_rate times the number of crypts that can undergo fission. 
##            3. Randomly choose a crypt from the grid to undergo fission and randomly choose one of it's neighbors
##            to die. 
##            4. Add two new lines to the dataframe, representing these new crypts. Update the death times of the original two
simulate_colon <- function(crypt_fission_rate=0.0333, largest_distance=30, patient_age=50, mySeed=42) {
  
  set.seed(mySeed)
  ## Define grid of size n x n, where n is three times the largest_distance.
  ## Assign coordinates and label cells that are on the edges.
  grid <- data.frame(c(1:((3*largest_distance)^2)), rep(1:(3*largest_distance), times=(3*largest_distance)), rep(1:(3*largest_distance), each=(3*largest_distance)))
  colnames(grid) <- c("cryptID", "x_coordinate", "y_coordinate")
  grid$onEdge <- 0
  grid$onEdge[grid$x_coordinate==1 | grid$x_coordinate==(3*largest_distance) | grid$y_coordinate==1 | grid$y_coordinate==(3*largest_distance) ] <- 1
  ## Parent of the original set of crypts is a pseudonode nr 0
  grid$parent <- 0
  grid$birth_time <- 0
  grid$death_time <- patient_age
  grid$pathString <- paste(grid$parent, grid$cryptID, sep="/")
  
  time <- 0
  Count <- nrow(grid)
  while(time < patient_age) {
    ## Sample the time until next fission event
    ## The time between ANY fissions is the fission rate times the number of crypts in the grid that can undergo fission.
    ## Remember that crypts on the edge can't undergo fission. 
    t <- rexp(1, rate=crypt_fission_rate*(((3*largest_distance)^2)-(4*3*largest_distance - 4)))
    time <- time + t
    
    ## Stop right away if too long time has passed. This will prevent the last two crypts from having birth time higher than the
    ## age of the patient. 
    if(time > patient_age) {
      break
    }
    
    ## Chose one crypt to divide (referred to as the parent) and one of its neighbors to die off.
    lineNr_of_dividing_crypt <- sample(grid$cryptID[grid$death_time==patient_age], 1)
    dividing_crypt_coord <- grid[lineNr_of_dividing_crypt, c("x_coordinate", "y_coordinate")]
    
    
    ## If the crypt chosen falls on the edge of the grid, choose again.
    while(grid$onEdge[lineNr_of_dividing_crypt]==1) {
      lineNr_of_dividing_crypt <- sample(grid$cryptID[grid$death_time==patient_age], 1)
      dividing_crypt_coord <- grid[lineNr_of_dividing_crypt, c("x_coordinate", "y_coordinate")]
    }
    crypt_to_die_coord <- choose_neighbor(dividing_crypt_coord[1], dividing_crypt_coord[2])
    ## The crypt in the dataframe that dies is the one that has the right coordinates and the highest cryptID, as that one
    ## is the only one that hasn't yet been replaced by another crypt.
    lineNr_of_dying_crypt <- max(which(grid$x_coordinate==crypt_to_die_coord[1] & grid$y_coordinate==crypt_to_die_coord[2]))

    ## Record the fission event by adding to the grid-history dataframe
    # Reminder: The columns are: cryptID x_coordinate y_coordinate onEdge parent birth_time death_time pathString
    tmp <- data.frame(c(Count+1, dividing_crypt_coord, 0, grid$cryptID[lineNr_of_dividing_crypt], time,patient_age, paste(grid$pathString[lineNr_of_dividing_crypt],Count+1, sep="/")))
    colnames(tmp) <- colnames(grid)
    grid <- rbind(grid, tmp)
    
    # The dying crypt could be on the edge, although the dividing crypt can not. 
    isEdgy <- ifelse((1 %in% crypt_to_die_coord) | ((3*largest_distance) %in% crypt_to_die_coord), 1, 0)
    tmp <- data.frame(c(Count+2, data.frame(t(crypt_to_die_coord)), isEdgy, grid$cryptID[lineNr_of_dividing_crypt], time,patient_age, paste(grid$pathString[lineNr_of_dividing_crypt],Count+2, sep="/")))
    colnames(tmp) <- colnames(grid)
    grid <- rbind(grid, tmp)
    
    ## Update time of death for both crypts
    grid$death_time[lineNr_of_dividing_crypt] <- time
    grid$death_time[lineNr_of_dying_crypt] <- time
    

    Count <- Count + 2
    
  }
  
  return(grid)
}





## Input: spat_rel_df is a spatial relationship matrix but one that exists as a dataframe object in R.
##        It should not have a header but it has an ID column which contains crypt IDs.
##        The matrix should contain 0s on the diagonal line, crypt distances in the upper triangle
##        and can contain anything (most likely text on the form LEFT-LEFT) in the lower triangle.
## Output: A propper matrix object with both triangles filled out. 
format_matrix <- function(spat_rel_df) {
  colnames(spat_rel_df) <- c("ID", as.character(spat_rel_df[, c(1)]))
  rownames(spat_rel_df) <- spat_rel_df$ID
  rn <- rownames(spat_rel_df)
  spat_rel_df$ID <- NULL
  
  ## Expect warning messages as strings can't be numeric. That's fine. NAs in the lower triangle
  ## will be overwritten by values in the upper triangle. The matrix is necessarily symmetrical. 
  spat_rel_df <- data.frame(lapply(spat_rel_df, function(x) {as.numeric(as.character(x))}))
  rownames(spat_rel_df) <- rn
  spat_rel_mat <- as.matrix(spat_rel_df)
  spat_rel_mat <- as.matrix(Matrix::forceSymmetric(spat_rel_mat,uplo="U"))
  
  return(spat_rel_mat)
}

## Input: spat_rel_mat is a numeric matrix of crypt distances - like the one produced by format_matrix()
##        get_order is a boolean parameter specifying if the user wants the function to return a character vector giving
##        the spatial ordering of the crypts in spat_rel_mat, or a dataframe containing the distances on the x and y axes. 
## Output: The output depends on the value of the get_order parameter. 
##          If this is TRUE, then a character vector is returned giving the spatial order of the crypts along the y-axis. 
##          Call this character vector crypt_order
##          If get_order=FALSE, a dataframe of two numeric vectors, x_dist and y_dist, is returned. This dataframe
##          contains one fewer rows than spat_rel_mat and contains the distance from the first crypt in crypt_order to all
##          other crypts. 
get_crypt_distances <- function(spat_rel_mat, get_order=FALSE) {
  
  if(nrow(spat_rel_mat)<3) {
    ## k = 1. Distances between 2 crypts only in 1 dimension. 
    ## Everyone likes 6-nested functions, right?
    y_dist <- cumsum(unname(round(diff(sort(cmdscale(spat_rel_mat, k = 1)[, 1])))))
    x_dist <- rep(0, (nrow(spat_rel_mat)-1))
    crypt_names <- names(sort(cmdscale(spat_rel_mat, k = 1)[, 1]))
  } else {
    ## Use k = 2 to find distances in two dimensions. 
    mds <- round(cmdscale(spat_rel_mat, k=2))
    crypt_order <- order(mds[,1])
    y_dist <- cumsum(unname(diff(mds[crypt_order,1])))
    x_dist <- cumsum(unname(diff(mds[crypt_order,2])))
    crypt_names <- rownames(mds)[crypt_order]
  }
  
  if(get_order) {
    return(crypt_names)
  } else {
    return(data.frame(x_dist, y_dist))
  }

}



## The purpose of this function is to sample crypts from a slice of the grid in such a way that
## it reflects Henry's sampling
## Input: slice is a dataframe with the same columns as grid. All the crypts in <slice> have the same x-coordinate.
##        spatial_rel_matrix is a spatial relationship matrix where the row numbers are the crypt IDs and the 
##        matrix contains integer values representing the one-dimensional distance between the crypts. 
##        Alternatively, spatial_rel_matrix can be left empty and the crypt_distances given as a numeric vector
##        of distances. You don't actually need the crypt_order vector for now... 
## Output: Returns a subset of slice that has the size nrow(spatial_rel_matrix)
## Note: Currently can't handle missing values in the spatial_rel_matrix
##        There are multiple instances in Henry‘s data where there are several groups of crypts within the biopsy. 
##        We know the spatial relationship between crypts within a group but not between groups. 
##        We don‘t know the relative positions or the orientation of the groups. 
take_representative_sample <- function(grid, ref_coordinates, crypt_distances_x, crypt_distances_y, mySeed=42) {
  
  set.seed(mySeed)
  
  ## To the left or the right of the "reference crypt"
  direction_x <- sample(c(1, -1), size=1)
  direction_y <- sample(c(1, -1), size=1)
  crypt_distances_x <- crypt_distances_x*direction_x
  crypt_distances_y <- crypt_distances_y*direction_y
  
  patientAge <- max(grid$death_time)
  grid <- subset(grid, death_time==patientAge)
  # Return the crypts that match these positions
  sampled_crypts <- grid[grid$x_coordinate==ref_coordinates[1] & grid$y_coordinate==ref_coordinates[2],]
  for(j in seq_along(crypt_distances_x)) {
    sampled_crypts <- rbind(sampled_crypts, grid[grid$x_coordinate==(ref_coordinates[1]+crypt_distances_x[j]) & grid$y_coordinate==(ref_coordinates[2]+crypt_distances_y[j]), ])
  }
  
  return(sampled_crypts)
}






## Input: This function takes a grid-history dataframe and a list of crypts that have been sampled
##        from the final-state (for example the output of take_representative_sample())
##        The columns of both grid and cryptList are:
##        cryptID x_coordinate y_coordinate onEdge parent birth_time death_time pathString
## Output: It returns a phylo object containing the phylogenic tree of the sampled crypts. 
##        The branch.length of this object contains the death_times of all the nodes, which
##        for internal nodes correspond to coalescence times. 
make_tree <- function(grid, cryptList) {
  
  patientAge=max(grid$death_time)
  ## Use the pathString to convert the cryptList into a data.tree object
  tree_obj <- as.Node(cryptList)
  ## That can then be converted to a phylo object. Game is won. 
  tree_phylo <- as.phylo(tree_obj)
  
  ## Collapse those nodes that link to only one branch. These are anchestors but don't represent
  ## coalescent events
  tree_collapsed <- ape::collapse.singles(tree_phylo)
  tree_df <- fortify(tree_collapsed)
  
  ## Add the death times as attributes to the tree and make the branch length
  ## reflect the time. 
  tree_df$death_time <- NA
  tree_df$death_time[1:nrow(cryptList)] <- patientAge
  tree_df$death_time[nrow(cryptList)+1] <- 0  ## Root
  
  if(nrow(tree_df)!=(nrow(cryptList)+1)) {
    ## There are some internal nodes in the tree and I should update their death times
    cryptIDs_internal_nodes <- as.numeric(tree_df$label[(nrow(cryptList)+2):nrow(tree_df)])
    index <- (nrow(cryptList)+2)
    for(id in cryptIDs_internal_nodes) {
      tree_df$death_time[index] <- grid[grid$cryptID==id, c("death_time")]
      index <- index+1
    }
  }

  
  tree_collapsed$edge.length=tree_df$death_time[tree_collapsed$edge[,2]]
  tree_collapsed$death_time <- tree_df$death_time
  
  return(tree_collapsed)
}



## Input: myTree is a phylogenic tree (for example the output of make_tree) that has the
##        death times of the nodes in the branch.length attribute and specifies tips with TRUE/FALSE in the isTip attribute.
##        timepoints is a numeric vector with 7 elements.
##        timeIncl_limit is a number >=0. Events happening before timeIncl_limit will not be counted. 
##        The point is to allow us to exclude early neonatal expansions. 
## Output: A numeric vector containing the number of coalescence events in each time-window
count_event_times <- function(myTree, timepoints=c(10,20,30,40,50,60,70), timeIncl_limit=0) {
  
  tree_df <- fortify(myTree)
  nrTips <- sum(fortify(tree_df)$isTip)
  # Remove the root and the tips
  tree_df <- tree_df[-c(1:(nrTips+1)), ]
  
  ## only count events that happen after the timeIncl_limit. Allows us to exclude
  ## early neonatal expansions
  tree_df <- subset(tree_df, branch.length>timeIncl_limit)
  
  # How many coalescence events do we observe in each time-slot
  t1 <- nrow(subset(tree_df, branch.length < timepoints[1]))
  t2 <- nrow(subset(tree_df, branch.length >= timepoints[1] & branch.length < timepoints[2]))
  t3 <- nrow(subset(tree_df, branch.length >= timepoints[2] & branch.length < timepoints[3]))
  t4 <- nrow(subset(tree_df, branch.length >= timepoints[3] & branch.length < timepoints[4]))
  t5 <- nrow(subset(tree_df, branch.length >= timepoints[4] & branch.length < timepoints[5]))
  t6 <- nrow(subset(tree_df, branch.length >= timepoints[5] & branch.length < timepoints[6]))
  t7 <- nrow(subset(tree_df, branch.length >= timepoints[6] & branch.length < timepoints[7]))
  t8 <- nrow(subset(tree_df, branch.length >= timepoints[7]))
  
  counts <- c(t1, t2, t3, t4, t5, t6, t7, t8)
  return(counts)
}








