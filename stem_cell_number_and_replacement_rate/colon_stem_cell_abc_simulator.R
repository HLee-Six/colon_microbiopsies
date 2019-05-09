# 2019.03.13_colon_abc_simulator.R

rootdir <- "~/colon_stem_cell_abc"
setwd(rootdir)

# sample the number of stem cells
stem <- sample(x=c(1:30),size = 1)
# bromodeoxyuridine measurements indicate that the cell cycle time at the base of human colonic crypts is of the order of 30 hours (Potten et al. 1992).
# but a proportion of these divisions could be asymmetrical, and so the time between symmetrical cell divisions (i.e. the Moran generation time) would be longer.
gendays <- sample(x=c(1:620),size = 1)

# get the total number of stem cells per crypt
cryptcells <- 2000

# decide the total number of mutations per doubling cell division (i.e. to go from stem cells to colonocytes). Try playing around with this to see if it changes anything
mutsperdoubling <- 1

# the mutation rate: the number of somatic mutations per year. 
myconstants <- apply(combn(c(letters, LETTERS), 2), 2, function(col) paste(as.vector(col), collapse=""))
add.mutns.fn <- function(cell, mutrate) {
  newmuts <- myconstants[1:rpois(1, lambda=mutrate)]
  withmuts <- paste(paste0(cell, newmuts), collapse="-")
  return(withmuts)
}
age <- 55

# create a unique id
myrand <- round((runif(1, min=0, max=1000)*100)+ runif(1, min=0, max=100))
print(paste0("Setting random seed: ", myrand))
# set this as the seed, and print out the seed. That way should be able to replicate any results, I hope.
set.seed(myrand)

id <- paste0(stem, "_stemcells_", gendays, "_gendays_", myrand, "_randseed_")
print(paste0("ID for this simulation is ", id))

dir.create(id)
setwd(paste0(rootdir,  "/", id))

genperyear <- (365/gendays)
fwgens <- genperyear*age
print(paste0("Number of fisher wright generations is ", fwgens))

# read in some unique ids. NB will have to change path if run this on the farm.
# uids are stored in a file with 1000 uids per row.
rowsofuids <- ceiling(stem/1000)+1

myuids <- read.csv("~/input_files/uids.txt",
                     header=F, sep=" ", stringsAsFactors = F, nrows = rowsofuids)

myuids <- as.vector(unlist(myuids))
myuids <- myuids[myuids != ""]
myuids <- myuids[!is.na(myuids)] # there's one NA in there that I think can screw everything up.
print(paste0(length(myuids), " unique ids made")) 

inputfiles <- list.files("~/input_files/",
                         pattern="input_data.txt", full.names = T)

######################
# from now on, it is crypt specific
allstats <- c(stem,gendays,mutsperdoubling)

for (inputfile in inputfiles) {
  
  # add a new id specific to this test crypt
  ifsplit <- strsplit(inputfile, "/")
  ifsplit <- as.vector(unlist(ifsplit))
  inputcrypt <- gsub("_input_data.txt", "", ifsplit[length(ifsplit)])
  tid <- paste0(id, inputcrypt, "_crypt_")
  
  # read in input data
  indata <- read.csv(inputfile, sep=" ", stringsAsFactors = F, header=F)
  
  # NB clonalmuts is the number of clonalmuts that COULD CALL - i.e. that ended up with at least 3 reads.
  # may need to adjust this for the proportion of callable mutations.
  clonalmuts <- indata[1,2] 
  mutrate <- clonalmuts/age
  
  # add a certain amount of contamination.
  stromaprop <- (0.5 - indata[1,1])*2
  
  # run the FW model for the specified number of generations.
  # with my new prefixing strategy I only need as many uids as there are hscs in the sample.
  # because crypt is probably only seeded by one cell, start off with all the stem cells sharing one mutation
  ngen <- rep(myuids[1], stem)
  
  for (f in 1:fwgens) {
    tgen <- sample(ngen, replace=T)
    newgen <- myuids[1:length(tgen)]
    
    # add a prefix related to this bin to all the ids in this bin.
    prefnewgen <- paste0(f, newgen)
    ngen <- paste0(tgen, "_", prefnewgen)
  }
  file <- paste0(tid, "final_gen.txt")
  write(ngen, file=file, ncolumns=1)
  
  # now add mutations.
  bygen <- read.csv(file, sep="_", header=F, stringsAsFactors = F)
  
  # work out how many generations are clonal
  clonalgens <- length(which(apply(bygen, 2, function(col) length(unique(col))==1)))
  # now calculate the mutation rate based on that.
  mutspergen <- clonalmuts/clonalgens
  
  adulthood <- bygen
  
  print("Adding mutations")
  admarkers <- (unique(unlist(adulthood)))
  print(paste0(length(admarkers), " unique generation markers from adulthood"))
  
  # function to add mutations:
  # can't just use apply for this bit for all the mutations as it demands too much memory, I think.
  # use apply in batches of 500 in a loop.
  batchsize <- 500
  batchstart <- seq(1,length(admarkers), by=batchsize)
  
  adkey <- data.frame()
  for (b in batchstart) {
    print(paste0("Up to adding mutations to adult generation ", b))
    batchend <- b + (batchsize-1)
    if (batchend > length(admarkers)) {
      batchend <- length(admarkers)
    }
    tadmarkers <- admarkers[b:batchend]
    tadmuts <- sapply(tadmarkers, function(marker) add.mutns.fn(marker, mutspergen))
    tadout <- cbind(tadmarkers, tadmuts)
    adkey <- as.data.frame(rbind(adkey, tadout))
  }
  colnames(adkey) <- c("gen_marker", "mutations")
  adkey$gen_marker <- as.character(adkey$gen_marker)
  adkey$mutations <- as.character(adkey$mutations)
  
  # now replace every cell in adulthood with its partner in adkey. 
  ahood <- unlist(adulthood)
  batchstart <- seq(1,length(ahood), by=batchsize)
  adultmuts <- c()
  for (b in batchstart) {
    batchend <- b + (batchsize-1)
    if (batchend > length(ahood)) {
      batchend <- length(ahood)
    }
    tahood <- ahood[b:batchend]
    tadultmuts <- sapply(tahood, function(marker) adkey[adkey$gen_marker==marker, "mutations"])
    adultmuts <- c(adultmuts, tadultmuts)
  }
  ad <- matrix(adultmuts, nrow=nrow(adulthood), ncol=ncol(adulthood), byrow=F)
  
  wimuts <- as.data.frame(ad)
  write.table(wimuts, paste0(tid, "_with_muts.txt"), sep="_", col.names = F, row.names = F, quote=F)
  
  bymuts <- apply(wimuts, 1, function(row) paste(as.character(row), collapse="-"))
  
  # keep all the mutations for now.
  asmat1 <- sapply(bymuts, function(hap) unlist(strsplit(hap, "-")))
  # turn into matrix
  tmaxmuts <- max(sapply(asmat1, function(element) length(element)))
  bymut <- matrix(nrow=length(asmat1), ncol=tmaxmuts)

  for (asm in 1:length(asmat1)) {
    tasm <- asmat1[[asm]]
    bymut[asm,1:length(tasm)] <- tasm
  }
  
  # so at this stage have a matrix where every row is a line of descent, and every column is a mutation.
  # then need to grow to full crypt: 2,000 cells.
  # and should be able to acquire a certain number of additional mutations in every generation.
  
  ### added 2019.03.15, as was failing with fully clonal crypts
    if (ncol(bymut)==1) {
      bymut <- t(bymut)
    }
  ### end of addition
  
  
  outmat <- matrix(nrow=(cryptcells*2), ncol=(ncol(bymut)+100)) # 2019.03.15 changed the width of the matrix to reflect the number of mutations.
  outmat[1:nrow(bymut),1:ncol(bymut)] <- bymut  
  i <- 1
  while(nrow(bymut)<2000) {
    # double the cells.
    bymut <- rbind(bymut, bymut)
    # add mutations
    for (j in 1:nrow(bymut)) {
      trow <- bymut[j,]
      trow <- trow[!is.na(trow)]
      toadd <- rpois(1, lambda=mutsperdoubling)
      if (toadd>0) {
        tdoubid <- paste0("doub", i, "row", j) 
        mutstoadd <- rep(tdoubid, toadd)
        trow <- c(trow, mutstoadd)
      }
      outmat[j,1:length(trow)] <- trow
    }
    bymut <- outmat[!is.na(outmat[,1]),]
    i <- i + 1
  }
  bm <- bymut[,apply(bymut,2,function(row) !all(is.na(row)==T))]
  
  # write out the mutuations
  write.table(bm, paste0(tid, "all_muts_in_crypt.txt"), sep="\t", col.names = T, row.names = F, quote=F)
  
  # now work out the frequency of mutations
  allmuts <- unlist(bm)
  allmuts <- allmuts[!is.na(allmuts)]
  truevaf <- table(allmuts)/nrow(bm) # so a clonal mutation has truevaf=1
  # account for diploidy
  dipvaf <- truevaf/2
  
  # I also want as output the time to most recent common ancestor for the crypt.
  polyclonal_gens <- (ncol(bygen) - clonalgens)
  years_to_MRCA <- polyclonal_gens/genperyear
  
  # add contamination
  vafs_with_contam <- (truevaf*(1-stromaprop))/2
  
  # resample vafs using binomial distribution.
  # read in a file of depths. This will be different for every crypt that I use.
  deps <- as.vector(unlist(t(indata)))[3:length(as.vector(unlist(indata)))]
  deps <- deps[!is.na(deps)]
  deps <- sample(deps, length(vafs_with_contam), replace=T)
  
  simpileup <- data.frame(cbind(dipvaf, vafs_with_contam, deps), stringsAsFactors = F)
  simpileup$mtrs <- apply(simpileup, 1, function(row) {
    rbinom(n = 1, size = row["deps"], prob = row["vafs_with_contam"])
  })
  simpileup$rbinom_vafs <- simpileup$mtrs/simpileup$deps
  
  write.table(simpileup, paste0(tid, "_simulated_vafs_with_resampling.txt"), sep="\t", col.names = T, row.names=F, quote=F)
  
  # now get summary statistics: the number of mutations in each VAF bin. Sample finely as can always lump back together again.
  # do once based on >3 reads, and once on >6 reads
  bincounts_6_mtrs <- table(cut(simpileup$rbinom_vafs[simpileup$mtrs>5], breaks=seq(from=0, to=0.7, by=0.05)))
  bincounts_3_mtrs <- table(cut(simpileup$rbinom_vafs[simpileup$mtrs>=3], breaks=seq(from=0, to=0.7, by=0.05)))
  
  # write out final summary statistics.
  sumstats <- c(stem,gendays,mutsperdoubling,years_to_MRCA, as.numeric(bincounts_6_mtrs), as.numeric(bincounts_3_mtrs))
  write(sumstats, paste0(tid, '_summary_statistics.txt'), ncolumns = length(sumstats))
  
  allstats <- c(allstats, years_to_MRCA, as.numeric(bincounts_6_mtrs), as.numeric(bincounts_3_mtrs))
}
write(allstats, paste0(id, '_summary_statistics.txt'), ncolumns=length(allstats))

# remove redundant files
system(paste0('rm ', paste0(list.files(pattern='_final_gen.txt'), collapse=" ")))
system(paste0('rm ', paste0(list.files(pattern='_all_muts_in_crypt.txt'), collapse=" ")))
system(paste0('rm ', paste0(list.files(pattern='__with_muts.txt'), collapse=" ")))

# compress remaining files
tospare <- list.files(pattern ="randseed__")
allfiles <- list.files()
tozip <- allfiles[allfiles != tospare]
system(paste0('gzip ', paste0(tozip, collapse=" ")))


