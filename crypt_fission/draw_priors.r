

## This script draws <nr_draws> crypt-fission rate values from a uniform prior distribution and writes the values to
## the output_dir


output_dir="/nfs/users/nfs_s/so11/phd/somatic_ibd_p1/ABC/normal_fission_rate/"

set.seed(42)

nr_draws <- 100000

## Crypt fission rates in the units fissions per year per crypt.
upper_limit <- 1
lower_limit <- 1/80


draws <- runif(nr_draws, min=lower_limit, max=upper_limit)

## Each simulation will have a number and a prior-draw
## the number can be used together with the biopsy number to create a unique ID
## for the simulation while still allowing us to run multiple jobs for the same
## prior value. 
draws <- data.frame(c(1:length(draws)), draws)
colnames(draws) <- c("Number", "CFR")

## Split up the file by crypt fission rate. Simulation will take different time to run
## depending on crypt fission rate and doing this will allow for more effective job
## scheduling
write.table(subset(draws, CFR<0.05), file=paste(output_dir, "crypt_fission_rate_lt_0.05.txt", sep=""), quote=F, row.names=F, col.names=F)
start=0.05
end=0.1
while(end < 1) {
  write.table(subset(draws, CFR>=start & CFR<end), file=paste(output_dir, "crypt_fission_rate_gt_", start, "_lt_",end,".txt", sep=""), 
              quote=F, row.names=F, col.names=F)
  start <- start + 0.05
  end <- end + 0.05
}


