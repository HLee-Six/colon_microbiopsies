
## This script uses the abc package to approximate the posterior distribution of the crypt
## fission rate parameter. 
## It expects event_counts for each simulation to be available. 

## The script was written by Sigurgeir Olafsson for the paper
## 'The landscape of somatic mutation in normal colorectal epithelial cells' by Henry Lee-Six et al.
## Please direct any questions or comments you may have to so11@sanger.ac.uk


library(pracma)
library(abc)
library(cowplot)

early_life_limit <- "4"
batches <- c("lt_0.0125","lt_0.05", "gt_0.05_lt_0.1", "gt_0.1_lt_0.15","gt_0.15_lt_0.2","gt_0.2_lt_0.25")
observed <- read.table(paste("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/observed_event_counts_excl", early_life_limit, "yr.txt", sep=""), h=T)
observed_m <- as.matrix(observed[, c(1:8)])
obs_sumStats <- apply(observed_m, 2, sum)



sim_cfr <- as.numeric()
sim_id <- as.numeric()
sim_sumStats <- data.frame()
j <- 1
for(batch in batches) {
  sim_list <- read.table(paste("/Users/so11/phd/so11_nfs/somatic_ibd_p1/ABC/normal_fission_rate/crypt_fission_rate_", batch, ".txt", sep=""))
  colnames(sim_list) <- c("ID", "CFR")
  print(batch)

  for(i in 1:nrow(sim_list)) {
    #sim_path <- paste("/Users/so11/phd/so11_lustre/somatic_ibd_p1/ABC/normal_fission_rate/event_times/gt_0.05_lt_0.1/event_counts_", sim_list2$ID[i], ".txt", sep="")
    sim_path <- paste("/Users/so11/phd/so11_abc/event_times/", batch, "/excl", early_life_limit, "yr/event_counts_", sim_list$ID[i], "_excl", early_life_limit, "yr.txt", sep="")
    if(!file.exists(sim_path)) {
      sim_cfr[j] <- NA
      sim_id[j] <- NA
      next
    }
    sim_df <- read.table(sim_path, comment.char="", h=T, nrows=104,
                         colClasses=c("numeric", "numeric","numeric", "numeric","numeric", "numeric","numeric", "numeric"))
    
    ## If some simulations were for some reason unfinished
    if(nrow(sim_df)<102) {
      sim_cfr[j] <- NA
      sim_id[j] <- NA
      print(paste("Simulation", sim_list$ID, "from batch:", batch, "had too few grids", nrow(sim_df)))
      next
    }
    sim_mat <- as.matrix(sim_df[, c(1:8)])

    sim_cfr[j] <- sim_list$CFR[i]
    sim_id[j] <- sim_list$ID[i]
    sim_sumStats <- rbind(sim_sumStats, apply(sim_mat, 2, sum))
    j <- j+1
  }
}


missing <- which(is.na(sim_cfr))
sim_cfr <- sim_cfr[-missing]
sim_id <- sim_id[-missing]
sim_sumStats <- sim_sumStats[-missing,]

test_reject <- abc(target=obs_sumStats, param=sim_cfr, sumstat=sim_sumStats, tol=0.05, method="rejection")
test_loclinear <- abc(target=obs_sumStats, param=sim_cfr, sumstat=sim_sumStats, tol=0.05, method="loclinear")
test_ridge <- abc(target=obs_sumStats, param=sim_cfr, sumstat=sim_sumStats, tol=0.05, method="ridge")
test_neural <- abc(target=obs_sumStats, param=sim_cfr, sumstat=sim_sumStats, tol=0.05, method="neuralnet")

plot(test_loclinear, param=sim_cfr)
plot(test_ridge, param=sim_cfr)
plot(test_neural, param=sim_cfr)

par(mfrow=c(2,2), mar=c(4, 4, 2, 2))
hist(test_reject$unadj.values, xlim=c(1/80, 0.2))
hist(test_loclinear$adj.values, xlim=c(1/80, 0.2))
hist(test_ridge$adj.values, xlim=c(1/80, 0.2))
hist(test_neural$adj.values, xlim=c(1/80, 0.2))


p1 <- ggplot(data.frame(sim_cfr), aes(x=sim_cfr)) + geom_histogram() + xlim(0,0.25) + labs(x="", y="Frequency")
p2 <- ggplot(data.frame(test_reject$unadj.values), aes(x=test_reject.unadj.values)) + geom_histogram() + xlim(0,0.25) + labs(x="", y="Frequency")
p3 <- ggplot(data.frame(test_neural$adj.values), aes(x=test_neural.adj.values)) + geom_histogram() + xlim(0,0.25) + labs(x="", y="Frequency")

plot_grid(p1, p2, p3, nrow=3, labels="AUTO")

x <- sim_sumStats
x$cfr <- sim_cfr
x$id <- sim_id

y <- data.frame(Reject=test_reject$unadj.values, LocLin=test_loclinear$adj.values, Ridge=test_ridge$adj.values, NeuralN=test_neural$adj.values)

write.table(x, file="/Users/so11/phd/normal_crypt_fission/sim_sumStats_excl4yr_July_Upd.txt", row.names = F, quote = F)
write.table(y, file="/Users/so11/phd/normal_crypt_fission/cfr_accepted_sims_excl4yr_July_Upd.txt", row.names=F, quote = F)
write.table(sim_cfr, file="/Users/so11/phd/normal_crypt_fission/prior_used_July_Upd.txt", row.names = F, quote=F, col.names = F)


