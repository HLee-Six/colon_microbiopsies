# expectation maximisation to divvy up hdp signatures into pcawg signatures.

library("lattice")
hist.cols <- rep(c("blue", "black", "red", "grey", "green", "pink"), each=16)

# i) For known signatures (e.g. when HDP calls something as “SBS 1”)
#   a.      we calculate the cosine similarity between the HDP version of the signature, and the PCAWG version of the signature.
#   b.      If the cosine similarity is below a certain threshold (0.95 seems to work well), we then apply an expectation maximisation (EM) algorithm to deconvolute this signature into its constituents. E.g. the preconditioned HDP SBS1 goes to 38.5% SBS1, 47.5% SBS5, and 15.0% SBS18. Reassuringly, 38% of the mutations that make up the HDP SBS1 are C to T at CpG, so that fits well.

# ii) For novel signatures (i.e. ones that Nicola’s algorithm calls as novel e.g. "SBS N1")
#   a.      We always apply EM to break an HDP signature down into a composite of PCAWG signatures.
#   b.      We then reconstitute the signature by adding up the PCAWG signatures in the proportion that EM gives us (e.g. IDN1 is 27% ID1, 14% ID2, 20% ID5 etc). I have tried only using signatures that contribute >10% of the mutations to avoid overfitting.
#   c.      We compute the cosine similarity of the reconstituted signature to the known signature. If that is less than a certain threshold, then we consider that the signature is truly novel. Otherwise, we break it down into its constituents and present the data in that way.

# read in the signature extraction results.
hsbs <- read.csv("Signature_category_counts.txt", sep="\t", header=T, stringsAsFactors = F)

# read in the pcawg sigs
psbs <- read.csv("PCAWG_sigProfiler_SBS_signatures.csv", header=T, stringsAsFactors = F, row.names = 1)

known <- grep("SBS", colnames(hsbs), value=T)
novel <- grep("N", colnames(hsbs), value=T)
# I think leave the residual as it is.

# calculate cosine similarity for all of them
cosmat <- data.frame(matrix(nrow=ncol(hsbs), ncol=ncol(psbs)))
rownames(cosmat) <- colnames(hsbs)
colnames(cosmat) <- colnames(psbs)

for (i in 1:nrow(cosmat)) {
  for (j in 1:ncol(cosmat)) {
    cosmat[i,j] <- cosine(x=hsbs[,rownames(cosmat)[i]], y=psbs[,colnames(cosmat)[j]])
  }
}

colnames(cosmat) <- paste0("pcawg_", colnames(cosmat))
rownames(cosmat) <- paste0("hdp_", rownames(cosmat))
s <- cosmat
write.table(cosmat, "Cosine_similarities.txt", sep="\t", col.names = T, row.names = T, quote=F)

pdf("Cosine_similarities.pdf", height=5, width=15)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot(t(s[dim(s)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()


# for known sigs, look at the cosine similarity between the hdp and the pcawg version. 
cosmat[paste0("hdp_", known),paste0("pcawg_", known)]

# only sbs 1 needs to be split.
# only split it into the good signatures.
gdsigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS25", "SBS28", "SBS30", "SBS37", "SBS40", "SBS41", "SBS43", "SBS45", "SBS49") 
sigs <- as.matrix(psbs[,gdsigs])
signatures <- t(sigs)

sample_list = c("X0", "SBS1")
mutations <- hsbs[,sample_list]
head(mutations)
signatures_names <- rownames(signatures)
signature_fraction = array(NA,dim=c(dim(signatures)[1], length(sample_list)))
rownames(signature_fraction) = signatures_names
colnames(signature_fraction) = sample_list
num_signatures = length(signatures_names)
maxiter <- 1000

for (j in 1:length(sample_list)) {
  
  mut_freqs = mutations[,j]
  #  as.numeric(table(mutations$channel[mutations$sampleID==sample_list[j]])[as.character(1:96)])
  mut_freqs[is.na(mut_freqs)] = 0
  
  # EM algowith to estimate the signature contribution
  alpha = runif(num_signatures); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  #alpha = rep(1/num_signatures,num_signatures) # Uniform start
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(num_signatures,96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * mut_freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

s <- signature_fraction
rownames(s) <- paste0("pcawg_", rownames(s))
colnames(s) <- paste0("hdp_", colnames(s))
dim(s)

pdf("EM_break_up_known_sigs.pdf", height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((s[dim(s)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(s, "hdp_known_sigs_broken_down_into_pcawg_gd_sigs.txt", sep="\t", col.names=T, row.names = T, quote=F)

constit <- rownames(signature_fraction[signature_fraction[,"SBS1"]>0.1,])

# now re-do the EM just on these three sigs.
gdsigs <- constit
sigs <- as.matrix(psbs[,gdsigs])
signatures <- t(sigs)

sample_list = c("X0", "SBS1")
mutations <- hsbs[,sample_list]
head(mutations)
signatures_names <- rownames(signatures)
signature_fraction = array(NA,dim=c(dim(signatures)[1], length(sample_list)))
rownames(signature_fraction) = signatures_names
colnames(signature_fraction) = sample_list
num_signatures = length(signatures_names)
maxiter <- 1000

for (j in 1:length(sample_list)) {
  
  mut_freqs = mutations[,j]
  #  as.numeric(table(mutations$channel[mutations$sampleID==sample_list[j]])[as.character(1:96)])
  mut_freqs[is.na(mut_freqs)] = 0
  
  # EM algowith to estimate the signature contribution
  alpha = runif(num_signatures); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  #alpha = rep(1/num_signatures,num_signatures) # Uniform start
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(num_signatures,96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * mut_freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

s <- signature_fraction
rownames(s) <- paste0("pcawg_", rownames(s))
colnames(s) <- paste0("hdp_", colnames(s))
dim(s)

pdf("EM_break_up_SBS1_into_pcawg_1-5-18.pdf", height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((s[dim(s)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(s, "hdp_SBS_broken_down_into_pcawg_1-5-18.txt", sep="\t", col.names=T, row.names = T, quote=F)
s
# HDP SBS 1 is: 38.5% psbs1, 47.5% psbs5, 14.0% psbs18.

# reconstitute HDP SBS1 using those PCAWG sigs. Calculate the cosine similarity between the original and the reconstituted pcawg sig.
constit
reconsbs <- rep(0, 96)
for (ct in constit) {
  reconsbs <- reconsbs + (psbs[,ct]*s[paste0("pcawg_", ct),"hdp_SBS1"])
}
reconsbs
sum(reconsbs)
cosine(x=reconsbs, y=hsbs$SBS1) # 0.99

# plot the original
# their broken down signatures
# the reconstituted sigs
pdf("HDP_SBS1_reconstitution.pdf")
par(mfrow=c(5,1))
par(mar=c(1,2,4,1))
barplot(hsbs$SBS1, col=hist.cols, main="HDP SBS1")
barplot(reconsbs, col=hist.cols, main=paste0("Reconstituted SBS1, cosine similarity to original: ", round(cosine(x=reconsbs, y=hsbs$SBS1), digits=2)))
for (ct in constit) {
  barplot(psbs[,ct], col=hist.cols, main=paste0("PCAWG ", ct, " accounts for ", round(s[paste0("pcawg_", ct),"hdp_SBS1"], digits=2)))
}
dev.off()



####################
####################
# now break down all novel signatures, rebuild them back up, and see how well each is explained.
# only split it into the good signatures.
gdsigs <- c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS25", "SBS28", "SBS30", "SBS37", "SBS40", "SBS41", "SBS43", "SBS45", "SBS49") 
sigs <- as.matrix(psbs[,gdsigs])
signatures <- t(sigs)

sample_list = novel
mutations <- hsbs[,sample_list]
head(mutations)
signatures_names <- rownames(signatures)
signature_fraction = array(NA,dim=c(dim(signatures)[1], length(sample_list)))
rownames(signature_fraction) = signatures_names
colnames(signature_fraction) = sample_list
num_signatures = length(signatures_names)
maxiter <- 10000

for (j in 1:length(sample_list)) {
  
  mut_freqs = mutations[,j]
  #  as.numeric(table(mutations$channel[mutations$sampleID==sample_list[j]])[as.character(1:96)])
  mut_freqs[is.na(mut_freqs)] = 0
  
  # EM algowith to estimate the signature contribution
  alpha = runif(num_signatures); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  #alpha = rep(1/num_signatures,num_signatures) # Uniform start
  for (iter in 1:maxiter) {
    contr = t(array(alpha,dim=c(num_signatures,96))) * t(signatures)
    probs = contr/array(rowSums(contr),dim=dim(contr))
    probs = probs * mut_freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  
  # Saving the signature contributions for the sample
  print(j/length(sample_list))
  signature_fraction[,j] = alpha
}

s <- signature_fraction
rownames(s) <- paste0("pcawg_", rownames(s))
colnames(s) <- paste0("hdp_", colnames(s))
dim(s)

pdf("EM_break_up_novel_sigs.pdf", height=5, width=10)
color.palette = colorRampPalette(c("white", "orange", "purple"))
levelplot((s[dim(s)[1]:1,]),col.regions=color.palette, aspect="fill", scales=list(x=list(rot=90)))
dev.off()

write.table(s, "hdp_novel_sigs_broken_down_into_pcawg_gd_sigs.txt", sep="\t", col.names=T, row.names = T, quote=F)


pdf("Breaking_down_novel_sigs_into_pcawg_sigs_by_nmf.pdf")
for (ss in colnames(s)) {
  print(ss)
  tcontrib <- names(which(s[,ss]>0.1))
  pcontribsigs <- sapply(strsplit(tcontrib, "_"), "[[", 2)
  thdpsig <- sapply(strsplit(ss, "_"), "[[", 2)
  
  if (length(pcontribsigs)<=1) {
    tcos <- s[tcontrib,ss]
    
    par(mfrow=c(3,1))
    par(mar=c(1,2,4,1))
    barplot(hsbs[,thdpsig], col=hist.cols, main=paste0("Original hdp ", thdpsig))
    barplot(psbs[,pcontribsigs], col=hist.cols, main=paste0("Reconstituted, cosine similarity to original: ", round(tcos, digits=3)))
    for (ct in pcontribsigs) {
      barplot(psbs[,ct], col=hist.cols, main=paste0("PCAWG ", ct, " accounts for 100%"))
    }
    
    next
  }
  
  # need to re-do the NMF just with the selected sigs.
  sigs <- as.matrix(psbs[,pcontribsigs])
  signatures <- t(sigs)
  
  sample_list = novel
  mutations <- hsbs[,sample_list]
  signatures_names <- pcontribsigs
  signature_fraction = array(NA,dim=c(dim(signatures)[1], length(sample_list)))
  rownames(signature_fraction) = signatures_names
  colnames(signature_fraction) = sample_list
  num_signatures = length(signatures_names)
  maxiter <- 10000
  
  for (j in 1:length(sample_list)) {
    mut_freqs = mutations[,j]
    #  as.numeric(table(mutations$channel[mutations$sampleID==sample_list[j]])[as.character(1:96)])
    mut_freqs[is.na(mut_freqs)] = 0
    # EM algowith to estimate the signature contribution
    alpha = runif(num_signatures); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
    #alpha = rep(1/num_signatures,num_signatures) # Uniform start
    for (iter in 1:maxiter) {
      contr = t(array(alpha,dim=c(num_signatures,96))) * t(signatures)
      probs = contr/array(rowSums(contr),dim=dim(contr))
      probs = probs * mut_freqs
      old_alpha = alpha
      alpha = colSums(probs)/sum(probs)
      if (sum(abs(alpha-old_alpha))<1e-5) {
        break
      }
    }
    # Saving the signature contributions for the sample
    print(j/length(sample_list))
    signature_fraction[,j] = alpha
  }
  
  ts <- signature_fraction[,sapply(strsplit(ss, "_"), "[[", 2)]
  trecon <- rep(0, 96)
  for (psig in pcontribsigs) {
    trecon <- trecon + (psbs[,psig]*as.numeric(ts[psig]))
  }
  
  thdpsig <- sapply(strsplit(ss, "_"), "[[", 2)
  tcos <- cosine(trecon, hsbs[,thdpsig])
  
  # plot the original
  # their broken down signatures
  # the reconstituted sigs
  par(mfrow=c(length(pcontribsigs)+2,1))
  par(mar=c(1,2,4,1))
  barplot(hsbs[,thdpsig], col=hist.cols, main=paste0("Original hdp ", thdpsig))
  barplot(trecon, col=hist.cols, main=paste0("Reconstituted, cosine similarity to original: ", round(tcos, digits=3)))
  for (ct in pcontribsigs) {
    barplot(psbs[,ct], col=hist.cols, main=paste0("PCAWG ", ct, " accounts for ", round(s[paste0("pcawg_", ct),ss], digits=3)))
  }
}
dev.off()

### I don't think that there is any evidence that we should be breaking down the novel ones further.
# just SBS1.