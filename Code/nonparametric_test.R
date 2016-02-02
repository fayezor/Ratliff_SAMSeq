library("samr")
# Let's say I want to edit this. 

###---------- Getting SULT2B1b expression ----------###

d.full <- read.csv("GeneCountMatrix.AllBatches_withannot_.csv")
rownames(d.full) <- d.full[,1]
d.full <- d.full[,-c(1,2)]
samples.full <- names(d.full)

# Remove the control and repeat samples (for now)
ctrl.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="trl", arr.ind=TRUE)
repeat.idx <- which(substr(samples.full, nchar(samples.full)-3+1, nchar(samples.full))=="eat", arr.ind=TRUE)
remove.idx <- c(ctrl.idx, repeat.idx)
ctrls <- d.full[, ctrl.idx]
repeats <- d.full[, repeat.idx]

d <- d.full[, -remove.idx] # 36135 genes, 399 samples

# Define groups by SULT2B1b
sult2b1b <- as.numeric(d[rownames(d)=="ENSG00000088002",])
sult2b1b.zero <- sult2b1b==0
sult2b1b.grp <- ifelse(sult2b1b.zero==TRUE, "Zero", "Nonzero")
table(sult2b1b.grp) # 352 zero vs. 47 nonzero

###---------- Getting filtered dataset ----------###

setwd("/Users/fayezheng/Dropbox/Faye/GitHub/Ratliff_Prostate")

# Read in the dataset, already filtered to remove low-count genes, including SULT2B1b
d <- read.table("prostatedata.txt")

samfit <- SAMseq(x=d, y=sult2b1b, resp.type="Quantitative", geneid=rownames(d))
samfit <- res

# Examine significant gene list
#print(samfit)

# Plot results
#plot(samfit)
