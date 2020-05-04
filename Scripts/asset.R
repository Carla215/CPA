#set working directory
#setwd("/rdsgpfs/general/user/cs3515/home/data")

#load asset
library(ASSET)

#load in p_sumstats
p_sumstats <- readRDS("/rdsgpfs/general/user/cs3515/home/cpa/p_sumstats2.Rda")

#assign variables
snp.vars <- p_sumstats$SNP
traits.lab <- c("monop", "insulin", "IL6", "tg","vitd")
beta.hat <- subset(p_sumstats, select = c("beta.monop", "beta.insulin", "beta.IL6", "beta.tg","beta.vitd_jiang"))
sigma.hat <- subset(p_sumstats, select = c("se.monop", "se.insulin", "se.IL6", "se.tg","se.vitd_jiang"))

#calculate case and control using beta values
ncase <- 1/(sigma.hat^2)
ncntl <- 1/(sigma.hat^2)

#convert dataframes to matrices
beta.hat <- as.matrix(beta.hat)
sigma.hat <- as.matrix(sigma.hat)
ncntl <- as.matrix(ncntl)
ncase <- as.matrix(ncase)

#calculate results and save locally
hres <- h.traits(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl)
saveRDS(hres, file="/rds/general/user/cs3515/home/cpa/data/hres2.Rda")
hsum <- h.summary(hres)
saveRDS(hsum, file="/rds/general/user/cs3515/home/cpa/data/hsum2.Rda")
