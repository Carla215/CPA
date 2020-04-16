library(GenomicSEM)

#get sumstats into correct format
munge(c("monop_locke","insulin_prokopenkol","IL6_ahola","tricglycerides_teslovich","fibrinogen_prins"),"w_hm3.noMHC.snplist",trait.names=c("monop","insulin","IL6","tg","fibrinogen"),c(131544,10831,95454,8293,96598,9762), info.filter = 0.9, maf.filter = 0.01)

#perform ldsc
traits <- c("monop.sumstats.gz", "insulin.sumstats.gz", "IL6.sumstats.gz","tg.sumstats.gz", "fibrinogen.sumstats.gz")
sample.prev <- c(NA,NA,NA,NA,NA)
population.prev <- c(NA,NA,NA,NA,NA)
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
trait.names<-c("monop","insulin","IL6","tg","fibrinogen")
LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, ld=ld, wld=wld, trait.names=trait.names)
saveRDS(LDSCoutput, file="/rds/general/user/cs3515/home/data/LDSCoutput.Rda")

#summary
files=c("monop_locke", "insulin_prokopenkol", "IL6_ahola", "triglycerides_teslovich", "fibrinogen_prins")
ref= "reference.1000G.maf.0.005.txt"
trait.names=c("monop","insulin","IL6","tg","fibrinogen")
se.logit=c(T,T,T,T,T)
info.filter=0.9
maf.filter=0.01
p_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=NULL,linprob=NULL,prop=NULL,N=NULL,info.filter=info.filter,maf.filter=maf.filter)
saveRDS(p_sumstats, file="/rds/general/user/cs3515/home/data/p_sumstats.Rda")


pfactor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = FALSE, Output = NULL)

#write to outfile
saveRDS(pfactor, file="/rds/general/user/cs3515/home/data/pfactor.Rda")
