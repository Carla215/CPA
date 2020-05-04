#load hyprcoloc
library(hyprcoloc)


#load in p_sumstats
p_sumstats <- readRDS("p_sumstats2.Rda")

#assign variables
rsid <- p_sumstats$SNP
traits <- c("monop", "insulin", "IL6", "tg","vitd")
betas <- subset(p_sumstats, select = c("beta.monop", "beta.insulin", "beta.IL6", "beta.tg","beta.vitd_jiang"))
ses <- subset(p_sumstats, select = c("se.monop", "se.insulin", "se.IL6", "se.tg","se.vitd_jiang"))

#convert dataframes to matrices
betas <- as.matrix(betas)
ses <- as.matrix(ses)


#run hyprcoloc
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
#saveRDS(hcl, file="/rds/general/user/cs3515/home/data/hcl.Rda")


#ALTERNATIVE - with SNP subset

#pick one subset for analysis
mono_il6_tg_ins <- subsets_one$`monop,insulin,IL6,tg`


#find matching SNPs from this subset and get betas/ses from p_sumstats
psumstats_miti <- subset(p_sumstats, SNP %in% mono_il6_tg_ins$SNP)


#assign variables
rsid <- psumstats_miti$SNP
traits <- c("monop", "insulin", "IL6", "tg", "vitd")
betas <- subset(psumstats_miti, select = c("beta.monop", "beta.insulin", "beta.IL6", "beta.tg","beta.vitd_jiang"))
ses <- subset(psumstats_miti, select = c("se.monop", "se.insulin", "se.IL6", "se.tg","se.vitd_jiang"))

#convert dataframes to matrices
betas <- as.matrix(betas)
ses <- as.matrix(ses)


#run hyprcoloc
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)
#saveRDS(hcl, file="/rds/general/user/cs3515/home/data/hcl.Rda")