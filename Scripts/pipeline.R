#!/usr/bin/env Rscript

#Carla Smith, first version completed 13th May 2020
#Run from pipeline.sh, so this script is not called by the user directly
#Inputs are traits and cohort sizes from command-line using pipeline.sh
#Outputs are GenomicSEM & ASSET subset R objects

#load libraries
library(GenomicSEM)
library(ASSET)
library(TwoSampleMR)

#print R console to screen
options(echo=TRUE)

#get command-line arguments, sort into traits and numbers
args=commandArgs(trailingOnly=TRUE)
num_args <- length(args)/2
ts <- args[1:num_args]
ns <- args[num_args+1:length(args)]
ns <- ns[!is.na(ns)]


# GenomicSEM

#perform ldsc to get genetic covariance matrix
traits <- paste0('"',ts,".sumstats.gz",'"', collapse=", ")
sample.prev <- c(rep(NA, length(args)))
sample.prev <- paste(as.character(sample.prev), collapse=", ")
population.prev <- sample.prev
ld <- "eur_w_ld_chr/"
dwld <- "eur_w_ld_chr/"
trait.names <- paste0('"',ts,'"', collapse=",")
LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, ld=ld, wld=wld, trait.names=trait.names)
saveRDS(LDSCoutput, file="LDSCoutput.Rda")


#summary
files <- trait.names
ref <- "reference.1000G.maf.0.005.txt"
se.logit <- c(rep(T,length(args)))
se.logit <- paste(as.character(se.logit), collapse=", ")
info.filter <- 0.9
maf.filter <- 0.01
p_sumstats <- sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=trait.names,linprob=NULL,prop=NULL,N=NULL,info.filter=info.filter,maf.filter=maf.filter)
saveRDS(p_sumstats, file="p_sumstats.Rda")

#find common factors
pfactor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = p_sumstats, estimation = "DWLS", cores = NULL, parallel = FALSE, Output = NULL)

#write to outfile
saveRDS(pfactor, file="pfactor.Rda")

#common factor without using SNPs
comfactor <- commonfactor(LDSCoutput,estimation="DWLS")
saveRDS(comfactor, file="comfactor.Rda")


# ASSET

#assign variables
snp.vars <- p_sumstats$SNP
traits.lab <- trait.names
beta.hat <- paste0('"',"beta.",ts,'"', collapse=", ")
sigma.hat <- paste0('"',"se.",ts,'"', collapse=", ")
#calculate case and control using beta values
ncase <- 1/(sigma.hat^2)
ncntl <- 1/(sigma.hat^2)

#convert dataframes to matrices
beta.hat <- as.matrix(beta.hat)
sigma.hat <- as.matrix(sigma.hat)
ncntl <- as.matrix(ncntl)
ncase <- as.matrix(ncase)

#calculate results and save locally
hres <- h.traits(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl, meth.pval="B")
saveRDS(hres, file="hres.Rda")
hsum <- h.summary(hres)
saveRDS(hsum, file="hsum.Rda")


#process ASSET results

#split results into separate dataframes from one- and two-sided subset search
onesided <- hsum$Subset.1sided
twosided <- hsum$Subset.2sided


#only keep SNPs with P-value < 0.05
one_sig <- subset(onesided, Pvalue < 0.05)
two_sig <- subset(twosided, Pvalue < 0.05)

                                     
