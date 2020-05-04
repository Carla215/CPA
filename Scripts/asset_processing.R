#load results from ASSET h.traits
asset_res <- readRDS("hsum.Rda")

#split results into separate dataframes from one- and two-sided subset search
onesided <- asset_res$Subset.1sided
twosided <- asset_res$Subset.2sided
rm(asset_res)

#only keep SNPs with P-value < 5e-8
one_sig <- subset(onesided, Pvalue < 5e-8)
two_sig <- subset(twosided, Pvalue < 5e-8)

#for two-sided, this goes for P-values in both directions too
two_really_sig <- subset(two_sig, Pvalue.1 < 5e-8)
two_really_sig <- subset(two_really_sig, Pvalue.2 < 5e-8)


#split df into list of dataframes based on unique subsets
two_sig[is.na(x = two_sig)] <- "" #remove na
subsets_one <- split(one_sig, f = one_sig$Pheno)
subsets_two <- split(two_sig, f = list(two_sig$Pheno.1, two_sig$Pheno.2))

#remove empty subsets
subsets_one <- subsets_one[sapply(subsets_one, nrow)>0]
subsets_two <- subsets_two[sapply(subsets_two, nrow)>0]

#extract a subset of interest, find corresponding SNPs in original dataframe



