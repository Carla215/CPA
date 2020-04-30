# load correlation matrix
rg_mat <- read.table("rg_matrix.txt", header=TRUE, sep='\t', row.names=1)

# convert to matrix
rg_mat <- as.matrix(rg_mat)

#plot
corrplot(rg_mat, method = "circle")
corrplot(rg_mat, method = "circle", type="upper")
corrplot(rg_mat, method = "number", type="upper")
