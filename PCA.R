## In this example, the data is in a matrix called
## data.matrix
## columns are individual samples (i.e. cells)
## rows are measurements taken for all the samples (i.e. genes)
## A matrix of data with 10 samples where we measured 100 genes with each sample. 

data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))
#for fake genes and fake read counts
rownames(data.matrix) <- paste("gene", 1:100, sep="")
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}

head(data.matrix)  #Displaying the first 6 rows.
dim(data.matrix)  #displaying the dimension of the matrix.

#Calling the prcomp() to do the PCA.
pca <- prcomp(t(data.matrix), scale=TRUE) 

#Goal is to draw a graph to show how the samples are related or not related to eachother.
##NOTE: By default, prcomp( expects the samples to be rows and the genes to be columns.
##Since the samples in our data matrix are columns and the genes(variables) are rows, we have to transpose
##the matrix using t() function.
##If we dont transpose the matrix, we will ultimately get a graph that shows how the genes are related
#to eachother.

#prcomp() returns three things i.e x, sdev, rotation.
#x contains the principal componenets for drawing a graph
##NOTE: the first PC accounts for the most variation in the original data(the gene expression across all 10 samples) and the 2nd PC ACCOUNTS
#for the second most variation and so on.

#Remember since there are 10 samples, there are 10 PCs.

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])
#Graph showing 5 samples on one side of the graph and 5 others on the other side.

#To get a sense of how meaningful these clusters are, let's see how much variation in the original data 
#PC1 accounts for

#To do this, we use the square of sdev(standard deviation) to calculate how much variation
#in the original data each principle componenet accounts for
#So we calculate the percentages using a bar plot.

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Screen Plot", xlab="Principal Component", ylab="Percent Variation")

#According to bar graph, PC1 accounts for almost all the variation in tha data

## now make a fancy looking plot that shows the PCs and the variation:



install.packages("ggplot2")

library(ggplot2)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data
##We obtain a dataframe.we have one row per sample.Each row has a sample ID and X/Y coordinates to that sample.

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

#######
##
## NOTE: Everything that follow is just bonus stuff.
## It simply demonstrates how to get the same
## results using "svd()" (Singular Value Decomposition) or using "eigen()"
## (Eigen Decomposition).
##
#######

############################################
##
## Now let's do the same thing with svd()
##
## svd() returns three things
## v = the "rotation" that prcomp() returns, this is a matrix of eigenvectors
##     in other words, a matrix of loading scores
## u = this is similar to the "x" that prcomp() returns. In other words,
##     sum(the rotation * the original data), but compressed to the unit vector
##     You can spread it out by multiplying by "d"
## d = this is similar to the "sdev" value that prcomp() returns (and thus
##     related to the eigen values), but not
##     scaled by sample size in an unbiased way (ie. 1/(n-1)).
##     For prcomp(), sdev = sqrt(var) = sqrt(ss(fit)/(n-1))
##     For svd(), d = sqrt(ss(fit))
##
############################################

svd.stuff <- svd(scale(t(data.matrix), center=TRUE))

## calculate the PCs
svd.data <- data.frame(Sample=colnames(data.matrix),
                       X=(svd.stuff$u[,1] * svd.stuff$d[1]),
                       Y=(svd.stuff$u[,2] * svd.stuff$d[2]))
svd.data

## alternatively, we could compute the PCs with the eigen vectors and the
## original data
svd.pcs <- t(t(svd.stuff$v) %*% t(scale(t(data.matrix), center=TRUE)))
svd.pcs[,1:2] ## the first to principal components

svd.df <- ncol(data.matrix) - 1
svd.var <- svd.stuff$d^2 / svd.df
svd.var.per <- round(svd.var/sum(svd.var)*100, 1)

ggplot(data=svd.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", svd.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", svd.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("svd(scale(t(data.matrix), center=TRUE)")

############################################
##
## Now let's do the same thing with eigen()
##
## eigen() returns two things...
## vectors = eigen vectors (vectors of loading scores)
##           NOTE: pcs = sum(loading scores * values for sample)
## values = eigen values
##
############################################
cov.mat <- cov(scale(t(data.matrix), center=TRUE))
dim(cov.mat)

## since the covariance matrix is symmetric, we can tell eigen() to just
## work on the lower triangle with "symmetric=TRUE"
eigen.stuff <- eigen(cov.mat, symmetric=TRUE)
dim(eigen.stuff$vectors)
head(eigen.stuff$vectors[,1:2])

eigen.pcs <- t(t(eigen.stuff$vectors) %*% t(scale(t(data.matrix), center=TRUE)))
eigen.pcs[,1:2]

eigen.data <- data.frame(Sample=rownames(eigen.pcs),
                         X=(-1 * eigen.pcs[,1]), ## eigen() flips the X-axis in this case, so we flip it back
                         Y=eigen.pcs[,2]) ## X axis will be PC1, Y axis will be PC2
eigen.data

eigen.var.per <- round(eigen.stuff$values/sum(eigen.stuff$values)*100, 1)

ggplot(data=eigen.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", eigen.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", eigen.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("eigen on cov(t(data.matrix))")
© 2021 GitHub, Inc.

