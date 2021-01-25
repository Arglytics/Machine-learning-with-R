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

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

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

