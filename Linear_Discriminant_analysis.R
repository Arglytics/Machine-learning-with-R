#LDA model to tell based on the reagents we have which group each should be in

require(MASS)
require(ggplot2)
 
#reading in the data iris
data("iris")
my.data <- iris
head(my.data)


#ploting mydata to see how the groups are..


## create the lda model
model <- lda(formula = Species ~ ., data = my.data) #note:the response(Species must be a categorical variable)
model    #calls up the model

#plot( my.data[ , c(2,3,4,5) ],col=data[ ,1 ] 

#Using a predict function, a very generic function,and we give it an object on which the predict will work on,
#so our object in this case is model and we asign this predict function to the data.lda.values to call the coordinates. 

data.lda.values <- predict(model)
data.lda.values$x   #it shows us the x coordinates
#mydatastruct<- predict(model)
#mydatastruct$x  #calling the object x of two columns LD1 and LD2

## create a dataframe that has all the info we need to draw a graph,[,1] means column one of the dataset
plot.data <- data.frame(X=data.lda.values$x[,1], Y=data.lda.values$x[,2], Species=my.data$Species)

head(plot.data)

table(X=data.lda.values$x[,1]) #prints column x

## draw a graph using ggplot2
p <- ggplot(data=plot.data, aes(x=X, y=Y)) +
  geom_point(aes(color=Species)) +
  theme_bw()

p


