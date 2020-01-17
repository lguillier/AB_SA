#' The function multinomilal_fit carries out the multinomial logistic regression, assesses the accuracy 
#' 
#' @param mnl_input The name of the input csv file (generated with create_input_mnl.r) 
#' @param percent_cross The fraction used for testing, default precent_cross=30  
#' @param nboot Number of bootstrap for assessing accuracy, default nboot=100
#' 
#' @return accuracy_trained The accuracy of the trained model (median, 90 and 95 credibility interval). A ksdensity of accuracy is also ploted
#' @return p Plot of the accuracy
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export


MNLTrainTest<-function(mnl_input,percent_cross,nboot)
{
  set.seed(123)  
 
 # ifelse( missing(percent_cross), percent_cross<-70, percent_cross<-percent_cross ))
#ifelse( missing(nboot), nboot<-100, nboot<-nboot ))

data<-read.table(mnl_input,sep=",",header = TRUE)
data2<-subset(data,select=-1)

trainingRows <- sample(1:nrow(data), percent_cross*nrow(data))
training <- data[trainingRows, ]
training2<-subset(training,select=-1)
testing <- data[-trainingRows, ]
testing2<-subset(testing,select=-1)


multinomModel <- multinom(Source ~ ., data=training2) # multinom Model
summary (multinomModel)

predicted_scores <- predict (multinomModel, testing2, "probs") 

predicted_class <- predict (multinomModel, testing2)

table(predicted_class, testing2$Source)


mean(as.character(predicted_class) != as.character(testing2$Source))

#pred<-prediction(predicted_class,testing2$Source)
#perf<-performance(pred)

k.folds1 <- function(k) {
  folds <- createFolds(data2$Source, k = k, list = TRUE, returnTrain = TRUE)
  for (i in 1:k) {
    model <- multinom(Source ~ .,
                      data = data2[folds[[i]],])
    predictions <- predict(object = model, newdata = data2[-folds[[i]],])
    predictions2 <- predict(object = model, newdata = data2[-folds[[i]],],type="probs")
    accuracies.dt <- c(accuracies.dt,
                       confusionMatrix(predictions, data2[-folds[[i]], ]$Source)$overall[[1]])
  
    }
  return(accuracies.dt)
}

k.folds2 <- function(k) {
  folds <- createFolds(data2$Source, k = k, list = TRUE, returnTrain = TRUE)
  for (i in 1:k) {
    model <- multinom(Source ~ .,
                      data = data2[folds[[i]],])
    predictions <- predict(object = model, newdata = data2[-folds[[i]],])
    predictions2 <- predict(object = model, newdata = data2[-folds[[i]],],type="probs")
    tpr.dt<-c(tpr.dt,confusionMatrix(predictions, data2[-folds[[i]], ]$Source)$byClass[,11])
    
  }
  return(tpr.dt)
}

accuracies.dt <- c()
accuracies.dt <- k.folds1(6)
accuracies.dt

tpr.dt<-c()
tpr.dt<-k.folds2(6)
tpr.dt

v<-replicate(nboot,k.folds1(6))
accuracies.dt <- c()
for (i in 1:nboot) { accuracies.dt <- c(accuracies.dt,v[,i])}

w<-replicate(nboot,k.folds2(6))
tpr.dt <- c()
for (i in 1:nboot) { tpr.dt <- c(tpr.dt,w[,i])}

p<-densityplot(accuracies.dt,bw=0.05)
print(p)
accuracy_trained<-quantile(accuracies.dt,c(.025,0.05,0.5,0.95,0.975))
names_sources<-unique(names(tpr.dt))

getmode<-function(v)
{
  uniqv<-unique(v)
  uniqv[which.max(tabulate(match(v,uniqv))) ]
}

self_attr<-c()
for (i in 1:length(names_sources))
{pos_source<-which(names(tpr.dt)==names_sources[i])
self_attribution<-tpr.dt[pos_source]
self_attr[i]<-median(self_attribution[ which(is.finite(self_attribution))])
                   plot(ecdf(self_attribution[ which(is.finite(self_attribution))] ))}

return(list(accuracy_trained,self_attr,p))


}
