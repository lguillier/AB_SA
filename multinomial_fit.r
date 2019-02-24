#' The function multinomilal_fit carries out the multinomial logistic regression, assesses the accuracy 
#' 
#' @param mnl_input The name trait file used as input by Scoary
#' @param percent_cross The fraction used for testing, default precent_cross=30  
#' @param nboot Number of bootstrap for assessing accuracy, default nboot=100
#' 
#' @return multinomModel_full The fitted multinomial logistic model 
#' 
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export


multinomial_fit<-function(mnl_input)
{
  set.seed(123)  

  data<-read.table(mnl_input,sep=",",header = TRUE)
data2<-subset(data,select=-1)

multinomModel_full <- multinom(Source ~ ., data=data2) # multinom Model
summary (multinomModel_full)
return(multinomModel_full)

}
