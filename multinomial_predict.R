#' The function multinomilal_predict returns the membership probabilities of strains with unknown source according to presence/absence matrix for enriched genes and the multinomial logistic model trained and tested on sources
#' 
#' @param sporadic_pres_abs The name of the csv file caracteryzing the presence absence matrix of enriched genes in the strains to attribute
#' @param fitted_mnl The multonomial fitted model (output from multinomial_fit.r
#' 
#' 
#' @return The parameter 
#' @return the full fitted model
#' 
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export


multinomial_predict<-function(sporadic_pres_abs,fitted_mnl)
{
data_pred<-read.table(sporadic_pres_abs,sep=";",header = TRUE)
data_pred2<-subset(data_pred,select=-1)
predicted_sporadic <- predict (fitted_mnl, data_pred2, "probs") 
return(predicted_sporadic)
}
