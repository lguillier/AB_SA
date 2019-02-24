#'For each individual, this function defines the category as text for the multinomial logistic regression from 0,1 of the scoary trait file  
#' 
#' @param filenametrait The name of trait file produced by Scoary
#'
#' @return cat_i The category in text
#' 
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export
#'  


DefStrainCategory<-function(filenametrait)
{
data<-read.csv(filenametrait,header=TRUE,sep=",")
categories<-colnames(data)
cat_i<-c()
for (i in 1:length(data$X))
{ cat_i[i]<-categories[c(which(data[i,]==1))]}
return(cat_i)
}



