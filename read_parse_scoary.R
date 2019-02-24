#' For each gene in the scoary trait file retrieve its position in roary output 
#' 
#' @param filenamescoary The name of trait file produced by Scoary
#' @param roary_Rtab The nam of gene presence absence Rtab file produced by Roary 
#' 
#' @return pos The position of Scoary genes in the roary output
#' 
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export
 

read_parse_scoary<-function(filenamescoary,roary_Rtab)
{roary_out<-read.table(roary_Rtab,header = TRUE)
data_<-read.csv(filenamescoary,header=TRUE,sep=",")
genes_<-c()
pos<-c()
for (i in 1:length(data_$Gene))
{
  genes_[i]<-sub("-\\b.*", "",data_$Gene[i])
  pos<-c(pos,which(genes_[i]==roary_out$Gene))
}
return(pos)
}
