#' The function prepares the input file for multinomial logistic regression
#' 
#' @param traitfile The name trait file used as input by Scoary
#' @param roary_Rtab The name of gene presence absence Rtab file produced by Roary 
#' @param max_genes The maximum of enriched genes to include per host group
#' 
#' @return The input csv file for fitting multinomial logistic model
#' @return the input csv file for prediction (of strains to attribute)
#' 
#'  @author Laurent Guillier, \email {guillier.laurent@gmail.com}
#'  
#'  @export


create_input_mnl<-function(traitfile,roary_Rtab,max_genes)
{


###### 1. Parse information from scoary trait file: reservoir names and association of a strain to its source

traits<-read.table(traitfile,header = TRUE,sep=",")
trait_names<-colnames(traits)
trait_names<-sort(trait_names[2:length(trait_names)])

categories<-def_strain_category(traitfile)
id_categories<-data.frame(id=traits$X,cat=categories) # will be used in step 3.

###### 2. Gather information from scoary enriched genes
# Load of output Rtab from roary
roary_Rtab<-c(roary_Rtab)

#List of files in the directory
files<-list.files() 

# List of csv files= output files of enriched genes from scoary
pos_scoaryfiles<-grep(pattern = "results.csv",files) 
filenamescoary<-files[pos_scoaryfiles]

# Extract from roary Rtab the list of all enriched genes (whatever the source) : create a list object
all_enriched<-mapply(read_parse_scoary,filenamescoary,roary_Rtab)
names(all_enriched)<-trait_names 


####### 3. Create input file for fitting multinomial logistic model, max_genes input is used to define the initial genes of interest

# Merge the genes in the list object in a data frame 
input<-c()
for (i in 1:length(all_enriched)) { potential<-all_enriched[[i]]
  chosen<-potential[1:min(max_genes,length(all_enriched[[i]]))]
  input<-c(input,chosen)}

# Keep one record by gene (a gene can be enriched in two or more source)
input<-unique(input)

# Preparation of gene presence/absence matrix
roary_out<-read.table(roary_Rtab,header = TRUE)
mnl_input<-roary_out[input,]
t_mnl <- transpose(mnl_input)
colnames(t_mnl) <- mnl_input$Gene
rownames(t_mnl) <- colnames(mnl_input)
t_mnl<-t_mnl[-c(1),]
row_mnl<-rownames(t_mnl)
pos<-c()
cat<-c()
for (i in (1:length(row_mnl)))
{if (isTRUE(which(row_mnl[i]==id_categories$id)>0))
     {val<-i
pos<-c(pos,val)
val-c()}}

t_mnl_input<-t_mnl[pos,]

# Addition of source to each row (strains)
Source<-id_categories$cat
t_mnl_input<-cbind(t_mnl_input,Source)

# Creation of the input for fitting
write.csv(t_mnl_input,file="mnl_input_0.csv")


####### 4. Create input file for prediction 

t_mnl_predict<-t_mnl[-pos,]

write.csv2(t_mnl_predict,file="predict_sporadic.csv")

}
