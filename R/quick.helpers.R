####################### QUICK.HELP #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################


#' Partial Name Search of Columns
#'
#' Helper function to quickly get list of column names
#' based on partial piece.
#'
#' @param myDF Dataframe
#' @param partial partial search as string
#' @return NULL
#' @keywords Explore
#' @export
#' @examples
#' quick.search(myDF,"CAS")

quick.search=function(myDF,partial){
  th1=grep(partial,colnames(myDF))
  print(colnames(myDF)[th1])
}


#' View SPSS Label
#'
#' Helper function to see SPSS label
#'
#' @param variable Variable from data set
#' @return String of label
#' @keywords Explore
#' @export
#' @examples
#' quick.lavaan(myfit)

quick.label=function(variable){
  print(attributes(variable)$label)
}
