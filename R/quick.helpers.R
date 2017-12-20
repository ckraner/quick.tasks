####################### QUICK.HELP #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################

#' Trace of a matrix
#'
#' Helper function to quickly get list of column names
#' based on partial piece.
#'
#' @param myDF Dataframe
#' @param partial partial search as string
#' @return NULL
#' @keywords Explore
#' @examples
#' quick.search(myDF,"CAS")

quick.tr=function(my.matrix){
  my.trace=NULL
  my.trace=my.matrix[1,1]
  for(i in 2:dim(my.matrix)[1]){
    my.trace=my.trace+my.matrix[i,i]
  }
  return(my.trace)
}

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

#' Get model type
#'
#' Helper function to see SPSS label
#'
#' @param variable Variable from data set
#' @return String of label
#' @keywords Explore
#' @export
#' @examples
#' quick.lavaan(myfit)

quick.type=function(model){
  my.call = as.character(my.model$call)
  my.split.call = strsplit(my.call, "\\(")
  my.reg.type2 = my.split.call[[1]][1]
  if (my.reg.type2 == "lm" | my.reg.type2 == "stats::lm") {
    my.reg.type = "lm"
    #ab.len=15
  } else if (my.reg.type2 == "glm" | my.reg.type2 == "stats::glm") {
    my.reg.type = "glm"
    #ab.len=15
  } else if (my.reg.type2 == "manova" |
             my.reg.type2 == "stats::manova") {
    my.reg.type = "manova"
    #ab.len=15
  } else if (my.reg.type2 == "clm" |
             my.reg.type2 == "ordinal::clm") {
    my.reg.type = "ord"
    # ab.len=30
    # library(ordinal)
  } else{
    stop("Type not supported")
  }
  return(my.reg.type)
}
