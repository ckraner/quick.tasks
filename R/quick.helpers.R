####################### QUICK.HELP #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################

#' Formula Dissecter
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.formula=function(my.model,my.envir){
  #### Split and get pieces
  for(i in 1:length(my.model$call)){
    my.formula.grep=grep("~",my.model$call[[i]])
    if(length(my.formula.grep)>0){

      #### Get vars
      my.vars=strsplit(as.character(my.model$call[[i]]),"~")
      if(length(grep("cbind",my.vars[[2]]))>0){
        my.dep.var="my.model$model[[1]]"
      }else{
        my.dep.var=my.vars[[2]]
      }
      my.ind.vars=strsplit(my.vars[[3]],"\\+")
      my.ind.vars=lapply(my.ind.vars,trimws)

      #### Null formula
      my.null.formula=eval(parse(text=paste(my.dep.var,"~1",sep="")),envir = my.envir)

      #### Other formulas
      my.ind.vars.perm=NULL
      if(length(my.ind.vars[[1]])>2){
        for(c in 1:length(my.ind.vars[[1]])){
          my.temp.vars=my.ind.vars[[1]][-c]

          #### build the line


          if(c>1){
          my.ind.vars.perm[[2]]=rbind(my.ind.vars.perm[[2]],c(my.temp.vars,my.ind.vars[[1]][c]))
          my.ind.vars.perm[[1]]=rbind(my.ind.vars.perm[[1]],my.temp.vars)
          }else{
            my.ind.vars.perm[[2]]=as.data.frame(matrix(ncol={length(my.ind.vars)},nrow=1))
            my.ind.vars.perm[[2]]=rbind(c(my.temp.vars,my.ind.vars[[1]][c]))
            my.ind.vars.perm[[1]]=as.data.frame(matrix(ncol={length(my.ind.vars)-1},nrow=1))
            my.ind.vars.perm[[1]]=rbind(my.temp.vars)
          }
        }

        #my.ind.vars.perm[[1]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]])-1,my.ind.vars[[1]])
        #my.ind.vars.perm[[2]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]]),my.ind.vars[[1]])
      }else{
        my.ind.vars.perm[[1]]=as.data.frame(matrix(ncol=1,nrow=2))
        my.ind.vars.perm[[1]]=rbind(names(my.model$model)[3],names(my.model$model)[2])
        my.ind.vars.perm[[2]]=rbind(c(names(my.model$model)[3],names(my.model$model)[2]),
                                    c(names(my.model$model)[2],names(my.model$model)[3]))
      }
    }
  }

  #### Make formulas
  my.pre.formula.list=NULL
  my.formula.list=NULL
  #### for pieces
  for(q in 1:1){
    #### for each row
    for(j in 1:dim(my.ind.vars.perm[[q]])[1]){
      my.temp.form=paste(my.dep.var,"~ 1")
      #### for each piece in the row
      for(f in 1:dim(my.ind.vars.perm[[q]])[2]){
        my.temp.form=paste(my.temp.form,"+",my.ind.vars.perm[[q]][j,f])
      }
      my.pre.formula.list=c(as.list(my.pre.formula.list),eval(parse(text=my.temp.form),envir = my.envir))
    }
  }
  for(q in 2:2){
    #### for each row
    for(j in 1:dim(my.ind.vars.perm[[q]])[1]){
      my.temp.form=paste(my.dep.var,"~ 1")
      #### for each piece in the row
      for(f in 1:dim(my.ind.vars.perm[[q]])[2]){
        my.temp.form=paste(my.temp.form,"+",my.ind.vars.perm[[q]][j,f])
      }
      my.formula.list=c(my.formula.list,eval(parse(text=my.temp.form),envir = my.envir))
    }
  }


  my.return.list=list(my.null.formula,my.pre.formula.list,my.formula.list)

  return(my.return.list)

}




#' Formula Dissecter for later...does ALL types.
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.formula.all=function(my.model,my.envir){
  library(gtools)
  #### Split and get pieces
  for(i in 1:length(my.model$call)){
    my.formula.grep=grep("~",my.model$call[[i]])
    if(length(my.formula.grep)>0){

      #### Get vars
      my.vars=strsplit(as.character(my.model$call[[i]]),"~")
      if(length(grep("cbind",my.vars[[2]]))>0){
        my.dep.var="my.model$model[[1]]"
      }else{
      my.dep.var=my.vars[[2]]
      }
      my.ind.vars=strsplit(my.vars[[3]],"\\+")

      #### Null formula
      my.null.formula=eval(parse(text=paste(my.dep.var,"~1",sep="")),envir = my.envir)

      #### Other formulas
      my.ind.vars.perm=NULL
      if(length(my.ind.vars[[1]])>2){

        my.ind.vars.perm[[1]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]])-1,my.ind.vars[[1]])
        my.ind.vars.perm[[2]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]]),my.ind.vars[[1]])
      }else{
        my.ind.vars.perm[[1]]=as.data.frame(matrix(ncol=1,nrow=2))
        my.ind.vars.perm[[1]]=rbind(names(my.model$model)[3],names(my.model$model)[2])
        my.ind.vars.perm[[2]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]]),my.ind.vars[[1]])
      }
    }
  }

  #### Make formulas
  my.pre.formula.list=NULL
  my.formula.list=NULL
  #### for pieces
  for(q in 1:1){
    #### for each row
    for(j in 1:dim(my.ind.vars.perm[[q]])[1]){
      my.temp.form=paste(my.dep.var,"~ 1")
      #### for each piece in the row
      for(f in 1:dim(my.ind.vars.perm[[q]])[2]){
        my.temp.form=paste(my.temp.form,"+",my.ind.vars.perm[[q]][j,f])
      }
      my.pre.formula.list=c(as.list(my.pre.formula.list),eval(parse(text=my.temp.form),envir = my.envir))
    }
  }
  for(q in 2:2){
    #### for each row
    for(j in 1:dim(my.ind.vars.perm[[q]])[1]){
      my.temp.form=paste(my.dep.var,"~ 1")
      #### for each piece in the row
      for(f in 1:dim(my.ind.vars.perm[[q]])[2]){
        my.temp.form=paste(my.temp.form,"+",my.ind.vars.perm[[q]][j,f])
      }
      my.formula.list=c(my.formula.list,eval(parse(text=my.temp.form),envir = my.envir))
    }
  }


  my.return.list=list(my.null.formula,my.pre.formula.list,my.formula.list)

  return(my.return.list)

}


#' Test statistics
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

quick.m.test=function(my.matrix,test.stat){

  if(test.stat=="Wilks"){
    my.eigen=eigen(my.matrix)$values[1]
    my.inv.part=1/{1+my.eigen}
    my.wilks=my.inv.part
    for(i in 2:dim(my.matrix)[1]){
      my.eigen=eigen(my.matrix)$values[i]
      my.inv.part=1/{1+my.eigen}
      my.wilks=my.wilks*my.inv.part
    }
    return(my.wilks)
  }else if(test.stat=="Pillai" | test.stat=="Roy"){
    my.eigen=eigen(my.matrix)$values[1]
    my.theta=my.eigen/{1+my.eigen}
    if(test.stat=="Roy"){
      return(my.theta)
    }
    my.pillai=my.theta
    for(i in 2:dim(my.matrix)[1]){
      my.eigen=eigen(my.matrix)$values[i]
      my.theta=my.eigen/{1+my.eigen}
      my.pillai=my.pillai+my.theta
    }
    return(my.pillai)
  }

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
