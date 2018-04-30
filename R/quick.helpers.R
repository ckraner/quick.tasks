####################### QUICK.HELP #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################



#' Check equality
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.eq.check=function(val.1,val.2){
  if(val.2=="&nbsp;"){
    val.2=NA
  }
  if(is.na(val.1) | is.na(val.2)){
    return(is.na(val.1) & is.na(val.2))
  }else{
    return(val.1==val.2)
  }
}


#' P-val Stuff
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.p.val=function(my.table3,p.val.row){
  for(i in 1:dim(my.table3)[1]){
    if(!is.na(my.table3[i,p.val.row])){
      my.table3[i,p.val.row]=round(as.numeric(my.table3[i,p.val.row]),digits=3)
      if(my.table3[i,p.val.row]==0){my.table3[i,p.val.row]="<.001"}
      if(my.table3[i,p.val.row]==1){my.table3[i,p.val.row]=">.999"}
    }
  }
  return(my.table3)
}

#' Formula Dissecter
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#'
#' @export
#'
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


#' Create latent variables
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param SSCP.list Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.latent=function(SSCP.list,marginality){
  library(purrr)
  temp.latent.SSCP.eigen=purrr::map(SSCP.list[[1]],quick.latent.reduce)
  if(!marginality){
    temp.latent.SSCP.error=quick.latent.reduce(SSCP.list[[2]])
  }else{
    temp.latent.SSCP.error=quick.latent.reduce(SSCP.list[[2]][[1]])
  }
  temp.latent=list(temp.latent.SSCP.eigen,temp.latent.SSCP.error)
  return(temp.latent)
}




#' Create latent variables
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.part Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.latent.reduce=function(my.part){
  temp.latent.SSCP.eigen=eigen(my.part)
  temp.latent.SSCP.values=as.matrix(temp.latent.SSCP.eigen$values)
  temp.latent=temp.latent.SSCP.eigen$vectors%*%temp.latent.SSCP.values
  return(temp.latent)
}



#' Put together SS from car::Anova package
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param my.model Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.SSP.matrix=function(my.summary,marginality){
  if(!marginality){
    library(car)
    car.sum=car::Anova(my.summary,type=3)
    sum.SSCP=list(car.sum$SSP,car.sum$SSPE)
  }else{
    sum.sum=summary(my.summary,intercept=T)
    sum.SSCP=list(sum.sum$SS[-length(sum.sum$SS)],sum.sum$SS[length(sum.sum$SS)])
  }
  return(sum.SSCP)
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


#' Quickly copy attributes from original data frame
#'
#' When you make a subset of data frames and do some other operations,
#' R unfortunately loses some of the attributes. This will help get them
#' back for use in label.explor.r
#'
#' @param my.new.df New data frame, normally subset of original
#' @param my.orig.df Original data frame
#'
#' @return New data frame with all relevant attributes from original df copied over
#' @export

quick.attribs=function(my.new.df,my.orig.df){
  for(i in 1:dim(my.new.df)[2]){
    if(length(grep(paste("^",tolower(colnames(my.new.df)[i]),"$",sep=""),tolower(colnames(my.orig.df))))>0){
      g=grep(paste("^",tolower(colnames(my.new.df)[i]),"$",sep=""),tolower(colnames(my.orig.df)))
      attributes(my.new.df[[i]])=attributes(my.orig.df[[g]])
    }else{
      print(paste(colnames(my.new.df)[i],"does not have a match."))
    }
  }
  return(my.new.df)
}

#' Remove rows with all NA
#'
#' @param my.df Dataframe
#'
#' @return Dataframe
#' @export
#'

quick.na.row=function(my.df){
  #### Remove rows with all NA
  my.df=my.df[rowSums(is.na(my.df))!=ncol(my.df),]
  return(my.df)
}

#' Change certain numbers to NA
#'
#' When NA do not come through from the SPSS dataset, apply value to all.
#'
#' @param my.df Dataframe
#' @param nums list of numeric numbers to turn to NA indiscriminantly
#'
#' @return Dataframe
#' @export
#'

quick.to.na=function(my.df,nums){

  for(i in 1:length(nums)){
    #### Remove 9's
    my.df=sapply(my.df, function(x) ifelse(x==nums[i],NA,x))
    my.df=as.data.frame(my.df)
  }
  return(my.df)
}


#' Quick VIM Aggregate Missing Plot
#'
#' Simple code, but used often enough that needed to be put in a function.
#' See the missing by calling the object. See graph with plot().
#'
#' @param my.df Dataframe
#'
#' @return VIM Object
#' @export
#'

quick.VIM=function(my.df){
#### VIM
library(VIM)
aggr_plot <- aggr(my.df, col=c('navyblue','red'), numbers=TRUE,
                  sortVars=TRUE, labels=names(my.df), cex.axis=.7, gap=3,
                  ylab=c("Histogram of missing data","Pattern"))
}


#' Quick ROC curve and Area Under ROC Curve
#'
#' Takes GLM and computes ROC curve using ROCR library.
#'
#' @param my.glm Object of class "glm"
#' @param title Main title for graph (default: "ROC Curve")
#' @param cutoff Cutoff point? (default: 0.5)
#' @param do.plot Plot ROC Curve? (default: T)
#'
#' @return List with area under curve and plot
#' @export
#'

quick.ROC=function(my.glm,title="ROC Curve",cutoff=0.5,do.plot=T){
  # ROC Curve
  library(ROCR)
  p <- predict(my.glm,newdata=my.glm$data, type="response")
  p=ifelse(p>cutoff,1,0)
  pr <- prediction(p, as.factor(my.glm$y))

  if(do.plot){
    prf <- performance(pr, measure = "tpr", x.measure = "fpr")
    plot(prf,main=title)
  }

  auc <- performance(pr, measure = "auc")
  sens=performance(pr, measure="sens",x.measure="spec")
  my.frame=as.data.frame(matrix(ncol=4))
  my.frame=c(cutoff,sens@y.values[[1]][2],sens@x.values[[1]][2],auc@y.values[[1]])
  names(my.frame)=c("Cutoff","Sensitivity","Specificity","AUC")
  return(my.frame)
}

#' ROC Diagnositics table
#'
#' Looks at Sensitivity, Specificity, and Area Under the Curve
#' for a range of cutoff values.
#'
#' @param my.glm Object of class "glm"
#' @param vals List of cutoff values to check.
#'
#' @return Matrix of values for each cutoff value.
#' @export
#'
#' @examples
#' quick.ROC.diag(my.glm,seq(0.4,0.6,by=0.05))

quick.ROC.diag=function(my.glm,vals){
  return(t(sapply(vals,quick.tasks::quick.ROC,my.glm=my.glm,do.plot=F,title="ROC Curve")))
}


#' Quick SAS Factor Labels
#'
#' Why should you have to type them? A dialog will come up to ask you for the file.
#'
#' @param my.df Data frame to get the new information
#'
#' @return Dataframe with label information
#' @export
#'

quick.SAS.labels=function(my.df){
  stupid.sas=file.choose()
  stupid.sas.commands=readLines(stupid.sas)
  sas.line=1
  while(sas.line<length(stupid.sas.commands)){
    if({length(grep(";",stupid.sas.commands[sas.line]))>0} | {length(grepl("^\\s*$",stupid.sas.commands[sas.line]))>0}){
      #Do nothing
      sas.line=sas.line+1
    }
    if(length(grep("value ",stupid.sas.commands[sas.line]))>0){
      my.temp.line=strsplit(stupid.sas.commands[sas.line]," ")
      my.var=my.temp.line[[1]][length(my.temp.line[[1]])]
      my.grep=grep(tolower(paste("^",substr(my.var,1,{nchar(my.var)-1}),"$",sep="")),tolower(names(my.df)))

      #Get labels
      sas.line=sas.line+1
      my.factor.num=NULL
      my.factor.label=NULL
      while(length(grep(";",stupid.sas.commands[sas.line]))!=1){
        my.temp.line=strsplit(stupid.sas.commands[sas.line],"\\=")
        my.factor.num=c(my.factor.num,as.numeric(my.temp.line[[1]][1]))
        my.factor.label=c(my.factor.label,as.character(substr(my.temp.line[[1]][2],ifelse(is.numeric(my.temp.line[[1]][2][5]),ifelse(is.numeric(my.temp.line[[1]][2][6]),7,6),5),{nchar(my.temp.line[[1]][2])-1})))
        sas.line=sas.line+1
      }
      names(my.factor.num)=my.factor.label
      if(length(my.grep)>0){
        attr(my.df[[my.grep]],"labels")=my.factor.num
      }
    }
  }
}
