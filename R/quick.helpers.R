####################### QUICK.HELP #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################

#' Table Printer
#'
#' Helper function to turn tables into html tables
#'
#' @param my.table Final table ready to be turned into html
#' @param type Type (lm, glm, manova, ordinal)
#' @param test.stat Test stat used in MANOVA
#' @param print.type "full" for full html, "part" for just <div>
#' @param the.caption Caption for top of table
#' @param the.footer
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.table=function(my.table,
                     type,
                     test.stat="Pillai",
                     print.type="full",
                     the.caption=NULL,
                     the.footer=NULL,
                     abbrev.length=ab.len,
                     SS.type=2,
                     new.rownames.int=NULL,
                     new.rownames.treat=NULL,
                     swap.na=NULL,
                     round.num=2,
                     col.names=my.colnames,
                     print.now=T){

  #### Inits ####

  if(type=="ord"){
    ab.len=30
    library(ordinal)
  }else{
    ab.len=15
  }

  attr(my.table,"quick.print.type")=print.type
  attr(my.table,"quick.abbrev.length")=abbrev.length
  attr(my.table,"class")=c(attr(my.table,"class"),"quick.table")


  if(type=="manova" | type=="stats::manova"){
    round.rows=c(2,3,4)
    p.val.row=8
    my.colnames=c("Variable",paste(test.stat,"<br />Statistic"),"F-Value",
                  paste("Type ",ifelse(SS.type==2,"II","III"),"<br />Sums of<br />Squares"),
                  "dF","Mult<br />dF","Resid<br />dF","Pr(>F)")
  }

  #### Find Intercept, Treatment, Total Change Locations ####
  int.loc=grep("Intercept",my.table[[1]])
  treat.loc=grep("Treatment Change",my.table[[1]])
  total.loc=grep("Total Change",my.table[[1]])
  end.loc=grep("^Total$",my.table[[1]])

  #### Swap value for NA (i.e. remove it) ####
  if(!is.null(swap.na)){
    na.vals=strsplit(swap.na)[[1]]
  }

  #### Replace rownames ####
  if(!is.null(new.rownames.int)){
    #### Do later
  }

  #### Round ####
  for(i in 1:length(round.rows)){
    my.table[[round.rows[i]]]=round(as.numeric(my.table[[round.rows[i]]]),digits =round.num)
  }

  #### P-val ####
  for(i in 1:end.loc){
    if(!is.na(my.table[i,p.val.row])){
      my.table[i,p.val.row]=round(as.numeric(my.table[i,p.val.row]),digits=3)
      if(my.table[i,p.val.row]==0){my.table[i,p.val.row]="<.001"}
      if(my.table[i,p.val.row]==1){my.table[i,p.val.row]=">.999"}
    }
  }

  #### Change NA to &nbsp; ####
  my.table=replace(my.table,is.na(my.table),"&nbsp;")

  #### Add style and basic tag structure####
  attr(my.table,"quick.doctype")=paste("<!DOCTYPE html --- Created with quick.tasks by Christopher Kraner>")
  attr(my.table,"quick.full.start")=paste("<html><head><style>table{border: 1px solid black;border-collapse: collapse;}",
                                          "th{padding: 15px;}td {padding: 5px;}#red {border: 2px solid red;}",
                                          "#black {border: 2px solid black;}#change {border-top: 1px solid black;text-align: left;}",
                                          "tr:hover {background-color: #f5f5f5;}#int {border-top: 1px solid black;}#col {border-bottom: 1px solid black;}</style></head>",sep="")
  attr(my.table,"quick.part.start")=paste("<div style=\"overflow-x:auto;\"><table style=\"width:100%\">")
  attr(my.table,"quick.caption")=ifelse(is.null(the.caption),NA,paste("<caption>",the.caption,"</caption>"))
  attr(my.table,"quick.part.end")=paste("</table></div>")
  attr(my.table,"quick.full.end")=paste("</body></html>")

  ##### Start table ####
  if(print.type=="full"){
    my.html.table=paste(attr(my.table,"quick.doctype"),attr(my.table,"quick.full.start"),attr(my.table,"quick.part.start"))
  }else{
    my.html.table=paste(attr(my.table,"quick.part.start"),attr(my.table,"quick.doctype"))
  }

  #### Put in caption
  if(!is.null(the.caption)){
    my.html.table=paste(my.html.table,attr(my.table,"quick.caption"))
  }

  #### Put in Column Headings
  if(type!="glm"){
    col.headings="<tr id=\"col\", align=\"center\">"
    for(i in 1:length(my.colnames)){
      col.headings=paste(col.headings,"<th>",my.colnames[i],"</th>")
    }
    col.headings=paste(col.headings,"</tr>")
  }else{
    stop("Sorry, not to GLM yet.")
  }
  my.html.table=paste(my.html.table,col.headings)

  #### Put in Rows
  for(i in 1:end.loc){

    #### Variable name
    if(i==1 & my.table[1,1]=="Intercept Change"){
      #### GLM stuff
      #### Make th, add id="change"
      my.line=paste("<tr id=\"int\"><th>",my.table[i,1],"</th>")
    }else if(i==treat.loc | i==total.loc){
      my.line=paste("<tr id=\"change\"><td align=\"left\"><b>",my.table[i,1],"</b></td>")
    }else if(i==1 & i==int.loc){
      my.line=paste("<tr id=\"int\"><td>",my.table[i,1],"</td>")
    }else if(i>total.loc){
      my.line=paste("<tr><td align=\"left\"><b>",my.table[i,1],"</b></td>")
    }else{
      my.line=paste("<tr><td>",my.table[i,1],"</td>")
    }

    #### Rest of row
    for(j in 2:p.val.row){
      my.line=paste(my.line,"<td>",my.table[i,j],"</td>")
    }

    my.line=paste(my.line,"</tr>")

    my.html.table=paste(my.html.table,my.line)
    if(i==1){
      attr(my.table,"quick.rows")=my.line
    }else{
      attr(my.table,"quick.rows")=paste(attr(my.table,"quick.rows"),my.line)
    }
  }

  #### Put in custom bottom


  #### Put in end
  my.html.table=paste(my.html.table,"</table></div>")
  if(print.type=="full"){
    my.html.table=paste(my.html.table,"</body></html>")
  }

  if(print.now){
    tempDir <- tempfile()
    dir.create(tempDir)
    htmlFile <- file.path(tempDir, "index.html")
    writeLines(my.html.table,htmlFile)
    viewer <- getOption("viewer")
    viewer(htmlFile)
  }

  return(my.table)
}



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


#' Create latent variables
#'
#' Helper function to make list of null model, pre models,
#' and full models for calculating type II SS (and eventually
#' type III SS).
#'
#' @param SSCP.list Model to be used
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.latent=function(SSCP.list){
  library(purrr)
  temp.latent.SSCP.eigen=purrr::map(SSCP.list[[1]],quick.latent.reduce)
  temp.latent.SSCP.error=quick.latent.reduce(SSCP.list[[2]])
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
quick.SSCP=function(my.summary){
  library(car)
  car.sum=car::Anova(my.summary,type=3)
  car.sum.SSCP=list(car.sum$SSP,car.sum$SSPE)
  return(car.sum.SSCP)
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
