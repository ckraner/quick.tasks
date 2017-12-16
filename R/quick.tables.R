###################### QUICK.TABLES ########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################


#' Contrast Tables in Pixiedust
#'
#' Beautiful tables using PHIA and PixieDust
#' for lm, glm, and mancova.
#'
#' @param my.model Model.
#' @param SS.type Type of sums of squares to report, either 2 or 3. (default = 3)
#' @param adjustment P-value adjustment sent to p.adjust (default = "bonferroni")
#' @param test.stat Character string giving the type of test statistic to be used (in MANOVA). (default="Wilks")
#' @param abbrev.length Integer telling how long of a label is too long. Longer than this and will be condensed (default=15)
#' @param pix.int Should this be viewable or is this for a document/dashboard? (default=T)
#' @param pix.method Print method. (default="html")
#' @param manova Is this a MAN(C)OVA?
#' @param my.factors If you only want some of the factors, use this. Otherwise, factors are pulled from the regression model
#' @param my.type If you have problems, you can specify the regression type. This is pulled from the model
#' @return Either pixiedust object or code (in HTML or LaTeX) for table
#' @keywords Explore
#' @export
#' @examples
#' quick.contrast(my.model)

quick.contrast=function(my.model, SS.type=3, adjustment="bonferroni",test.stat="Wilks", abbrev.length=15, pix.int=T,pix.method="html",
                        my.factors=my.contrasts,my.type=my.reg.type){

  #### Find type
  my.call=as.character(my.model$call)
  my.split.call=strsplit(my.call,"\\(")
  my.reg.type=my.split.call[[1]][1]


  #### Find factors
  my.contrasts=names(my.model$contrasts)

  if(my.type!="manova" | my.type!="stats::manova"){
    #### Find levels
    my.x.levels=NULL
    for(i in 1:length(my.factors)){
      my.x.levels=c(my.x.levels,my.model$xlevels[[i]])
    }

    #### Find num non var
    my.non.var=length(my.model$coefficients)-length(my.x.levels)-1+length(my.contrasts)
  }

  library(pixiedust)
  library(phia)
  if(my.type=="manova" | my.type=="stats::manova"){
    my.phia.print=as.data.frame(matrix(ncol=8,nrow=1))
    my.lengths=NULL
    this.table.var=1

    while(this.table.var<{length(my.factors)+1}){

      my.phia=phia::testInteractions(my.model,fixed=my.factors[this.table.var],adjustment = adjustment, abbrev.levels = abbrev.length)
      my.phia$names=attr(my.phia,"row.names")
      my.phia=my.phia[c("names","Df","test stat","approx F","num Df","den Df","Pr(>F)")]
      attr(my.phia,"class")=attr(my.phia,"class")[-1]
      my.lengths=c(my.lengths,nrow(my.phia))
      this.stuff=c(my.factors[this.table.var],NA)

      if(my.lengths[this.table.var]>2){

        for(i in 1:{my.lengths[this.table.var]-2}){

          this.stuff=c(this.stuff,NA)

        }

      }

      my.phia=cbind(this.stuff,my.phia)

      if(this.table.var==1){

        my.phia.print=my.phia

      }else{

        my.phia.print=rbind(my.phia.print,my.phia)

      }

      this.table.var=this.table.var+1
      #print(my.phia.print)

    }

    rownames(my.phia.print)=NULL
    phia.length=dim(my.phia.print)[1]

    options(pixie_interactive = pix.int)

    my.phia.pixie=pixiedust::dust(my.phia.print)%>%
      sprinkle_print_method(pix.method)%>%
      sprinkle(cols="Pr(>F)",fn=quote(pvalString(value,digits=3,format="default")))%>%
      sprinkle(cols="test stat",round=3)%>%
      sprinkle(cols="approx F",round=3)%>%
      sprinkle_colnames("","Levels","df","Pillai <br /> Statistic","approx <br /> F-value","num <br /> df","den <br /> df","Pr(>F)")%>%
      sprinkle(cols=1:8,rows={sum(my.lengths)},border="bottom")%>%
      sprinkle(cols=1:8,pad=10)%>%
      sprinkle(cols=1,rows=1:{sum(my.lengths)},border="left")%>%
      sprinkle(cols=3:8,rows=1:{sum(my.lengths)},border=c("right","left"))%>%
      sprinkle(cols=1:8,rows=1,border=c("top","bottom"),part="head")%>%
      sprinkle(cols=1,rows=1,border="left",part="head")%>%
      sprinkle(cols=3:8,rows=1,border=c("right","left"),part="head")%>%
      sprinkle_na_string(na_string="")%>%
      sprinkle_width(cols=1,rows=1:2,width=70,width_units="pt")%>%
      sprinkle_width(cols=2,rows=1:2,width=70,width_units="pt")%>%
      sprinkle_width(cols=3,width=30,width_units="pt")%>%
      sprinkle_width(cols=4,width=60,width_units="pt")%>%
      sprinkle_width(cols=5,width=50,width_units="pt")%>%
      sprinkle_width(cols=6,width=50,width_units="pt")%>%
      sprinkle_width(cols=7,width=50,width_units="pt")%>%
      sprinkle_width(cols=8,width=70,width_units="pt")%>%
      sprinkle(rows=1,halign="center",part="head")

    adj.method=as.data.frame(matrix(ncol=8,nrow=1))
    adj.method[1,]=c(paste("p adjustment method: ",adjustment,sep=""),NA,NA,NA,NA,NA,NA,NA)
    my.phia.pixie=redust(my.phia.pixie,adj.method,part="foot")%>%
      sprinkle(merge=T,halign="center",part="foot")

    if(pix.int){
      return(my.phia.pixie)
    }else{
      my.phia.pixie=print(my.phia.pixie,quote=F)[1]
      return(my.phia.pixie)
    }
  }else if(my.type=="lm" | my.type=="stats::lm"){
    my.phia.print=as.data.frame(matrix(ncol=7,nrow=1))
    my.lengths=NULL
    this.table.var=1

    while(this.table.var<{length(my.factors)+1}){

      my.phia=phia::testInteractions(my.model,pairwise=my.factors[this.table.var],adjustment = adjustment, abbrev.levels = abbrev.length)
      my.phia=my.phia[-{dim(my.phia)[1]},]
      my.phia$names=attr(my.phia,"row.names")
      my.phia=my.phia[c("names","Value","Df","Sum of Sq","F","Pr(>F)")]
      attr(my.phia,"class")=attr(my.phia,"class")[-1]
      my.lengths=c(my.lengths,{nrow(my.phia)})
      this.stuff=c(my.factors[this.table.var],NA)

      if(my.lengths[this.table.var]>2){

        for(i in 1:{my.lengths[this.table.var]-2}){

          this.stuff=c(this.stuff,NA)

        }

      }

      my.phia=cbind(this.stuff,my.phia)

      if(my.lengths[this.table.var]==1){
        my.phia=my.phia[-2,]
      }
      if(this.table.var==1){

        my.phia.print=my.phia

      }else{

        my.phia.print=rbind(my.phia.print,my.phia)

      }

      this.table.var=this.table.var+1
      #print(my.phia.print)

    }
    #print(my.phia.print)
    #my.phia.print=my.phia.print[-{dim(my.phia.print)[1]-1},]
    rownames(my.phia.print)=NULL
    phia.length=dim(my.phia.print)[1]

    options(pixie_interactive = pix.int)

    my.phia.pixie=pixiedust::dust(my.phia.print)%>%
      sprinkle_print_method(pix.method)%>%
      sprinkle(cols="Pr(>F)",fn=quote(pvalString(value,digits=3,format="default")))%>%
      sprinkle(cols="test stat",round=3)%>%
      sprinkle(cols="approx F",round=3)%>%
      sprinkle_colnames("","Levels","Value","df","Sums of <br /> Squares","F-value","Pr(>F)")%>%
      sprinkle(cols=1:7,rows={sum(my.lengths)},border="bottom")%>%
      sprinkle(cols=1:7,pad=10)%>%
      sprinkle(cols=1,rows=1:{sum(my.lengths)},border="left")%>%
      sprinkle(cols=3:7,rows=1:{sum(my.lengths)},border=c("right","left"))%>%
      sprinkle(cols=1:7,rows=1,border=c("top","bottom"),part="head")%>%
      sprinkle(cols=1,rows=1,border="left",part="head")%>%
      sprinkle(cols=3:7,rows=1,border=c("right","left"),part="head")%>%
      sprinkle_na_string(na_string="")%>%
      sprinkle_width(cols=1,rows=1,width=70,width_units="pt")%>%
      sprinkle_width(cols=2,rows=1,width=70,width_units="pt")%>%
      sprinkle_width(cols=3,width=50,width_units="pt")%>%
      sprinkle_width(cols=4,width=30,width_units="pt")%>%
      sprinkle_width(cols=5,width=60,width_units="pt")%>%
      sprinkle_width(cols=6,width=50,width_units="pt")%>%
      sprinkle_width(cols=7,width=50,width_units="pt")%>%
      sprinkle(rows=1,halign="center",part="head")

    adj.method=as.data.frame(matrix(ncol=7,nrow=1))
    adj.method[1,]=c(paste("p adjustment method: ",adjustment,sep=""),NA,NA,NA,NA,NA,NA)
    my.phia.pixie=redust(my.phia.pixie,adj.method,part="foot")%>%
      sprinkle(merge=T,halign="center",part="foot")

    if(pix.int){
      return(my.phia.pixie)
    }else{
      my.phia.pixie=print(my.phia.pixie,quote=F)[1]
      return(my.phia.pixie)
    }
  }else{
    my.phia.print=as.data.frame(matrix(ncol=6,nrow=1))
    my.lengths=NULL
    this.table.var=1

    while(this.table.var<{length(my.factors)+1}){

      my.phia=phia::testInteractions(my.model,pairwise=my.factors[this.table.var],adjustment = adjustment, abbrev.levels = abbrev.length)
      my.phia=my.phia[-{dim(my.phia)[1]},]
      my.phia$names=attr(my.phia,"row.names")
      my.phia=my.phia[c("names","Value","Df","Chisq","Pr(>Chisq)")]
      attr(my.phia,"class")=attr(my.phia,"class")[-1]
      my.lengths=c(my.lengths,{nrow(my.phia)})
      this.stuff=c(my.factors[this.table.var],NA)

      if(my.lengths[this.table.var]>2){

        for(i in 1:{my.lengths[this.table.var]-2}){

          this.stuff=c(this.stuff,NA)

        }

      }

      my.phia=cbind(this.stuff,my.phia)

      if(my.lengths[this.table.var]==1){
        my.phia=my.phia[-2,]
      }
      if(this.table.var==1){

        my.phia.print=my.phia

      }else{

        my.phia.print=rbind(my.phia.print,my.phia)

      }

      this.table.var=this.table.var+1
      #print(my.phia.print)

    }
    #print(my.phia.print)
    #my.phia.print=my.phia.print[-{dim(my.phia.print)[1]-1},]
    rownames(my.phia.print)=NULL
    phia.length=dim(my.phia.print)[1]
    my.phia.print[[6]]=as.numeric(my.phia.print[[6]])


    options(pixie_interactive = pix.int)
    my.phia.pixie=pixiedust::dust(my.phia.print)%>%
      sprinkle_print_method(pix.method)%>%
      sprinkle(cols="Pr(>Chisq)",fn=quote(pvalString(value,digits=3,format="default")))%>%
      sprinkle(cols="value",round=3)%>%
      sprinkle(cols="Chisq",round=3)%>%
      sprinkle_colnames("","Levels","Value","df","Chi-Sq","Pr(>Chisq)")%>%
      sprinkle(cols=1:6,rows={sum(my.lengths)},border="bottom")%>%
      sprinkle(cols=1:6,pad=10)%>%
      sprinkle(cols=1,rows=1:{sum(my.lengths)},border="left")%>%
      sprinkle(cols=3:6,rows=1:{sum(my.lengths)},border=c("right","left"))%>%
      sprinkle(cols=1:6,rows=1,border=c("top","bottom"),part="head")%>%
      sprinkle(cols=1,rows=1,border="left",part="head")%>%
      sprinkle(cols=3:6,rows=1,border=c("right","left"),part="head")%>%
      sprinkle_na_string(na_string="")%>%
      sprinkle_width(cols=1,rows=1,width=70,width_units="pt")%>%
      sprinkle_width(cols=2,rows=1,width=70,width_units="pt")%>%
      sprinkle_width(cols=3,width=50,width_units="pt")%>%
      sprinkle_width(cols=4,width=30,width_units="pt")%>%
      sprinkle_width(cols=5,width=60,width_units="pt")%>%
      sprinkle_width(cols=6,width=50,width_units="pt")%>%
      sprinkle(rows=1,halign="center",part="head")

    adj.method=as.data.frame(matrix(ncol=6,nrow=1))
    adj.method[1,]=c(paste("p adjustment method: ",adjustment,sep=""),NA,NA,NA,NA,NA)
    my.phia.pixie=redust(my.phia.pixie,adj.method,part="foot")%>%
      sprinkle(merge=T,halign="center",part="foot")
    if(pix.int){
      return(my.phia.pixie)
    }else{
      my.phia.pixie=print(my.phia.pixie,quote=F)[1]
      return(my.phia.pixie)
    }
  }
}


#' Regression Tables in Pixiedust
#'
#' Beautiful tables adding sums of squares
#' and p-value formatting, then giving html or
#' latex output. If interactive and HTML, will show up
#' in viewer.
#'
#' @param my.model Model. NOTE: If have factors, please place them first in the regression.
#' @param myDF Dataframe
#' @param my.factor If there are any factors, list them here.
#' @param SS.type Type of sums of squares to report (default = 3)
#' @param pix.int Should this be viewable or is this for a document/dashboard? (default=T)
#' @param pix.method Print method. (default="html")
#' @param type Type of regression? Currently supported: lm, glm (binary), manova
#' @return Either pixiedust object or code (in HTML or LaTeX) for table
#' @keywords Explore
#' @export
#' @examples
#' quick.reg(my.model, myDF)

quick.reg=function(my.model, myDF=my.found.df, SS.type=3, abbrev.length=15, pix.int=T,pix.method="html",
                   type=my.reg.type,test.stat="Wilks",my.factor=NULL,part.eta=F,VIF=F){
  library(pixiedust)
  library(broom)
  library(car)

  #### Find type
  my.call=as.character(my.model$call)
  my.split.call=strsplit(my.call,"\\(")
  my.reg.type2=my.split.call[[1]][1]
  if(my.reg.type2=="lm" | my.reg.type2 == "stats::lm"){
    my.reg.type="lm"
  }else if(my.reg.type2=="glm" | my.reg.type2== "stats::glm"){
    my.reg.type="glm"
  }else if(my.reg.type2=="manova" | my.reg.type2=="stats::manova"){
    my.reg.type="manova"
  }else if(my.reg.type2=="clm" | my.reg.type2=="ordinal::clm"){
    my.reg.type="ord"
    library(ordinal)
  }else{
    stop("Type not supported")
  }
  my.found.df=my.model$call$data
  if(is.null(my.found.df)){
    stop(paste("No data frame found"))
  }

  #### Make factor list
  if(type=="manova" | type=="stats::manova"){
    x3=capture.output(car::Anova(my.model,type=SS.type,test=test.stat))
    my.manova.test=data.frame(matrix(ncol=7,nrow=1))
    my.var.temp=4

    while(my.var.temp<{length(x3)-1}){

      test=strsplit(x3[my.var.temp],"\\s+")

      if(length(test[[1]])==9){

        test2=test[[1]][-9]
        test2=test2[-7]

      }else if(length(test[[1]])==8){

        test2=test[[1]][-8]

      }else{

        test2=test[[1]]
      }

      my.manova.test[{my.var.temp-3},]=test2
      my.var.temp=my.var.temp+1

    }

    my.manova.test[[2]]=as.numeric(my.manova.test[[2]])
    my.manova.test[[3]]=as.numeric(my.manova.test[[3]])
    my.manova.test[[4]]=as.numeric(my.manova.test[[4]])
    my.manova.test[[5]]=as.numeric(my.manova.test[[5]])
    my.manova.test[[6]]=as.numeric(my.manova.test[[6]])
    my.manova.test[[7]]=as.numeric(my.manova.test[[7]])

    options(pixie_interactive = pix.int)
    my.manova.pixie=pixiedust::dust(my.manova.test)%>%
      sprinkle_print_method(pix.method)%>%
      sprinkle(cols="X7",fn=quote(pvalString(value,digits=3,format="default")))%>%
      sprinkle(cols="X3",round=3)%>%
      sprinkle_colnames("","df",paste(test.stat," <br /> Statistic"),"approx <br /> F-value","num df","den df","Pr(>F)")%>%
      sprinkle(cols=1:7,rows={length(x3)-5},border=c("bottom","left","right"))%>%
      sprinkle(cols=1:7,pad=10)%>%
      sprinkle(cols=1:7,rows=1:{length(x3)-5},border=c("left","right"))%>%
      sprinkle(cols=1:7,rows=1,border=c("top","bottom","left","right"),part="head")%>%
      sprinkle(cols=1,rows=1,border="left",part="head")%>%
      sprinkle(cols=7,rows=1,border="right",part="head")%>%
      sprinkle_width(cols=1,rows=1:2,width=90,width_units="pt")%>%
      sprinkle_width(cols=2,rows=1:2,width=30,width_units="pt")%>%
      sprinkle_width(cols=3,width=60,width_units="pt")%>%
      sprinkle_width(cols=4,width=60,width_units="pt")%>%
      sprinkle_width(cols=5,width=50,width_units="pt")%>%
      sprinkle_width(cols=6,width=50,width_units="pt")%>%
      sprinkle_width(cols=7,width=70,width_units="pt")%>%
      sprinkle(rows=1,halign="center",part="head")


    if(pix.int){
      return(my.manova.pixie)
    }else{
      my.manova.pixie=print(my.manova.pixie,quote=F)[1]
      return(my.manova.pixie)
    }
  }else{
    #### Use car::Anova to get SS Type 3

    my.summary=summary(my.model)
    my.coefficients=my.summary$coefficients
    my.coefficients=as.data.frame(my.coefficients)
    if(type=="ord" & length(my.model$model)==1){
      my.III.summary=NULL
      the.length=length(my.model$y.levels)-1
    }else{
    my.III.summary=car::Anova(my.model,type=SS.type)

    if(type=="glm" & is.null(my.factor)){
      the.length=dim(my.III.summary)[1]+1
    }else{
      the.length=dim(my.III.summary)[1]
    }
}

    if(type=="lm"){
      #### Calculate total SS
      my.total=sum(my.III.summary$`Sum Sq`[2:length(my.III.summary$`Sum Sq`)])
      my.df.total=sum(my.III.summary$Df[2:length(my.III.summary$Df)])
      total.intercepts=1
      my.rownames=c(abbreviate(rownames(my.summary$coefficients),minlength = abbrev.length),"Residuals","Total")
    }else if(type=="glm"){
      #### Calculate model deviance stats
      ddeviance2=my.model$null.deviance-my.model$deviance
      ddf2=my.model$df.null-my.model$df.residual
      fit2=1-pchisq(ddeviance2,ddf2)
      total.intercepts=1
      my.rownames=c(abbreviate(rownames(my.summary$coefficients),minlength = abbrev.length),"Total")
    }else if(type=="ord"){
      my.temp.ord=update(my.model,~1)
      ddeviance2=my.model$logLik-my.temp.ord$logLik
      ddf2=my.model$edf-my.temp.ord$edf
      fit2=1-pchisq(ddeviance2,ddf2)
      total.intercepts=my.temp.ord$edf
      my.rownames=c(abbreviate(rownames(my.summary$coefficients),minlength=abbrev.length),"Total")
    }else{
      print("Error")
      return()
    }


    if(!VIF & !part.eta){
      v.p.len=7
      v.p.rep=0
    }else if(!VIF){
      v.p.len=8
      v.p.rep=1
    }else if(!part.eta){
      v.p.len=8
      v.p.rep=1
    }else{
      v.p.len=9
      v.p.rep=2
    }
    my.std.error=c(my.coefficients$`Std. Error`,NA)
    my.estimate=c(my.coefficients$Estimate,NA)
    my.tables.df=as.data.frame(matrix(ncol=v.p.len,nrow=1))

    if(type=="lm"){
      if(!VIF & !part.eta){
      names(my.tables.df)=c("rownames","sumsq","df","est","std.err","f.val","p.val")
      }else if(!part.eta){
        my.VIF=car::vif(my.model)
        names(my.tables.df)=c("rownames","sumsq","df","est","std.err","f.val","p.val","VIF")
      }else if(!VIF){
        names(my.tables.df)=c("rownames","sumsq","df","est","std.err","f.val","p.val","p.eta")
      }else{
        my.VIF=car::vif(my.model)
        names(my.tables.df)=c("rownames","sumsq","df","est","std.err","f.val","p.val","p.eta","VIF")
      }
    }else{
      names(my.tables.df)=c("var","od.rat","std.err","z.val","dev","df","p.val")
    }

    #### Make the double table entries ####

    #### Was very annoying...took my frustration out on
    #### Variable names

    factor.stupid=NULL
    factor.rownames=NULL
    num.of.levels=NULL
    ordinal.temp=0
    if(!is.null(my.factor)){

      for(i in 1:length(my.factor)){
        factor.stupid=c(factor.stupid,grep(paste("^",my.factor[i],"$",sep=""),names(myDF)))

        if(type=="lm"){
          factor.rownames=c(factor.rownames,grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary)))
          num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
          # }else if(type=="glm"){
          #   factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+1})
          #   num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
        }else{
          factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+total.intercepts+sum(num.of.levels)-ordinal.temp})
          num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
          ordinal.temp=ordinal.temp+1
        }
      }

    }else{

      factor.rownames=0L

    }


    my.factor.var=1
    this.temp.var=1
    this.shift.temp=1
    yet.another.var=1
    my.shift=0
    other.temp=2
    ord.temp=1
    if(type=="lm"){
      dang.length=length(rownames(my.III.summary))
    }else if (type=="glm"){
      dang.length=total.intercepts+length(rownames(my.III.summary))+max({sum(num.of.levels)-length(my.factor)},0)+1
    }else{
      dang.length=total.intercepts+length(rownames(my.III.summary))+max({sum(num.of.levels)-length(my.factor)},0)+1
      #dang.length=7
    }



    while(this.shift.temp<dang.length){

      if(is.na(factor.rownames[my.factor.var])){

        my.factor.rownames=1

      }else{

        my.factor.rownames=factor.rownames[my.factor.var]

      }
      if(this.shift.temp==1){
        i=1
        while(i<=total.intercepts){
          if(type=="lm"){
            my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
            my.df=my.III.summary$Df[this.shift.temp]
            my.est=my.estimate[this.shift.temp]
            my.std.err=my.std.error[this.shift.temp]
            my.f.val=my.III.summary$`F value`[this.shift.temp]
            my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
            my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val,rep(NA,v.p.rep))
          }else{
            my.or=exp(my.estimate[this.shift.temp])
            my.est=my.estimate[this.shift.temp]
            my.z.val=my.summary$coefficients[this.shift.temp,3]
            my.std.err=my.std.error[this.shift.temp]
            #my.dev=my.III.summary$`LR Chisq`[this.shift.temp]
            #my.df=my.III.summary$Df[this.shift.temp]
            my.p.val=my.summary$coefficients[{this.shift.temp},4]
            my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
          }
          this.shift.temp=this.shift.temp+1
          this.temp.var=this.temp.var+1
          i=i+1
        }
      }else if(this.shift.temp==my.factor.rownames){
        if(type=="lm"){
          my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
          my.df=my.III.summary$Df[this.shift.temp]
          my.est=NA
          my.std.err=NA
          my.z.val=NA
          my.f.val=my.III.summary$`F value`[this.shift.temp]
          my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
          my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.sumsq,my.df,my.std.err,my.z.val,my.f.val,my.p.val,rep(NA,v.p.rep))
        }else if(type=="glm"){
          my.or=NA
          my.est=NA
          my.std.err=NA
          my.dev=my.III.summary$`LR Chisq`[ord.temp]
          my.df=my.III.summary$Df[ord.temp]
          my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
          my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.or,my.est,my.std.err,my.dev,my.df,my.p.val)
          ord.temp=ord.temp+1
        }else{
          my.or=NA
          my.est=NA
          my.std.err=NA
          my.dev=my.III.summary$Chisq[ord.temp]
          my.df=my.III.summary$Df[ord.temp]
          my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
          my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.or,my.est,my.std.err,my.dev,my.df,my.p.val)
          ord.temp=ord.temp+1

        }
        yet.another.var=yet.another.var+1
        this.temp.var=this.temp.var+1
        this.shift.temp=this.shift.temp+1
        other.other.temp=2


        if(length(grepl(":",my.summary$coefficients[other.temp,1]))>0){

          while(other.other.temp<{num.of.levels[my.factor.var]+1}){
            if(type=="lm"){
              my.sumsq=NA
              my.df=NA
              my.est=my.estimate[other.temp]
              my.std.err=my.std.error[other.temp]
              #### NEED TO FIX ####
              my.f.val={my.summary$coefficients[other.temp,3]^2}
              my.p.val=my.summary$coefficients[other.temp,4]
              my.tables.df[this.temp.var,]=c(my.rownames[other.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val,rep(NA,v.p.rep))
            }else if(type=="glm"){
              my.or=exp(my.estimate[other.temp])
              my.est=my.estimate[other.temp]
              my.std.err=my.std.error[other.temp]
              my.z.val=my.summary$coefficients[other.temp,3]
              my.dev=my.III.summary$`LR Chisq`[other.temp]
              my.df=my.III.summary$Df[other.temp]
              my.p.val=my.summary$coefficients[other.temp,4]
              my.tables.df[this.temp.var,]=c(my.rownames[other.temp],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
              this.shift.temp=this.shift.temp+1
            }else{
              my.or=exp(my.estimate[other.temp+total.intercepts-1])
              my.est=my.estimate[other.temp+total.intercepts-1]
              my.std.err=my.std.error[other.temp+total.intercepts-1]
              my.z.val=my.summary$coefficients[{other.temp+total.intercepts-1},3]
              #my.dev=my.III.summary$Chisq[other.temp]
              #my.df=my.III.summary$Df[other.temp]
              my.p.val=my.summary$coefficients[{other.temp+total.intercepts-1},4]
              my.tables.df[this.temp.var,]=c(my.rownames[other.temp+total.intercepts-1],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
              this.shift.temp=this.shift.temp+1
            }
            this.temp.var=this.temp.var+1
            other.temp=other.temp+1
            other.other.temp=other.other.temp+1
            the.length=the.length+1

          }

        }else{

        }

        if(my.factor.var==1){

          my.shift={my.shift+other.temp-2}

        }else{

          my.shift=my.shift+other.temp

        }

        my.factor.var=my.factor.var+1

      }else{
        if(type=="lm"){
          my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
          my.df=my.III.summary$Df[this.shift.temp]
          my.est=my.estimate[this.shift.temp]
          my.std.err=my.std.error[this.shift.temp]
          my.f.val=my.III.summary$`F value`[this.shift.temp]
          my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
          if(!VIF & !part.eta){
          my.tables.df[this.temp.var,]=c(rownames(my.III.summary)[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val)
          }else if (!part.eta){
            my.tables.df[this.temp.var,]=c(rownames(my.III.summary)[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val,my.VIF[this.shift.temp-1,1])
          }else if(!VIF){
            my.p.eta=my.sumsq/my.total
            my.tables.df[this.temp.var,]=c(rownames(my.III.summary)[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val,my.p.eta)
          }else{
            my.p.eta=my.sumsq/my.total
            my.tables.df[this.temp.var,]=c(rownames(my.III.summary)[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val,my.p.eta,my.VIF[this.shift.temp-1,1])
          }
        }else if(type=="glm"){
          if(!is.null(my.factor)){
            this.shift.temp=this.shift.temp-ordinal.temp
          }
          my.or=exp(my.estimate[this.shift.temp])
          my.est=my.estimate[this.shift.temp]
          my.std.err=my.std.error[this.shift.temp]
          my.z.val=my.summary$coefficients[other.temp,3]
          my.dev=my.III.summary$`LR Chisq`[ord.temp]
          my.df=my.III.summary$Df[ord.temp]
          my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
          my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,NA,my.dev,my.df,my.p.val)
          if(!is.null(my.factor)){
            this.shift.temp=this.shift.temp+ordinal.temp
          }
          ord.temp=ord.temp+1
        }else{
          if(!is.null(my.factor)){
            this.shift.temp=this.shift.temp-ordinal.temp
          }
          my.or=exp(my.estimate[this.shift.temp])
          my.est=my.estimate[this.shift.temp]
          my.std.err=my.std.error[this.shift.temp]
          my.z.val=my.summary$coefficients[this.shift.temp,3]
          my.dev=my.III.summary$Chisq[ord.temp]
          my.df=my.III.summary$Df[ord.temp]
          my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
          my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,NA,my.dev,my.df,my.p.val)
          if(!is.null(my.factor)){
            this.shift.temp=this.shift.temp+ordinal.temp
          }
          ord.temp=ord.temp+1
        }
        this.shift.temp=this.shift.temp+1
        this.temp.var=this.temp.var+1

      }
    }

    if(type=="lm"){

      my.tables.df[this.temp.var,]=c("Residuals",my.III.summary$`Sum Sq`[this.shift.temp],my.III.summary$Df[this.shift.temp],NA,NA,NA,NA,rep(NA,v.p.rep))
      my.tables.df[this.temp.var+1,]=c("Total",my.total,my.df.total,NA,NA,NA,NA,rep(NA,v.p.rep))
      if(!VIF & !part.eta){
      }else if(!VIF){
        my.tables.df$p.eta=as.numeric(my.tables.df$p.eta)
      }else if(!part.eta){
        my.tables.df$VIF=as.numeric(my.tables.df$VIF)
      }else{
        my.tables.df$p.eta=as.numeric(my.tables.df$p.eta)
        my.tables.df$VIF=as.numeric(my.tables.df$VIF)
      }
      my.tables.df$f.val=as.numeric(my.tables.df$f.val)
      my.tables.df$est=as.numeric(my.tables.df$est)
      my.tables.df$sumsq=as.numeric(my.tables.df$sumsq)
    }else{
      my.tables.df[this.temp.var,]=c("Change from Null",NA,NA,NA,ddeviance2,ddf2,fit2)
      my.tables.df$od.rat=as.numeric(my.tables.df$od.rat)
      my.tables.df$z.val=as.numeric(my.tables.df$z.val)
      my.tables.df$dev=as.numeric(my.tables.df$dev)
    }

    my.tables.df$df=as.numeric(my.tables.df$df)
    my.tables.df$std.err=as.numeric(my.tables.df$std.err)
    my.tables.df$p.val=as.numeric(my.tables.df$p.val)

    #### Make custom glance stats ####

    #### Can eventually make it options
    if(type=="lm"){
      glance_stats=broom::glance(my.model)
      glance_stats=tidyr::gather(glance_stats)


      glance_stats[3:{5+v.p.rep}]=NA
      glance_stats[6+v.p.rep]=c(glance_stats$key[7:8],NA,glance_stats$key[9:11],NA,NA,NA,NA,NA)
      glance_stats[7+v.p.rep]=c(glance_stats$value[7:8],NA,glance_stats$value[9:11],NA,NA,NA,NA,NA)
      glance_stats=glance_stats[-7:-11,]
      glance_stats=glance_stats[-3,]
      glance_stats=glance_stats[-5,]
    }else{
      # glance_stats2=as.data.frame(matrix(ncol=7,nrow=1))
      # glance_stats2[1,]=c(glance_stats$key[1],glance_stats$value[1],NA,NA,NA,glance_stats$key[6],glance_stats$value[6])
      # glance_stats2[2,]=c(glance_stats$key[2],glance_stats$value[2],NA,NA,NA,glance_stats$key[7],glance_stats$value[7])
      # glance_stats2[3,]=c(glance_stats$key[3],glance_stats$value[3],NA,NA,NA,glance_stats$key[4],glance_stats$value[4])
      # glance_stats=glance_stats2
    }
    #### For total
    the.length=the.length+1
    this.temp.var=this.temp.var+1
    #### Make table ####

    options(pixie_interactive = pix.int,pixie_na_string="")


    if(type=="lm"){

      my.dust=pixiedust::dust(my.tables.df)%>%
        sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
        sprinkle_print_method(pix.method)%>%
        sprinkle_na_string()%>%
        sprinkle(rows=1:the.length,cols=1:v.p.len,
                 border="right",border_color="black")%>%
        sprinkle(rows=1:the.length,cols=1,
                 border="left",border_color="black")%>%
        sprinkle(rows=1,cols=1:v.p.len,
                 border=c("top","bottom","left","right"),border_color="black",part="head")%>%
        sprinkle(rows=this.temp.var,cols=1:v.p.len,
                 border="bottom",border_color="black")

      if(!VIF & !part.eta){
        my.dust=my.dust%>%
          sprinkle(cols=c("sumsq","est","std.err","f.val"),round=2)%>%
          sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
          sprinkle_colnames("Variable",paste("Type ",SS.type,"<br /> Sums of Sq",sep=""),"df","Estimate","Std. Error","F-value","Pr(>F)")%>%
          sprinkle_align(rows=1,halign="center",part="head")%>%
          sprinkle_pad(rows=1,pad=5,part="head")

      }else if(!VIF){
        my.dust=my.dust%>%
          sprinkle(cols=c("sumsq","est","std.err","f.val","p.eta"),round=2)%>%
          sprinkle(cols=c(2,3,4,5,6,7,8),pad=5)%>%
          sprinkle_colnames("Variable",paste("Type ",SS.type,"<br /> Sums of Sq",sep=""),"df","Estimate","Std. Error","F-value","Pr(>F)","Part <br /> eta")%>%
          sprinkle_align(rows=1,halign="center",part="head")%>%
          sprinkle_pad(rows=1,pad=5,part="head")

      }else if(!part.eta){
        my.dust=my.dust%>%
          sprinkle(cols=c("sumsq","est","std.err","f.val","VIF"),round=2)%>%
          sprinkle(cols=c(2,3,4,5,6,7,8),pad=5)%>%
          sprinkle_colnames("Variable",paste("Type ",SS.type,"<br /> Sums of Sq",sep=""),"df","Estimate","Std. Error","F-value","Pr(>F)","VIF")%>%
          sprinkle_align(rows=1,halign="center",part="head")%>%
          sprinkle_pad(rows=1,pad=5,part="head")

      }else{
        my.dust=my.dust%>%
          sprinkle(cols=c("sumsq","est","std.err","f.val","p.eta","VIF"),round=2)%>%
          sprinkle(cols=c(2,3,4,5,6,7,8),pad=5)%>%
          sprinkle_colnames("Variable",paste("Type ",SS.type,"<br /> Sums of Sq",sep=""),"df","Estimate","Std. Error","F-value","Pr(>F)","Part <br /> eta","VIF")%>%
          sprinkle_align(rows=1,halign="center",part="head")%>%
          sprinkle_pad(rows=1,pad=5,part="head")

      }

    }else if(type=="glm"){
      my.dust=pixiedust::dust(my.tables.df)%>%
        sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
        sprinkle_print_method(pix.method)%>%
        sprinkle_na_string()%>%
        sprinkle(cols=2:5,round=2)%>%
        sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
        sprinkle(cols=2:7,pad=5,part="head")%>%
        sprinkle(rows={this.temp.var-2},cols=1:7,border="bottom",border_color="black")%>%
        sprinkle(rows=1:the.length,cols=1:7,
                 border="right",border_color="black")%>%
        sprinkle(rows=1:the.length,cols=1,
                 border="left",border_color="black")%>%
        sprinkle(rows=1,cols=1:7,
                 border=c("top","bottom","left","right"),border_color="black",part="head")%>%
        sprinkle(rows=this.temp.var-1,cols=1:7,
                 border="bottom",border_color="black")%>%
        sprinkle_colnames("Variable","Odds Ratio","Std. Error","z-Value","Deviance","df","p-Value")

    }else{
      my.dust=pixiedust::dust(my.tables.df)%>%
        sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
        sprinkle_print_method(pix.method)%>%
        sprinkle_na_string()%>%
        sprinkle(cols=2:5,round=2)%>%
        sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
        sprinkle(cols=2:7,pad=5,part="head")%>%
        sprinkle(rows=total.intercepts,cols=1:7,border="bottom",border_color="black")%>%
        sprinkle(rows={this.temp.var-2},cols=1:7,border="bottom",border_color="black")%>%
        sprinkle(rows=1:{the.length+total.intercepts},cols=1:7,
                 border="right",border_color="black")%>%
        sprinkle(rows=1:{the.length+total.intercepts},cols=1,
                 border="left",border_color="black")%>%
        sprinkle(rows=1,cols=1:7,
                 border=c("top","bottom","left","right"),border_color="black",part="head")%>%
        sprinkle(rows=this.temp.var-1,cols=1:7,
                 border="bottom",border_color="black")%>%
        sprinkle_colnames("Variable","Odds Ratio","Std. Error","z-Value","Deviance","df","p-Value")
    }

    if(type=="lm"){
      my.dust=pixiedust::redust(my.dust,glance_stats,part="foot")%>%
        sprinkle(cols=c(2,{7+v.p.rep}),round=3,part="foot")%>%
        sprinkle(cols=3:{5+v.p.rep},replace=c("","","","","","","","","","","","",rep("",4*v.p.rep)),part="foot")%>%
        sprinkle(cols=1, replace=c("R-Square","Adj R-Sq","F-Statistic","P-Value"),part="foot")%>%
        sprinkle(cols=2,rows=4,fn=quote(pvalString(value,digits=3,format="default")),part="foot")%>%
        sprinkle(cols=1:v.p.len,rows=1,halign="center",part="head")%>%
        sprinkle_width(cols=1,width=90,width_units="pt")%>%
        sprinkle_width(cols=2,width=108,width_units="pt")%>%
        sprinkle_width(cols=4,width=62,width_units="pt")%>%
        sprinkle_width(cols=5,width=68,width_units="pt")%>%
        sprinkle_width(cols=6,width=68,width_units="pt")%>%
        sprinkle_width(cols=7,width=71,width_units="pt")%>%
        sprinkle(cols=2,halign="left",part="foot")
    }else{
      # my.dust=pixiedust::redust(my.dust,glance_stats,part="foot")%>%
      #   sprinkle(cols=c(2,7),round=2,part="foot")%>%
      #   sprinkle(cols=3:5,replace=c("","","","","","","","",""),part="foot")
    }
    if(pix.int){
      return(my.dust)
    }else{
      my.dust.print=print(my.dust,quote=F)[1]
      return(my.dust.print)
    }

  }
}



#' PixieDust Tables for Lavaan
#'
#' Places html tables in viewer for different areas of Lavaan output
#' for SEM. Based on NIU class 11/2017
#'
#' @param myfit Fit from Lavaan package
#' @return NULL
#' @keywords Explore
#' @export
#' @examples
#' quick.lavaan(myfit)

quick.lavaan=function(myfit){
  require(lavaan)
  require(pixiedust)
  require(tibble)

  prev.width=getOption("width")
  options(width=80)

  mysummary=capture.output(lavaan::summary(myfit,
                                           standardized = TRUE,rsq=T))
  #### Fit Table ####
  my.fit.table=as.data.frame(matrix(nrow=7,ncol=4))
  my.fit.table[1,]=c("Number of Iterations",lavInspect(myfit,what="iterations"),NA,NA)
  my.fit.table[2,]=c("Total Observations",lavInspect(myfit,what="ntotal"),NA,NA)
  my.fit.table[3,]=c("Chi-Sq Test of Fit",round(fitMeasures(myfit)[3],2),fitMeasures(myfit)[4],pixiedust::pval_string(fitMeasures(myfit)[5]))
  my.fit.table[4,]=c("Comparitive Fit Index",round(fitMeasures(myfit)[9],3),NA,NA)
  my.fit.table[5,]=c("Tucker-Lewis Index",round(fitMeasures(myfit)[10],3),NA,NA)
  my.fit.table[6,]=c("RMSEA",round(fitMeasures(myfit)[23],3),NA,pixiedust::pval_string(fitMeasures(myfit)[26]))
  my.fit.table[7,]=c("SRMR",round(fitMeasures(myfit)[29],3),NA,NA)
  colnames(my.fit.table)=c("Name","Value","df","p-val")
  #my.fit.table

  #### R Sq ####
  my.r2=lavaan::lavInspect(myfit,what="r2")
  my.r2=as.data.frame(my.r2)
  my.r2=tibble::rownames_to_column(my.r2)
  my.r2$my.r2=round(my.r2$my.r2,3)
  #my.r2

  #### Other tables ####
  lavaan.latent=NULL
  lavaan.covar=NULL
  lavaan.vari=NULL
  lavaan.reg=NULL
  lavaan.latent.temp=1
  lavaan.temp=1
  while(lavaan.temp<{length(mysummary)+1}){
    ### Latent Variables matrix
    if(length(grep("Latent Variables:",mysummary[lavaan.temp]))>0){
      lavaan.temp=lavaan.temp+2
      while((length(grep("Covariances:",mysummary[lavaan.temp]))==0) &&
            (length(grep("Regressions:",mysummary[lavaan.temp]))==0) &&
            (length(grep("Variances:",mysummary[lavaan.temp]))==0) &&
            (length(grep("R-Square:",mysummary[lavaan.temp]))==0)){
        str.temp=strsplit(mysummary[lavaan.temp],"  ")
        if(length(str.temp[[1]])>0){
          lavaan.latent[[lavaan.latent.temp]]=list()
          for(i in 1:(length(str.temp[[1]]))){
            if(nchar(str.temp[[1]][i])>0){
              lavaan.latent[[lavaan.latent.temp]]=c(lavaan.latent[[lavaan.latent.temp]],str.temp[[1]][i])
            }
          }
          lavaan.latent.temp=lavaan.latent.temp+1
        }
        lavaan.temp=lavaan.temp+1
      }

      if(length(grep("Regressions:",mysummary[lavaan.temp]))!=0){
        #### Regressions Matrix
        lavaan.temp=lavaan.temp+2
        lavaan.latent.temp=1
        while((length(grep("Covariances:",mysummary[lavaan.temp]))==0) &&
              (length(grep("Variances:",mysummary[lavaan.temp]))==0) &&
              (length(grep("R-Square:",mysummary[lavaan.temp]))==0)){
          str.temp=strsplit(mysummary[lavaan.temp],"  ")
          if(length(str.temp[[1]])>0){
            lavaan.reg[[lavaan.latent.temp]]=list()
            for(i in 1:(length(str.temp[[1]]))){
              if(nchar(str.temp[[1]][i])>0){
                lavaan.reg[[lavaan.latent.temp]]=c(lavaan.reg[[lavaan.latent.temp]],str.temp[[1]][i])
              }
            }
            lavaan.latent.temp=lavaan.latent.temp+1
          }
          lavaan.temp=lavaan.temp+1
        }
      }

      #### Covariance matrix
      if((length(grep("Covariances:",mysummary[lavaan.temp]))!=0)){
        lavaan.temp=lavaan.temp+2
        lavaan.latent.temp=1
        while((length(grep("Variances:",mysummary[lavaan.temp]))==0) &&
              (length(grep("R-Square:",mysummary[lavaan.temp]))==0)){
          str.temp=strsplit(mysummary[lavaan.temp],"  ")
          if(length(str.temp[[1]])>0){
            lavaan.covar[[lavaan.latent.temp]]=list()
            for(i in 1:(length(str.temp[[1]]))){
              if(nchar(str.temp[[1]][i])>0){
                lavaan.covar[[lavaan.latent.temp]]=c(lavaan.covar[[lavaan.latent.temp]],str.temp[[1]][i])
              }
            }
            lavaan.latent.temp=lavaan.latent.temp+1
          }
          lavaan.temp=lavaan.temp+1
        }
      }

      #### Variance matrix
      if((length(grep("Variances:",mysummary[lavaan.temp]))!=0)){
        lavaan.temp=lavaan.temp+2
        lavaan.latent.temp=1
        while(length(grep("R-Square:",mysummary[lavaan.temp]))==0){
          str.temp=strsplit(mysummary[lavaan.temp],"  ")
          if(length(str.temp[[1]])>0){
            lavaan.vari[[lavaan.latent.temp]]=list()
            for(i in 1:(length(str.temp[[1]]))){
              if(nchar(str.temp[[1]][i])>0){
                lavaan.vari[[lavaan.latent.temp]]=c(lavaan.vari[[lavaan.latent.temp]],str.temp[[1]][i])
              }
            }
            lavaan.latent.temp=lavaan.latent.temp+1
          }
          lavaan.temp=lavaan.temp+1
        }
      }

    }
    lavaan.temp=lavaan.temp+1
  }


  #### lavaan.covar
  if(length(lavaan.covar)>0){
    covar.table=as.data.frame(matrix(nrow=length(lavaan.covar),ncol=7))
    table.count=1
    for(i in 1:length(lavaan.covar)){
      if(length(lavaan.covar[[i]])<7){
        covar.table[table.count,]=c(lavaan.covar[[i]][1],NA,NA,NA,NA,NA,NA)
      }
      else{
        covar.table[table.count,]=c(lavaan.covar[[i]][1],lavaan.covar[[i]][2],lavaan.covar[[i]][3],lavaan.covar[[i]][6],lavaan.covar[[i]][7],lavaan.covar[[i]][4],lavaan.covar[[i]][5])
      }
      table.count=table.count+1
    }
    colnames(covar.table)=c("Name","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")
    #covar.table
  }else{
    covar.table=NULL
  }

  #### lavaan.latent
  if(length(lavaan.latent)>0){
    latent.table=as.data.frame(matrix(nrow=length(lavaan.latent),ncol=7))
    latent.table.count=1
    for(i in 1:length(lavaan.latent)){
      if(length(lavaan.latent[[i]])<3){
        latent.table[latent.table.count,]=c(lavaan.latent[[i]][[1]],NA,NA,NA,NA,NA,NA)
      }else if(length(lavaan.latent[[i]])==4){
        latent.table[latent.table.count,]=c(lavaan.latent[[i]][[1]],lavaan.latent[[i]][[2]],NA,lavaan.latent[[i]][[3]],lavaan.latent[[i]][[4]],NA,NA)
      }else{
        latent.table[latent.table.count,]=c(lavaan.latent[[i]][[1]],lavaan.latent[[i]][[2]],lavaan.latent[[i]][[3]],lavaan.latent[[i]][[6]],lavaan.latent[[i]][[7]],lavaan.latent[[i]][[4]],lavaan.latent[[i]][[5]])
      }
      latent.table.count=latent.table.count+1
    }
    colnames(latent.table)=c("Name","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")
    #latent.table
  }else{
    latent.table=NULL
  }

  #### lavaan.vari
  if(length(lavaan.vari)>0){
    vari.table=as.data.frame(matrix(nrow=length(lavaan.vari),ncol=7))
    vari.table.count=1
    for(i in 1:length(lavaan.vari)){
      vari.table[vari.table.count,]=c(lavaan.vari[[i]][[1]],lavaan.vari[[i]][[2]],lavaan.vari[[i]][[3]],lavaan.vari[[i]][[6]],lavaan.vari[[i]][[7]],lavaan.vari[[i]][[4]],lavaan.vari[[i]][[5]])
      vari.table.count=vari.table.count+1
    }
    colnames(vari.table)=c("Name","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")
  }else{
    vari.table=NULL
  }

  #### lavaan.reg
  if(length(lavaan.reg)>0){
    reg.table=as.data.frame(matrix(nrow=length(lavaan.reg),ncol=7))
    reg.table.count=1
    for(i in 1:length(lavaan.reg)){
      if(length(lavaan.reg[[i]])<7){
        reg.table[reg.table.count,]=c(lavaan.reg[[i]][1],NA,NA,NA,NA,NA,NA)
      }
      else{
        reg.table[reg.table.count,]=c(lavaan.reg[[i]][1],lavaan.reg[[i]][2],lavaan.reg[[i]][3],lavaan.reg[[i]][6],lavaan.reg[[i]][7],lavaan.reg[[i]][4],lavaan.reg[[i]][5])
      }
      reg.table.count=reg.table.count+1
    }
    colnames(reg.table)=c("Name","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")
  }else{
    reg.table=NULL
  }

  dusted.fit.table=pixiedust::dust(my.fit.table)%>%
    pixiedust::sprinkle_na_string(na_string = "")%>%
    pixiedust::sprinkle_print_method(print_method = "html")%>%
    pixiedust::sprinkle_colnames("","Value","df","P-Val")%>%
    pixiedust::sprinkle_border(border="all")%>%
    pixiedust::sprinkle_border(border="all",part="head")%>%
    pixiedust::sprinkle_pad(pad=7)%>%
    pixiedust::sprinkle_align(halign="center",part="head")
  #dusted.fit.table

  latent.table$`P(>|z|)`=as.numeric(latent.table$`P(>|z|)`)
  dusted.latent.table=pixiedust::dust(latent.table)%>%
    pixiedust::sprinkle_na_string(na_string = "")%>%
    pixiedust::sprinkle_print_method(print_method = "html")%>%
    pixiedust::sprinkle_colnames("Latent","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")%>%
    pixiedust::sprinkle_border(border="all")%>%
    pixiedust::sprinkle_border(border="all",part="head")%>%
    pixiedust::sprinkle_pad(pad=7)%>%
    pixiedust::sprinkle_pad(pad=7,part="head")%>%
    pixiedust::sprinkle(cols=7,fn = quote(pvalString(value)))%>%
    pixiedust::sprinkle_align(halign="center",part="head")
  #dusted.latent.table

  if(length(reg.table)>0){
    reg.table$`P(>|z|)`=as.numeric(reg.table$`P(>|z|)`)
    dusted.reg.table=pixiedust::dust(reg.table)%>%
      pixiedust::sprinkle_na_string(na_string = "")%>%
      pixiedust::sprinkle_print_method(print_method = "html")%>%
      pixiedust::sprinkle_colnames("Regression","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")%>%
      pixiedust::sprinkle_border(border="all")%>%
      pixiedust::sprinkle_border(border="all",part="head")%>%
      pixiedust::sprinkle_pad(pad=7)%>%
      pixiedust::sprinkle_pad(pad=7,part="head")%>%
      pixiedust::sprinkle(cols=7,fn = quote(pvalString(value)))%>%
      pixiedust::sprinkle_align(halign="center",part="head")
    #dusted.reg.table
  }else{
    dusted.reg.table=NULL
  }

  if(length(covar.table)>0){
    covar.table$`P(>|z|)`=as.numeric(covar.table$`P(>|z|)`)
    dusted.covar.table=pixiedust::dust(covar.table)%>%
      pixiedust::sprinkle_na_string(na_string = "")%>%
      pixiedust::sprinkle_print_method(print_method = "html")%>%
      pixiedust::sprinkle_colnames("Covariance","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")%>%
      pixiedust::sprinkle_border(border="all")%>%
      pixiedust::sprinkle_border(border="all",part="head")%>%
      pixiedust::sprinkle_pad(pad=7)%>%
      pixiedust::sprinkle_pad(pad=7,part="head")%>%
      pixiedust::sprinkle(cols=7,fn = quote(pvalString(value)))%>%
      pixiedust::sprinkle_align(halign="center",part="head")
    #dusted.covar.table
  }else{
    dusted.covar.table=NULL
  }

  if(length(vari.table)){
    vari.table$`P(>|z|)`=as.numeric(vari.table$`P(>|z|)`)
    dusted.vari.table=pixiedust::dust(vari.table)%>%
      pixiedust::sprinkle_na_string(na_string = "")%>%
      pixiedust::sprinkle_print_method(print_method = "html")%>%
      pixiedust::sprinkle_colnames("Variance","Estimate","Std.Err","Std.lv","Std.all","z-value","P(>|z|)")%>%
      pixiedust::sprinkle_border(border="all")%>%
      pixiedust::sprinkle_border(border="all",part="head")%>%
      pixiedust::sprinkle_pad(pad=7)%>%
      pixiedust::sprinkle_pad(pad=7,part="head")%>%
      pixiedust::sprinkle(cols=7,fn = quote(pvalString(value)))%>%
      pixiedust::sprinkle_align(halign="center",part="head")
    #dusted.vari.table
  }else{
    dusted.vari.table=NULL
  }

  dusted.r2=pixiedust::dust(my.r2)%>%
    pixiedust::sprinkle_na_string(na_string = "")%>%
    pixiedust::sprinkle_print_method(print_method = "html")%>%
    pixiedust::sprinkle_border(border="all")%>%
    pixiedust::sprinkle_border(border="all",part="head")%>%
    pixiedust::sprinkle_pad(pad=7)%>%
    pixiedust::sprinkle_pad(pad=7,part="head")%>%
    pixiedust::sprinkle_align(halign="center",part="head")%>%
    pixiedust::sprinkle_colnames("Variable","R^2")
  #dusted.r2

  # mydustlist=list(dusted.fit.table,dusted.latent.table)
  # if(!is.null(dusted.reg.table)){
  #   mydustlist=list(mydustlist,dusted.reg.table)
  # }
  # if(!is.null(dusted.covar.table)){
  #   mydustlist=list(mydustlist,dusted.covar.table)
  # }
  # if(!is.null(dusted.vari.table)){
  #   mydustlist=list(mydustlist,dusted.vari.table)
  # }
  # mydustlist=list(mydustlist,dusted.r2)
  mydustlist=list(dusted.fit.table,dusted.latent.table, dusted.reg.table, dusted.covar.table, dusted.vari.table,dusted.r2)
  #print(mydustlist[[1]])
  dust.num=1
  while(dust.num<{length(mydustlist)+1}){
    print(mydustlist[[dust.num]])
    x=readline(prompt = "Press z and enter to go back, or enter to go forward\n")
    if(x!="z"){
      dust.num=dust.num+1
    }else{
      if(dust.num>1){
        dust.num=dust.num-1
      }
    }
  }
  options(width=prev.width)
}

