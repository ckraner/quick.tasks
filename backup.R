#' Regression Tables in Pixiedust
#'
#' Beautiful tables adding sums of squares
#' and p-value formatting, then giving html or
#' latex output. If interactive and HTML, will show up
#' in viewer.
#'
#' @param my.model Model. NOTE: If have factors, please place them first in the regression.
#' @param VIF include Variable Inflation Factor? (calculated by car::vif) (default=F)
#' @param part.eta If lm, include partial eta square by calculating SS_part/SS_total? (default=F)
#' @param myDF Dataframe, not needed if use data= in call
#' @param my.factor If there are any factors, list them here.
#' @param SS.type Type of sums of squares to report (default = 3)
#' @param pix.int Should this be viewable or is this for a document/dashboard? (default=T)
#' @param pix.method Print method. (default="html")
#' @param type Type of regression? Currently supported: lm, glm (binary), manova
#' @param do.glance Include glance stats?
#' @return Either pixiedust object or code (in HTML or LaTeX) for table
#' @keywords Explore
#' @export
#' @examples
#' quick.reg(my.model, myDF)

quick.reg = function(my.model,
                     part.eta = F,
                     VIF = F,
                     myDF = my.found.df,
                     SS.type = 3,
                     abbrev.length = ab.len,
                     pix.int = T,
                     pix.method = "html",
                     type = my.reg.type,
                     test.stat = "Pillai",
                     my.factor = NULL,
                     do.glance=T,
                     show.footer=T,
                     adjustment = "bonferroni",
                     show.contrasts=F,
                     show.y.contrasts=F,
                     show.latent=F,
                     show.intercepts=F,
                     real.names=T,
                     do.return=T) {
  library(pixiedust)
  library(broom)
  library(car)
  library(tidyr)
  library(phia)
  library(quick.tasks)
  library(dplyr)

  #### Find type ####
  my.reg.type=quick.type(my.model)

  #### Set Inits ####
  if(type=="ord"){
    ab.len=30
    library(ordinal)
  }else{
    ab.len=15
  }
  #### Get data frame from parent environment ####
  my.found.df = eval(parse(text=capture.output(my.model$call$data)),envir = parent.frame())
  #print(dim(my.found.df))
  if (is.null(my.found.df)) {
    stop(paste("No data frame found"))
  }

  #### for expansion
  if (!VIF & !part.eta) {
    v.p.len = 8
    v.p.rep = 0
  } else if (!VIF) {
    v.p.len = 9
    v.p.rep = 1
  } else if (!part.eta) {
    v.p.len = 9
    v.p.rep = 1
  } else{
    v.p.len = 10
    v.p.rep = 2
  }

  #### Make factor list ####
  #### Use car::Anova to get SS Type 3
  if (type == "lm") {
    #### ANOVA TABLES ####
    my.summary = summary(my.model)
    my.coefficients = my.summary$coefficients
    my.coefficients = as.data.frame(my.coefficients)
    if (type == "ord" & length(my.model$model) == 1) {
      my.III.summary = NULL
      the.length = length(my.model$y.levels) - 1
    } else{
      my.III.summary = car::Anova(my.model, type = SS.type)

      if (type == "glm" & is.null(my.factor)) {
        the.length = dim(my.III.summary)[1] + 1
      } else{
        the.length = dim(my.III.summary)[1]
      }
    }

    #### Calculate total SS ####
    my.total = sum(my.III.summary$`Sum Sq`[2:length(my.III.summary$`Sum Sq`)])
    my.df.total = sum(my.III.summary$Df[2:length(my.III.summary$Df)])
    total.intercepts = 1
    my.rownames = c(abbreviate(rownames(my.summary$coefficients), minlength = abbrev.length),
                    "Residuals",
                    "Total")

    treat.SS=sum(my.III.summary$`Sum Sq`[2:{length(my.III.summary$`Sum Sq`)-1}])
    treat.df=sum(my.III.summary$Df[2:{length(my.III.summary$Df)-1}])
  } else if(type == "manova"){


    #### Begin MANOVA ####
    #### Inits ####


    my.envir=environment()

    my.nested.table=quick.SSCP(my.model, myDF, SS.type, show.contrasts, show.latent,my.envir)
    #### Get treatment ####
    treat.model=my.nested.table[dim(my.nested.table)[1],4]

    #### Regular totals
    treat.total.temp=lapply(lapply(treat.model[[1]][[1]],diag),sum)
    treat.total=0
    for(W in 2:length(treat.total.temp)){
      treat.total=treat.total+treat.total.temp[[W]]
    }
    treat.total.df=sum(as.numeric(my.nested.table[,7]),na.rm = T)


    #### Partial totals
    part.treat.total=NULL
    for(X in 2:length(treat.model[[1]][[1]])){
      if(X==2){
        part.treat.total=treat.model[[1]][[1]][[X]]
      }else{
        part.treat.total=part.treat.total+treat.model[[1]][[1]][[X]]
      }
    }

    #### Latent totals
    if(show.latent){
      latent.treat.model=my.nested.table[dim(my.nested.table)[1],8]
      latent.part.treat.total=NULL
      latent.part.total=NULL
      for(X in 2:length(latent.treat.model[[1]][[1]])){
        if(X==2){
          latent.part.treat.total=latent.treat.model[[1]][[1]][[X]]
        }else{
          latent.part.treat.total=latent.part.treat.total+latent.treat.model[[1]][[1]][[X]]
        }
      }
    }



    #### Get residuals ####
    #### Regular residuals
    if(SS.type==3){
      the.resid=sum(diag(treat.model[[1]][[2]]))
      #### Partial residuals
      part.resid.total=treat.model[[1]][[2]]
    }else{
      the.resid=sum(diag(treat.model[[1]][[2]][[1]]))
      part.resid.total=treat.model[[1]][[2]][[1]]
    }
    the.resid.df=my.model$df.residual





    #### Latent totals
    if(show.latent){
      if(SS.type==3){
        latent.part.resid.total=latent.treat.model[[1]][[2]]
      }else{
        latent.part.resid.total=latent.treat.model[[1]][[2]][[1]]
      }
    }

    #### Get totals ####
    the.total=total.resid+treat.total+sum(diag(treat.model[[1]][[1]][[1]]))
    the.total.df=total.resid.df+treat.df+my.y.levels

    #### Partial totals
    partial.total=part.resid.total+part.treat.total+treat.model[[1]][[1]][[1]]

    #### Latent totals
    if(show.latent){
      latent.part.total=latent.part.treat.total+latent.part.treat.total+latent.treat.model[[1]][[1]][[1]]
    }



    #### Make dependent variable rownames ####
    my.dv.rownames=rownames(part.resid.total)

    #### Make basic table ####

    #### NEEDS TO BE FIXED ####
    my.table.names=c("var","test.stat","f.val","SS","df","mult.df","resid df","p.val")


    my.manova.table=as.data.frame(matrix(ncol=v.p.len,nrow=1))
    names(my.manova.table)=my.table.names
    my.line.var=1
    for(i in 1:length(my.SSP.treat)){
      if(i==1){
        my.i=i
      }else{
        my.i=2*i-1
      }
      #### Put in basic line ####
      my.treat.err=solve(my.SSP.err)%*%treat.model[[1]][[1]][[i]]
      my.test.stat=quick.m.test(my.treat.err,test.stat)
      my.SS=sum(diag(treat.model[[1]][[1]][[i]]))
      my.df=my.y.levels*ifelse(my.i!=1,as.numeric(my.nested.table[my.i,7]),1)
      my.resid.df=ifelse(my.i==1,{my.y.levels*the.resid.df},min({the.resid.df*as.numeric(my.nested.table[my.i,5])-as.numeric(my.nested.table[my.i,5])},my.y.levels*the.resid.df-as.numeric(my.nested.table[my.i,5])))
      my.f.val={my.SS/my.df}/{the.resid/the.resid.df}
      my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

      my.manova.table[my.line.var,]=c(ifelse(i==1,"(Intercept)",my.nested.table[my.i,1]),my.test.stat,my.f.val,my.SS,ifelse(my.i==1,NA,my.df/my.y.levels),ifelse(my.i==1,my.y.levels,my.df),my.resid.df,my.p.val)
      my.line.var=my.line.var+1


      #### Put in decomposed factors (y constrasts) ####
      #### Intercept
      if({show.intercepts & my.i==1}){
        for(y in 1:my.y.levels){
          my.name=my.dv.rownames[y]
          my.test.stat=NA
          my.SS=as.numeric(treat.model[[1]][[1]][[i]][y,y])
          my.df=1
          my.resid.df=the.resid.df
          my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
          my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

          my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
          my.line.var=my.line.var+1
        }
      }
      #### Other
      if({show.y.contrasts & my.i!=1}){
        for(y in 1:my.y.levels){
          my.name=paste(ifelse(real.names,my.dv.rownames[y],y),"|",my.nested.table[my.i,1],sep="")
          my.test.stat=NA
          my.SS=treat.model[[1]][[1]][[i]][y,y]
          my.df=as.numeric(my.nested.table[my.i,7])
          my.resid.df=min(abs(my.y.levels*the.resid.df-my.df*the.resid.df-my.df),abs(my.y.levels*the.resid.df-my.df))
          my.resid=part.resid.total[y,y]
          my.f.val={my.SS/my.df}/{my.resid/my.resid.df}
          my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

          my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
          my.line.var=my.line.var+1


          #### Put in latents (ANOVA) ####
          my.counter=NULL
          for(r in 2:my.y.levels){
            my.counter=c(my.counter,r)
          }
          my.counter=c(my.counter,1)
          if(show.latent &my.i!=1){
            # my.my.i.temp=my.counter[my.i-1]
            my.my.i.temp=my.i
            my.y=y
            my.name=NA
            my.df=as.numeric(my.nested.table[my.i,7])
            my.test.stat=NA
            my.resid.df=the.resid.df-my.df
            my.SS=as.numeric(latent.treat.model[[1]][[1]][[i]][y,1])
            my.f.val={my.SS/my.df}/{my.latent.SSP.err[my.y,my.y]/my.resid.df}
            my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

            my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
            my.line.var=my.line.var+1
          }

          #### Put in latent contrasts ####
          if(show.contrasts & show.latent & my.i!=1){
            #### Check length
            if({!is.na(my.nested.table[my.i,11])}){
              #other.manova.grep=grep(paste("^",names(my.SSP.treat)[my.i],"$",sep=""),names(my.model$xlevels))
              for(k in 1:{my.nested.table[my.i,7]}){
                my.name=my.nested.table[my.i,11][[1]][k,1]
                my.f.val=as.numeric(my.nested.table[my.i,11][[1]][k,2])
                my.SS=as.numeric(my.nested.table[my.i,11][[1]][k,3])
                my.test.stat=NA
                my.df=1
                my.mult.df=NA
                my.resid.df=the.resid.df-my.y.levels+1
                my.p.val=as.numeric(my.nested.table[my.i,11][[1]][k,4])

                my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,my.p.val)
                my.line.var=my.line.var+1
              }
            }
          }

          #### Put in contrasts ####
          #### Not right. Don't have it decomposed this way.
          if(show.contrasts & !show.latent & my.i!=1){
            #### Check length
            if(!is.na(my.nested.table[my.i,9])){
              #other.manova.grep=grep(paste("^",names(my.SSP.treat)[my.i],"$",sep=""),names(my.model$xlevels))
              for(k in 1:as.numeric(my.nested.table[my.i,7])){
                my.name=my.nested.table[my.i,9][[1]][k,1]
                my.f.val=as.numeric(my.nested.table[my.i,9][[1]][k,2])
                my.SS=as.numeric(my.nested.table[my.i,9][[1]][k,3])
                my.test.stat=NA
                my.df=NA
                my.mult.df=my.y.levels
                my.resid.df=the.resid.df-my.y.levels+1
                my.p.val=as.numeric(my.nested.table[my.i,9][[1]][k,4])

                my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,my.p.val)
                my.line.var=my.line.var+1
              }
            }
          }
        }
      }

      #### Put in contrasts ####
      if(show.contrasts & !show.latent & !show.y.contrasts){
        #### Check length
        if(!is.na(my.nested.table[my.i,8])){
          #other.manova.grep=grep(paste("^",names(my.SSP.treat)[my.i],"$",sep=""),names(my.model$xlevels))
          for(k in 1:as.numeric(my.nested.table[my.i,7])){
            my.name=my.nested.table[my.i,8][[1]][k,1]
            my.f.val=as.numeric(my.nested.table[my.i,8][[1]][k,2])
            my.SS=as.numeric(my.nested.table[my.i,8][[1]][k,3])
            my.test.stat=NA
            my.df=NA
            my.mult.df=my.y.levels
            my.resid.df=the.resid.df-my.y.levels+1
            my.p.val=as.numeric(my.nested.table[my.i,8][[1]][k,4])

            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,my.p.val)
            my.line.var=my.line.var+1
          }
        }
      }


      #### Put in latents (ANOVA) ####
      my.counter=NULL
      for(r in 2:my.y.levels){
        my.counter=c(my.counter,r)
      }
      my.counter=c(my.counter,1)
      if(show.latent &my.i!=1 & !show.y.contrasts){
        #my.my.i.temp=my.counter[my.i-1]
        my.my.i.temp=my.i
        for(y in 1:{my.y.levels}){
          my.y=y
          my.name=paste(ifelse(real.names,my.dv.rownames[y],y),"|",my.nested.table[my.i,1],sep="")
          my.df=as.numeric(my.nested.table[my.i,7])
          my.test.stat=NA
          my.resid.df=the.resid.df-my.df
          my.SS=as.numeric(my.nested.table[my.i,8][[1]][[1]][[i]][my.y])
          my.f.val={my.SS/my.df}/{my.latent.SSP.err[my.y,my.y]/my.resid.df}
          my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

          my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
          my.line.var=my.line.var+1

          #### Put in latent contrasts ####
          if(show.contrasts & !{show.contrasts & show.latent}){
            #### Check length
            if({my.SSP.treat.df[my.i]>1}){
              other.manova.grep=grep(paste("^",names(my.SSP.treat)[my.i],"$",sep=""),names(my.model$xlevels))
              for(k in 1:{my.SSP.treat.df[my.i]}){
                my.name=my.contrasts.table[[other.manova.grep]][k,1]
                my.f.val=my.latent.contrasts.table[[other.manova.grep]][[y]][k,2]
                my.SS=my.latent.contrasts.table[[other.manova.grep]][[y]][k,3]
                my.test.stat=NA
                my.df=1
                my.mult.df=NA
                my.resid.df=the.resid.df-my.y.levels+1
                my.p.val=my.latent.contrasts.table[[other.manova.grep]][[y]][k,4]

                my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,my.p.val)
                my.line.var=my.line.var+1
              }
            }
          }
        }
      }


      #### Put in treatment ####
      #### From Type my.i statistics
      if(my.i==1){
        my.test.stat=quick.m.test(part.treat.total,test.stat)
        my.SS=as.numeric(treat.total)
        my.df=my.y.levels*as.numeric(treat.total.df)
        #### Should really be a min statement, but for later...
        my.resid.df=my.y.levels*as.numeric(the.resid.df)
        my.f.val={my.SS/my.df}/{as.numeric(the.resid.SS)/the.resid.df}
        my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

        my.manova.table[my.line.var,]=c("Treatment Change",my.test.stat,my.f.val,my.SS,{my.df/my.y.levels},my.df,my.resid.df,my.p.val)
        my.line.var=my.line.var+1

        #### Treatment Y Contrasts ####
        if(show.y.contrasts){
          for(b in 1:my.y.levels){
            my.test.stat=NA
            my.name=paste(ifelse(real.names,my.dv.rownames[b],b),"|Treatment",sep="")
            my.SS=latent.part.treat.total[b]
            my.df=the.resid.df
            #### Should really be a min statement, but for later...
            my.resid.df=my.y.levels*the.resid.df
            my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
            my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,{my.df/my.y.levels},NA,{my.resid.df/my.y.levels},my.p.val)
            my.line.var=my.line.var+1

            #### Treatment Latents ####
            if(show.latent){
              #### Show latent treatments (ANOVAs)
              #### Need to change latents to type II
              my.SS=sum(weighted.residuals(my.null.model)[,b]^2)-sum(weighted.residuals(my.model)[,b]^2)
              my.df=as.numeric(treat.total.df)
              #### Should really be a min statement, but for later...
              my.resid.df=the.resid.df-my.df
              my.f.val={my.SS/my.df}/{my.latent.SSP.err[y,y]/my.resid.df}
              my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
              my.manova.table[my.line.var,]=c(NA,NA,my.f.val,my.SS,my.df,NA,{my.resid.df},my.p.val)
              my.line.var=my.line.var+1
            }
          }
        }

        #### Latents ####
        my.counter=NULL
        for(r in 2:my.y.levels){
          my.counter=c(my.counter,r)
        }
        my.counter=c(my.counter,1)
        if(show.latent & !show.y.contrasts){
          for(b in 1:{my.y.levels}){
            my.name=paste(ifelse(real.names,my.dv.rownames[b],b),"|Treatment",sep="")
            my.SS=as.numeric(part.treat.total[b])
            my.df=as.numeric(treat.total.df)
            #### Should really be a min statement, but for later...
            my.resid.df=the.resid.df-my.df
            my.f.val={my.SS/my.df}/{as.numeric(latent.part.resid.total[b])/the.resid.df}
            my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),NA,my.f.val,my.SS,my.df,NA,{my.resid.df},my.p.val)
            my.line.var=my.line.var+1
          }
        }
      }
    }

    #### Put in residuals, total change, total ss ####
    my.manova.table[my.line.var,]=c("Total Change",NA,NA,as.numeric(treat.total),treat.total.df,{my.y.levels*treat.total.df},NA,NA,rep(NA,v.p.rep))
    my.line.var=my.line.var+1
    my.manova.table[my.line.var,]=c("Residuals",NA,NA,the.resid,the.resid.df,my.y.levels*the.resid.df,NA,NA,rep(NA,v.p.rep))
    my.line.var=my.line.var+1
    if(show.y.contrasts){
      for(i in 1:my.y.levels){
        my.manova.table[my.line.var,]=c(paste(ifelse(real.names,my.dv.rownames[i],i),"|Residual",sep=""),NA,NA,part.resid.total[i,i],the.resid.df,NA,NA,NA)
        my.line.var=my.line.var+1
        if(show.latent){
          my.manova.table[my.line.var,]=c(NA,NA,NA,latent.part.resid.total[i],the.resid.df,NA,NA,NA)
          my.line.var=my.line.var+1
        }
      }
    }
    if(show.latent & !show.y.contrasts){
      for(i in 1:my.y.levels){
        my.manova.table[my.line.var,]=c(paste(ifelse(real.names,my.dv.rownames[i],i),"|Residual",sep=""),NA,NA,latent.part.resid.total[i],the.resid.df,NA,NA,NA)
        my.line.var=my.line.var+1
      }
    }
    my.manova.table[my.line.var,]=c("Total",NA,NA,the.total,the.total.change.df+1,my.y.levels*the.total.change.df+2,NA,NA,rep(NA,v.p.rep))
    if(show.footer){
      the.footer=paste(ifelse(dim(my.new.df)[1]==dim(myDF)[1],"Data have same number of rows <br />",paste({dim(myDF)[1]-dim(my.new.df)[1]}," cases deleted due to missingness <br />")),"Method: QR decomposition",if(show.contrasts){paste(" <br />Adjustment: ", adjustment,sep="")},if(show.latent){paste(" <br /> Latent Contrasts")})
    }else{
      the.footer=NULL
    }
    my.html.table=quick.table(my.manova.table,"manova",test=test.stat,SS.type = SS.type, abbrev.length = abbrev.length,the.footer = the.footer)
    if(do.return){
      return(my.html.table)
    }else{}
    #### ETA-SQ Stuff to finish.... ####
    # #### Make eta-sq
    # my.SSP.err.t=t(my.SSP.err)
    #
    # eta.sq=NULL
    # for(i in 1:length(my.SSP)){
    #   eta.sq[[i]]=my.SSP.err.t*my.SSP[[i]]
    # }
    #
    #       x3 = capture.output(car::Anova(my.model, type = SS.type, test = test.stat))
    #       my.manova.test = data.frame(matrix(ncol = 7, nrow = 1))
    #       my.var.temp = 4
    #       while (my.var.temp < {
    #         length(x3) - 1
    #       }) {
    #         test = strsplit(x3[my.var.temp], "\\s+")
    #
    #         if (length(test[[1]]) == 9) {
    #           test2 = test[[1]][-9]
    #           test2 = test2[-7]
    #
    #         } else if (length(test[[1]]) == 8) {
    #           test2 = test[[1]][-8]
    #
    #         } else{
    #           test2 = test[[1]]
    #         }
    #
    #         my.manova.test[{
    #           my.var.temp - 3
    #         }, ] = test2
    #         my.var.temp = my.var.temp + 1
    #
    #       }
    #       my.summary=summary(my.model)

    #### END MANOVA table ####
  }else if (type == "glm") {

    #### DevAn Table ####

    #### Inits ####
    new.df=my.model$model
    new.model=update(my.model,data=new.df)
    null.model=update(new.model,~1)

    total.intercepts = 1

    #### Check if null model ####
    vars.df=new.model$df.null-new.model$df.residual
    if(vars.df==0){
      stop("This is the null model.")
    }

    resid.dev=new.model$deviance
    resid.df=new.model$df.residual
    vars.dev=anova(new.model,test="Chi")$Deviance[-1]
    vars.dev.df=anova(new.model,test="Chi")$Df[-1]
    vars.dev.p=anova(new.model,test="Chi")$`Pr(>Chi)`[-1]
    total.dev=new.model$null.deviance-new.model$deviance
    vars.dev.total=sum(vars.dev)
    if(length(grep("*",new.model$formula))>0){
      my.full.model=new.model
    }else{
      my.vars=strsplit(as.character(new.model$formula),"~")
      my.dep.var=my.vars[[2]]
      my.other.vars=names(new.model$model)
      my.var.grep=grep(paste("^",my.dep.var,"$",sep=""),my.other.vars)
      my.other.vars=my.other.vars[-my.var.grep]
      my.formula=paste("~",my.other.vars[1])
      for(i in 2:length(my.other.vars)){
        my.formula=paste(my.formula,"*",my.other.vars[i],sep="")}

      #### FULL MODEL THING ####
      my.full.model=update(new.model,my.formula)
    }



    my.full.model.z=summary(my.full.model)$coefficients[1:total.intercepts,3]^2
    my.null.model.z=summary(null.model)$coefficients[1:total.intercepts,3]^2
    my.full.dev=sum(my.full.model.z)
    my.null.dev=sum(my.null.model.z)
    my.int.dev.total=abs(total.intercepts*{my.full.dev-my.null.dev})
    my.int.dev=summary(new.model)$coefficients[1:total.intercepts,2]^2

    treat.dev=anova(null.model,new.model)$`Deviance`[2]

    #vars.dev.df=drop1(new.model,test="Chi")$Df[-1]
    #vars.dev.p=drop1(new.model,test="Chi")$`Pr(>Chi)`[-1]
    #total.dev=-2*{as.integer(levels(null.model$info$logLik)[1])-as.integer(levels(new.model$info$logLik)[1])}
    total.dev.change=treat.dev+my.int.dev.total
    total.dev.change.df=total.intercepts*vars.df


    total.dev=new.model$null.deviance
    resid.dev=total.dev-total.dev.change

    vars.dev.df=anova(new.model,test="Chi")$Df
    resid.df=new.model$df.residual
    total.df=resid.df+vars.df




    my.int.dev.df=max({total.intercepts*sum(vars.dev.df,na.rm = T)-sum(vars.dev.df,na.rm=T)},1)
    my.int.dev.or=exp(coef(new.model))[1:{total.intercepts}]
    my.int.dev.or.confint=exp(confint(new.model,type="Wald"))[1:total.intercepts,]
    my.int.names=abbreviate(names(my.int.dev.or),minlength=abbrev.length)


    #resid.dev=new.model$deviance

    vars.dev=anova(new.model,test="Chi")$Deviance

    vars.dev.p=anova(new.model,test="Chi")$`Pr(>Chi)`
    vars.names=rownames(anova(new.model,test="Chi"))
    #total.dev=new.model$null.deviance-new.model$deviance
    vars.dev.total=sum(vars.dev)
    #total.df=sum(vars.dev.df,na.rm = T)

    vars.or=exp(coef(new.model))
    vars.or.confint=exp(confint(new.model,type="Wald"))
    my.names=abbreviate(names(vars.or),minlength=abbrev.length)
    #vars.or=exp(coef(new.model))[-1]
    #vars.or.confint=exp(confint(new.model,type="Wald"))[-1,]




  } else if (type == "ord" | type =="glm2") {
    new.df=my.model$model
    new.model=update(my.model,data=new.df)
    null.model=update(new.model,~1)

    total.intercepts = null.model$edf
    my.summary=summary(new.model)

    if(length(grep("*",new.model$formula))>0){
      my.full.model=new.model
    }else{
      my.vars=strsplit(as.character(new.model$formula),"~")
      my.dep.var=my.vars[[2]]
      my.other.vars=names(new.model$model)
      if(length(my.other.vars)>1){
        my.var.grep=grep(paste("^",my.dep.var,"$",sep=""),my.other.vars)
        my.other.vars=my.other.vars[-my.var.grep]
        my.formula=paste("~",my.other.vars[1])
        for(i in 2:length(my.other.vars)){
          my.formula=paste(my.formula,"*",my.other.vars[i],sep="")}

        #### FULL MODEL THING ####
        my.full.model=update(new.model,my.formula)
      }else{
        my.full.model=new.model
      }
    }

    my.full.model.z=summary(my.full.model)$coefficients[1:total.intercepts,3]^2
    my.null.model.z=summary(null.model)$coefficients[1:total.intercepts,3]^2
    my.full.dev=sum(my.full.model.z)
    my.null.dev=sum(my.null.model.z)
    my.int.dev.total=abs(total.intercepts*{my.full.dev-my.null.dev})
    my.int.dev=summary(new.model)$coefficients[1:total.intercepts,2]^2

    vars.df=new.model$edf-null.model$edf

    if(vars.df==0){
      stop("This is the null model.")
    }

    vars.dev=summary(new.model)$coefficients[{total.intercepts+1}:dim(summary(new.model)$coefficients)[1],3]^2
    vars.dev.p=summary(new.model)$coefficients[{total.intercepts+1}:dim(summary(new.model)$coefficients)[1],4]

    vars.dev.df=NULL
    weird.var=2
    track.var=1
    while(weird.var<=dim(new.model$model)[2]){
      if(is.factor(new.model$model[[weird.var]])){
        vars.dev.df[track.var]=length(levels(new.model$model[[weird.var]]))-1
      }else{
        vars.dev.df[track.var]=1
      }
      weird.var=weird.var+1
      track.var=track.var+1
    }

    ### interaction effects
    weird.var=weird.var+total.intercepts-1
    while(weird.var<=dim(summary(new.model)$coefficients)[1]){
      vars.dev.df[track.var]=1
      weird.var=weird.var+1
      track.var=track.var+1
    }

    vars.df.total=sum(vars.dev.df,na.rm = T)
    treat.dev=anova(null.model,new.model)$`LR.stat`[2]

    #vars.dev.df=drop1(new.model,test="Chi")$Df[-1]
    #vars.dev.p=drop1(new.model,test="Chi")$`Pr(>Chi)`[-1]
    #total.dev=-2*{as.integer(levels(null.model$info$logLik)[1])-as.integer(levels(new.model$info$logLik)[1])}
    total.dev.change=treat.dev+my.int.dev.total
    if(total.intercepts>1){
      total.dev.change.df=total.intercepts*vars.df.total
    }else{
      total.dev.change.df=vars.df.total+total.intercepts
    }

    resid.dev=-2*new.model$logLik
    total.dev=total.dev.change+resid.dev


    total.df=dim(myDF)[1]
    resid.df=total.df-total.dev.change.df

    vars.or=exp(coef(new.model))[{null.model$edf+1}:length(coef(new.model))]
    vars.or.confint=exp(confint(new.model,type="Wald"))[{null.model$edf+1}:length(coef(new.model)),]
    my.names=abbreviate(names(vars.or),minlength=abbrev.length)


    my.int.dev.df=max({total.intercepts*sum(vars.df,na.rm = T)-sum(vars.df,na.rm=T)},1)
    my.int.dev.or=exp(coef(new.model))[1:{null.model$edf}]
    my.int.dev.or.confint=exp(confint(new.model,type="Wald"))[1:null.model$edf,]
    my.int.names=abbreviate(names(my.int.dev.or),minlength=abbrev.length)


  } else{
    print("Error")
    return()
  }


  if (!VIF & !part.eta) {
    v.p.len = 7
    v.p.rep = 0
  } else if (!VIF) {
    v.p.len = 8
    v.p.rep = 1
  } else if (!part.eta) {
    v.p.len = 8
    v.p.rep = 1
  } else{
    v.p.len = 9
    v.p.rep = 2
  }

  #### Make table if not MANOVA ####
  my.tables.df = as.data.frame(matrix(ncol = v.p.len, nrow = 1))

  if (type == "lm") {
    if (!VIF & !part.eta) {
      my.std.error = c(my.coefficients$`Std. Error`, NA)
      my.estimate = c(my.coefficients$Estimate, NA)
      names(my.tables.df) = c("rownames",
                              "sumsq",
                              "df",
                              "est",
                              "std.err",
                              "f.val",
                              "p.val")
    } else if (!part.eta) {
      my.VIF = car::vif(my.model)
      names(my.tables.df) = c("rownames",
                              "sumsq",
                              "df",
                              "est",
                              "std.err",
                              "f.val",
                              "p.val",
                              "VIF")
    } else if (!VIF) {
      names(my.tables.df) = c("rownames",
                              "sumsq",
                              "df",
                              "est",
                              "std.err",
                              "f.val",
                              "p.val",
                              "p.eta")
    } else{
      my.VIF = car::vif(my.model)
      names(my.tables.df) = c("rownames",
                              "sumsq",
                              "df",
                              "est",
                              "std.err",
                              "f.val",
                              "p.val",
                              "p.eta",
                              "VIF")
    }
  } else if (type == "glm" & VIF) {
    my.VIF = car::vif(my.model)
    names(my.tables.df) = c("var",
                            "p.odd",
                            "p.odd.2.5",
                            "p.odd.97.5",
                            "deviance",
                            "df",
                            "p.val",
                            "VIF")
  } else{
    names(my.tables.df) = c("var", "p.odd", "p.odd.2.5", "p.odd.97.5", "deviance", "df", "p.val")
  }

  #### Make the double table entries ####

  #### Was very annoying...took my frustration out on
  #### Variable names

  factor.stupid = NULL
  factor.rownames = NULL
  num.of.levels = NULL
  ord.temp = 0
  if (!is.null(my.factor)) {
    for (i in 1:length(my.factor)) {
      factor.stupid = c(factor.stupid, grep(paste("^", my.factor[i], "$", sep =
                                                    ""), names(myDF)))

      if (type == "lm") {
        factor.rownames = c(factor.rownames, grep(
          paste("^", my.factor[i], "$", sep = ""),
          rownames(my.III.summary)
        ))
        num.of.levels = c(num.of.levels, length(levels(myDF[[factor.stupid[i]]])))
        # }else if(type=="glm"){
        #   factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+1})
        #   num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
      } else{
        my.dev.grep=grep(paste("^", my.factor[i], sep = ""),names(new.model$coefficients))[1]
        factor.rownames=c(factor.rownames,my.dev.grep)
        num.of.levels = c(num.of.levels, length(levels(myDF[[factor.stupid[i]]])))

        # factor.rownames = c(factor.rownames, {
        #   grep(paste("^", my.factor[i], "$", sep = ""),
        #        rownames(my.III.summary)) + total.intercepts + sum(num.of.levels) - ordinal.temp
        # })
        # num.of.levels = c(num.of.levels, length(levels(myDF[[factor.stupid[i]]])))
        # ordinal.temp = ordinal.temp + 1
      }
    }

  } else{
    factor.rownames = 0L

  }



  #### Make phia stuff ####
  if(show.contrasts){
    if(type=="glm" & !is.null(my.factor)){
      my.phia.reg=quick.contrast(new.model,skip.me=T,adjustment = adjustment,SS.type = SS.type,abbrev.length = abbrev.length,my.factors = my.factor)
      my.phia.rownames=my.phia.reg$names
      phia.dev=my.phia.reg$Chisq

    }else if(type=="lm" & !is.null(my.factor)){
      my.phia.reg=quick.contrast(my.model,skip.me=T,adjustment = adjustment,SS.type = SS.type,abbrev.length = abbrev.length)
      my.j=0
      my.big.phia=NULL
      phia.shift=0
      real.shift=0
      for(i in 1:dim(my.phia.reg)[1]){
        if(!is.na(my.phia.reg[i,1])){
          my.j=my.j+1
          real.shift=real.shift+phia.shift
          my.big.phia[[my.j]]=my.phia.reg[i,2:7]
        }else{
          my.big.phia[[my.j]][i-real.shift,]=my.phia.reg[i,2:7]
        }
        phia.shift=phia.shift+1
      }

      my.phia.rownames=NULL
      my.phia.SS=NULL
      my.phia.value=NULL
      my.phia.F=NULL
      my.phia.p=NULL
      my.phia.err=NULL
      for(i in 1:{my.j}){
        my.phia.rownames[[i]]=my.big.phia[[i]]$names
        my.phia.SS[[i]]=my.big.phia[[i]]$`Sum of Sq`
        my.phia.value[[i]]=my.big.phia[[i]]$Value
        my.phia.F[[i]]=my.big.phia[[i]][,5]
        my.phia.p[[i]]=my.big.phia[[i]][,6]

      }
    }else{

    }
  }



  my.factor.var = 1
  this.temp.var = 1
  this.shift.temp = 1
  yet.another.var = 1
  my.shift = 0
  other.temp = 2
  ord.temp = 0
  phia.temp=1
  if (type == "lm") {
    dang.length = length(rownames(my.III.summary))
  } else if (type == "glm") {
    # dang.length = total.intercepts + length(rownames(my.III.summary)) + max({
    #   sum(num.of.levels) - length(my.factor)
    # }, 0) + 1
    dang.length=length(new.model$coefficients)-sum(num.of.levels)+length(num.of.levels)+total.intercepts+1
  } else{
    # dang.length = total.intercepts + length(rownames(my.III.summary)) + max({
    #   sum(num.of.levels) - length(my.factor)
    # }, 0) + 1
    #dang.length=7
    dang.length=length(new.model$coefficients)+length(my.factor)+1
  }



  while (this.shift.temp < dang.length) {
    if (is.na(factor.rownames[my.factor.var])) {
      my.factor.rownames = 1

    } else{
      my.factor.rownames = factor.rownames[my.factor.var]

    }
    #### LOOP ####
    if (this.shift.temp == 1) {
      i = 1
      if(type=="ord" | type=="glm"){
        my.tables.df[this.temp.var, ] = c(
          "Intercept Change",
          NA,
          NA,
          NA,
          my.int.dev.total,
          my.int.dev.df,
          pchisq(my.int.dev.total,my.int.dev.df,lower.tail = F),
          rep(NA, v.p.rep))

        this.temp.var=this.temp.var+1
      }
      while (i <= total.intercepts) {
        if (type == "lm") {
          my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
          my.df = my.III.summary$Df[this.shift.temp]
          my.est = my.estimate[this.shift.temp]
          my.std.err = my.std.error[this.shift.temp]
          my.f.val = my.III.summary$`F value`[this.shift.temp]
          my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
          my.tables.df[this.temp.var, ] = c(
            my.rownames[this.shift.temp],
            summary(my.model)[[4]][this.shift.temp,1],
            summary(my.model)[[4]][this.shift.temp,2],
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            rep(NA, v.p.rep)
          )
        } else if(type=="glm"){
          #### HERE IS WHERE I LEFT OFF!! ####
          # my.or = exp(my.estimate[this.shift.temp])
          # my.est = my.estimate[this.shift.temp]
          # my.z.val = my.summary$coefficients[this.shift.temp, 3]
          # my.std.err = my.std.error[this.shift.temp]
          #my.dev=my.III.summary$`LR Chisq`[this.shift.temp]
          #my.df=my.III.summary$Df[this.shift.temp]
          # my.p.val = my.summary$coefficients[{
          #   this.shift.temp
          #}, 4]

          my.tables.df[this.temp.var, ] = c(
            paste(names(attr(my.model$model[[i]],"labels"))[1],"-",names(attr(my.model$model[[i]],"labels"))[2],sep=""),
            my.int.dev.or[i],
            my.int.dev.or.confint[1],
            my.int.dev.or.confint[2],
            my.int.dev[i],
            1,
            NA,
            rep(NA, v.p.rep))

          #this.temp.var=this.temp.var+1
          #my.tables.df[this.temp.var,]=c("Treatment",NA,NA,NA,total.dev,total.df,dchisq(total.dev,total.df),rep(NA,v.p.rep))

        }else{
          my.tables.df[this.temp.var, ] = c(
            my.int.names[i],
            my.int.dev.or[i],
            my.int.dev.or.confint[i,1],
            my.int.dev.or.confint[i,2],
            my.int.dev[i],
            NA,
            NA,
            rep(NA, v.p.rep))

          #this.temp.var=this.temp.var+1

        }

        this.shift.temp = this.shift.temp + 1

        this.temp.var = this.temp.var + 1
        i = i + 1

      }
      if(type=="ord" | type=="glm"){
        my.tables.df[this.temp.var,]=c("Treatment Change",NA,NA,NA,treat.dev,vars.df,pchisq(treat.dev,vars.df,lower.tail = F))
        this.temp.var=this.temp.var+1
      }else{
        my.tables.df[this.temp.var,]=c("Treatment",NA,NA,treat.SS,treat.df,glance(my.model)[4],glance(my.model)[5],rep(NA,v.p.rep))
        this.temp.var=this.temp.var+1
      }
    } else if (this.shift.temp == my.factor.rownames) {
      if (type == "lm") {
        my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
        my.df = my.III.summary$Df[this.shift.temp]
        my.est = my.estimate[this.shift.temp]
        my.std.err = my.std.error[this.shift.temp]
        my.z.val = NA
        my.f.val = my.III.summary$`F value`[this.shift.temp]
        my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
        if(!VIF & !part.eta){
          my.tables.df[this.temp.var, ] = c(
            my.factor[yet.another.var],
            NA,
            NA,
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            rep(NA, v.p.rep)
          )
        }else if(part.eta & !VIF){
          my.p.eta = my.sumsq / my.total
          my.tables.df[this.temp.var, ] = c(
            my.factor[yet.another.var],
            NA,
            NA,
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.p.eta
          )
        }else if(!part.eta & VIF){

          my.tables.df[this.temp.var, ] = c(
            my.factor[yet.another.var],
            NA,
            NA,
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.VIF[this.shift.temp - 1, 1]
          )
        }else{
          my.p.eta = my.sumsq / my.total
          my.tables.df[this.temp.var, ] = c(
            my.factor[yet.another.var],
            NA,
            NA,
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.p.eta,
            my.VIF[this.shift.temp - 1, 1]
          )
        }
      }else{

        # my.or = NA
        # my.est = NA
        # my.std.err = NA
        # my.dev = my.III.summary$`LR Chisq`[ord.temp]
        # my.df = my.III.summary$Df[ord.temp]
        # my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
        my.tables.df[this.temp.var,]=c(my.factor[this.shift.temp-1],
                                       NA,
                                       NA,
                                       NA,
                                       drop1(new.model,test="Chi")$`LRT`[this.shift.temp],
                                       drop1(new.model,test="Chi")$`Df`[this.shift.temp],
                                       drop1(new.model,test="Chi")$`Pr(>Chi)`[this.shift.temp],
                                       rep(NA,v.p.rep))
        #
        # my.tables.df[this.temp.var, ] = c(my.factor[yet.another.var],
        #                                   my.or,
        #                                   my.est,
        #                                   my.std.err,
        #                                   my.dev,
        #                                   my.df,
        #                                   my.p.val)
        #ord.temp = ord.temp + 1
      }
      # else{
      #   my.or = NA
      #   my.est = NA
      #   my.std.err = NA
      #   my.dev = my.III.summary$Chisq[ord.temp]
      #   my.df = my.III.summary$Df[ord.temp]
      #   my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
      #   my.tables.df[this.temp.var, ] = c(my.factor[yet.another.var],
      #                                     my.or,
      #                                     my.est,
      #                                     my.std.err,
      #                                     my.dev,
      #                                     my.df,
      #                                     my.p.val)
      #   ord.temp = ord.temp + 1
      #
      # }
      yet.another.var = yet.another.var + 1
      this.temp.var = this.temp.var + 1
      this.shift.temp = this.shift.temp + 1
      if(type=="lm"){
        other.other.temp = 2
      }else{
        other.other.temp=1
      }

      #### INTERACTION EFFECTS? NOT WORRIED YET ####
      if ({length(grepl(":", my.summary$coefficients[other.temp, 1])) >
          0}) {

        #### I think the while should not be +1
        if(type=="glm" | type=="ord"){
          num.of.levels[my.factor.var]=num.of.levels[my.factor.var]-1
        }
        while (other.other.temp < {
          num.of.levels[my.factor.var] + 1
        }) {
          if (type == "lm") {
            if(show.contrasts){
              my.sumsq = NA
              my.df = NA
              my.est = my.estimate[other.temp]
              my.std.err = my.std.error[other.temp]
              #### NEED TO FIX ####
              my.f.val = {
                my.summary$coefficients[other.temp, 3] ^ 2
              }
              my.p.val = my.summary$coefficients[other.temp, 4]
              my.tables.df[this.temp.var, ] = c(
                my.phia.rownames[[phia.temp]][other.temp-1],
                my.phia.value[[phia.temp]][other.temp-1],
                my.est,
                my.phia.SS[[phia.temp]][other.temp-1],
                1,
                my.phia.F[[phia.temp]][other.temp-1],
                my.phia.p[[phia.temp]][other.temp-1],
                rep(NA, v.p.rep)
              )
            }
            ord.temp=ord.temp + 1
          } else if (type == "glm") {
            # my.or = exp(my.estimate[other.temp])
            # my.est = my.estimate[other.temp]
            # my.std.err = my.std.error[other.temp]
            # my.z.val = my.summary$coefficients[other.temp, 3]
            # my.dev = my.III.summary$`LR Chisq`[other.temp]
            # my.df = my.III.summary$Df[other.temp]
            # my.p.val = my.summary$coefficients[other.temp, 4]
            # my.tables.df[this.temp.var, ] = c(my.rownames[other.temp],
            #                                   my.or,
            #                                   my.std.err,
            #                                   my.z.val,
            #                                   NA,
            #                                   NA,
            #                                   my.p.val)
            if(show.contrasts){
              my.tables.df[this.temp.var, ] = c(my.phia.rownames[other.temp-1],
                                                vars.or[other.temp-1],
                                                vars.or.confint[other.temp-1,1],
                                                vars.or.confint[other.temp-1,2],
                                                my.phia[other.temp-1,3],
                                                my.phia[other.temp-1,4],
                                                my.phia[other.temp-1,6],rep(NA,v.p.rep))
              #this.shift.temp = this.shift.temp + 1
            }
            ord.temp=ord.temp + 1
          } else{
            my.or = exp(my.estimate[other.temp + total.intercepts - 1])
            my.est = my.estimate[other.temp + total.intercepts - 1]
            my.std.err = my.std.error[other.temp + total.intercepts -
                                        1]
            my.z.val = my.summary$coefficients[{
              other.temp + total.intercepts - 1
            }, 3]
            #my.dev=my.III.summary$Chisq[other.temp]
            #my.df=my.III.summary$Df[other.temp]
            my.p.val = my.summary$coefficients[{
              other.temp + total.intercepts - 1
            }, 4]
            if (!VIF) {
              my.tables.df[this.temp.var, ] = c(my.rownames[other.temp +
                                                              total.intercepts - 1],
                                                my.or,
                                                my.std.err,
                                                my.z.val,
                                                NA,
                                                NA,
                                                my.p.val)
            } else{
              my.tables.df[this.temp.var, ] = c(my.rownames[other.temp +
                                                              total.intercepts - 1],
                                                my.or,
                                                my.std.err,
                                                my.z.val,
                                                NA,
                                                NA,
                                                my.p.val,
                                                my.VIF[this.shift.temp - 1, 1])

            }
            this.shift.temp = this.shift.temp + 1
          }
          this.temp.var = this.temp.var + 1
          other.temp = other.temp + 1
          other.other.temp = other.other.temp + 1
          the.length = the.length + 1

        }
        phia.temp=phia.temp+1
        if(!show.contrasts){this.temp.var=this.temp.var-ord.temp}
      } else{

      }

      if (my.factor.var == 1) {
        my.shift = {
          my.shift + other.temp - 2
        }

      } else{
        my.shift = my.shift + other.temp

      }

      my.factor.var = my.factor.var + 1

    } else{
      if (type == "lm") {
        my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
        my.df = my.III.summary$Df[this.shift.temp]
        my.est = my.estimate[this.shift.temp]
        my.std.err = my.std.error[this.shift.temp]
        my.f.val = my.III.summary$`F value`[this.shift.temp]
        my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
        if (!VIF & !part.eta) {
          my.tables.df[this.temp.var, ] = c(
            rownames(my.III.summary)[this.shift.temp],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val
          )
        } else if (!part.eta) {
          my.tables.df[this.temp.var, ] = c(
            rownames(my.III.summary)[this.shift.temp],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.VIF[this.shift.temp - 1, 1]
          )
        } else if (!VIF) {
          my.p.eta = my.sumsq / my.total
          my.tables.df[this.temp.var, ] = c(
            rownames(my.III.summary)[this.shift.temp],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.p.eta
          )
        } else{
          my.p.eta = my.sumsq / my.total
          my.tables.df[this.temp.var, ] = c(
            rownames(my.III.summary)[this.shift.temp],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
            summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
            my.sumsq,
            my.df,
            my.f.val,
            my.p.val,
            my.p.eta,
            my.VIF[this.shift.temp - 1, 1]
          )
        }
        ord.temp=ord.temp+1
      } else if (type == "glm") {
        # if (!is.null(my.factor)) {
        #   this.shift.temp = this.shift.temp - ord.temp-1
        # }
        # my.or = exp(my.estimate[this.shift.temp])
        # my.est = my.estimate[this.shift.temp]
        # my.std.err = my.std.error[this.shift.temp]
        # my.z.val = my.summary$coefficients[other.temp, 3]
        # my.dev = my.III.summary$`LR Chisq`[ord.temp]
        # my.df = my.III.summary$Df[ord.temp]
        # my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
        if (!VIF) {
          # my.tables.df[this.temp.var, ] = c(my.rownames[this.shift.temp],
          #                                   my.or,
          #                                   my.std.err,
          #                                   NA,
          #                                   my.dev,
          #                                   my.df,
          #                                   my.p.val)
          my.tables.df[this.temp.var,]=c(names(vars.or)[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
                                         vars.or[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
                                         vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},1],
                                         vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},2],
                                         vars.dev[this.shift.temp],
                                         vars.dev.df[this.shift.temp],
                                         vars.dev.p[this.shift.temp],rep(NA,v.p.len))
        } else{
          my.tables.df[this.temp.var, ] = c(names(vars.or)[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
                                            vars.or[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
                                            vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},1],
                                            vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},2],
                                            vars.dev[this.shift.temp],
                                            vars.dev.df[this.shift.temp],
                                            vars.dev.p[this.shift.temp],
                                            rep(NA,v.p.len))
        }
        # if (!is.null(my.factor)) {
        #   this.shift.temp = this.shift.temp + ord.temp+1
        # }
        yet.another.var=yet.another.var+1
        ord.temp = ord.temp + 1
      } else{
        if (!is.null(my.factor)) {
          this.shift.temp = this.shift.temp - ord.temp
        }else{
          this.shift.temp=this.shift.temp-total.intercepts
        }
        my.tables.df[this.temp.var, ] = c(names(vars.or)[this.shift.temp],
                                          vars.or[{this.shift.temp}],
                                          if(!is.null(dim(vars.or.confint))){
                                            vars.or.confint[{this.shift.temp},1]
                                          }else{
                                            vars.or.confint[1]
                                          },
                                          if(!is.null(dim(vars.or.confint))){
                                            vars.or.confint[{this.shift.temp},2]
                                          }else{
                                            vars.or.confint[2]
                                          },
                                          vars.dev[this.shift.temp-ord.temp],
                                          vars.dev.df[this.shift.temp-ord.temp],
                                          vars.dev.p[this.shift.temp-ord.temp],
                                          rep(NA,v.p.len))
        if (!is.null(my.factor)) {
          this.shift.temp = this.shift.temp + ord.temp
        }else{
          this.shift.temp=this.shift.temp+total.intercepts
        }
        #ord.temp = ord.temp + 1
      }
      this.shift.temp = this.shift.temp + 1
      this.temp.var = this.temp.var + 1

    }
  }

  if (type == "lm") {
    my.tables.df[this.temp.var, ] = c(
      "Residuals",
      NA,
      NA,
      my.III.summary$`Sum Sq`[this.shift.temp],
      my.III.summary$Df[this.shift.temp],
      NA,
      NA,
      rep(NA, v.p.rep)
    )
    my.tables.df[this.temp.var + 1, ] = c("Total Change",
                                          NA,
                                          NA,
                                          my.total,
                                          my.df.total,
                                          rep(NA, v.p.rep + 2))
    my.tables.df[this.temp.var+2,]=c("Total SS",NA,NA,my.total+my.III.summary$`Sum Sq`[1],my.df.total+1,rep(NA,v.p.rep+2))
    if (!VIF & !part.eta) {

    } else if (!VIF) {
      my.tables.df$p.eta = as.numeric(my.tables.df$p.eta)
    } else if (!part.eta) {
      my.tables.df$VIF = as.numeric(my.tables.df$VIF)
    } else{
      my.tables.df$p.eta = as.numeric(my.tables.df$p.eta)
      my.tables.df$VIF = as.numeric(my.tables.df$VIF)
    }
    my.tables.df$f.val = as.numeric(my.tables.df$f.val)
    my.tables.df$est = as.numeric(my.tables.df$est)
    my.tables.df$sumsq = as.numeric(my.tables.df$sumsq)
    my.tables.df$df = as.numeric(my.tables.df$df)
    my.tables.df$std.err = as.numeric(my.tables.df$std.err)
    my.tables.df$p.val = as.numeric(my.tables.df$p.val)
  } else{

    my.tables.df[this.temp.var,]=c("Total Change",NA,NA,NA,total.dev.change,total.dev.change.df,pchisq(total.dev.change,total.dev.change.df,lower.tail = F),rep(NA,v.p.rep))
    my.tables.df[this.temp.var+1,]=c("Residuals",NA,NA,NA,resid.dev,resid.df,NA,rep(NA,v.p.rep))
    my.tables.df[this.temp.var+2,]=c("Total",NA,NA,NA,total.dev,total.df,NA,rep(NA,v.p.rep))


    #my.tables.df[this.temp.var,] = c("Change from Null", NA, NA, NA, ddeviance2, ddf2, fit2)
    my.tables.df$p.odd = as.numeric(my.tables.df$p.odd)
    my.tables.df$p.odd.2.5 = as.numeric(my.tables.df$p.odd.2.5)
    my.tables.df$p.odd.97.5 = as.numeric(my.tables.df$p.odd.97.5)
    my.tables.df$deviance=as.numeric(my.tables.df$deviance)
    my.tables.df$p.val=as.numeric(my.tables.df$p.val)
    if (VIF) {
      my.tables.df$VIF = as.numeric(my.tables.df$VIF)
    }
  }



  #### Make custom glance stats ####
  if(do.glance){
    #### Can eventually make it options
    if (type == "lm") {
      glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
      glance_stats[1,]=c(paste("Method: ","QR Decomposition",if(show.contrasts){paste("<br />Adjustment Method: ",adjustment,sep="")},sep=""),rep(NA,6),rep(NA,v.p.rep))

      # glance_stats = broom::glance(my.model)
      # glance_stats = tidyr::gather(glance_stats)
      #
      #
      # glance_stats[3:{
      #   5 + v.p.rep
      # }] = NA
      # glance_stats[6 + v.p.rep] = c(glance_stats$key[7:8],
      #                               NA,
      #                               glance_stats$key[9:11],
      #                               NA,
      #                               NA,
      #                               NA,
      #                               NA,
      #                               NA)
      # glance_stats[7 + v.p.rep] = c(glance_stats$value[7:8],
      #                               NA,
      #                               glance_stats$value[9:11],
      #                               NA,
      #                               NA,
      #                               NA,
      #                               NA,
      #                               NA)
      # glance_stats = glance_stats[-7:-11, ]
      # glance_stats = glance_stats[-3, ]
      # glance_stats = glance_stats[-5, ]

    } else if (type == "glm") {
      glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
      glance_stats[1,]=c(paste("Family: ",new.model$family$family," <br /> Link: ",new.model$family$link,if(show.contrasts){paste(" <br />Adjustment: ",adjustment,sep="")},sep=""),rep(NA,6),rep(NA,v.p.rep))

      # glance_stats = broom::glance(my.model)
      # glance_stats = tidyr::gather(glance_stats)
      #
      #
      # glance_stats[[3]] = c("Family: ", "Link: ", rep(NA, 5))
      # glance_stats[[4]] = c(my.model$family$family,
      #                       my.model$family$link,
      #                       rep(NA, 5))
      # glance_stats[5:{
      #   5 + v.p.rep
      # }] = NA
      # glance_stats[6 + v.p.rep] = c(glance_stats$key[3:6], NA, NA, NA)
      # glance_stats[7 + v.p.rep] = c(glance_stats$value[3:6], NA, NA, NA)
      # glance_stats[[1]] = c("Null.dev", "Chi-Sq", "dF", "Pr(>Chisq)", NA, NA, NA)
      # glance_stats[[2]] = c(
      #   glance_stats$value[1],
      #   ddeviance2,
      #   ddf2,
      #   pvalString(fit2, digits = 3, format = "default"),
      #   NA,
      #   NA,
      #   NA
      # )
      # glance_stats = glance_stats[-5:-7, ]
      # glance_stats2=as.data.frame(matrix(ncol=7,nrow=1))
      # glance_stats2[1,]=c(glance_stats$key[1],glance_stats$value[1],NA,NA,NA,glance_stats$key[6],glance_stats$value[6])
      # glance_stats2[2,]=c(glance_stats$key[2],glance_stats$value[2],NA,NA,NA,glance_stats$key[7],glance_stats$value[7])
      # glance_stats2[3,]=c(glance_stats$key[3],glance_stats$value[3],NA,NA,NA,glance_stats$key[4],glance_stats$value[4])
      # glance_stats=glance_stats2
    } else if(type=="ord"){
      glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
      glance_stats[1,]=c(paste("Family: Ordinal <br /> Link: ",new.model$info$link,sep=""),rep(NA,6),rep(NA,v.p.rep))
      # glance_stats[[1]]=c("Null.dev","Chi-Sq","dF","Pr(>Chisq)")
      # glance_stats[[2]]=c({{-2*my.model$logLik}-ddeviance2},ddeviance2,ddf2,pvalString(fit2, digits = 3, format = "default"))
      # glance_stats[[3]]=c("Family: ","Link: ",NA,NA)
      # glance_stats[[4]]=c("ordinal",levels(my.model$info$link),NA,NA)
      # glance_stats[[5:{5+v.p.rep}]]=c(NA,NA,NA,NA)
      # glance_stats[[{6+v.p.rep}]]=c("logLik","AIC","BIC","deviance")
      # glance_stats[[{7+v.p.rep}]]=c(my.model$logLik,as.numeric(levels(my.model$info$AIC)),BIC(my.model),{-2*my.model$logLik})
    }else{

    }
  }
  #### For total
  the.length = the.length + 1
  this.temp.var = this.temp.var + 1
  #### Make table ####

  options(pixie_interactive = pix.int,
          pixie_na_string = "")


  if (type == "lm") {
    the.length=the.length+2
    my.dust = pixiedust::dust(my.tables.df) %>%
      sprinkle(cols = "p.val", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle_na_string() %>%
      sprinkle(
        rows = 1:the.length,
        cols = v.p.len,
        border = "right",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1:the.length,
        cols = 1,
        border = "left",
        border_color = "black"
      ) %>%
      sprinkle(
        rows=2,
        cols=1:v.p.len,
        border=c("top","bottom")
      )%>%
      sprinkle(rows=this.temp.var-1,
               border="top")%>%
      sprinkle(
        rows = 1,
        cols = 1:v.p.len,
        border = c("top", "bottom"),
        border_color = "black",
        part = "head"
      ) %>%
      sprinkle(rows=1,cols=1,border="left",part="head")%>%
      sprinkle(rows=1,cols=v.p.len,border="right",part="head")%>%
      sprinkle(
        rows = this.temp.var+1,
        cols = 1:v.p.len,
        border = "bottom",
        border_color = "black"
      )

    if (!VIF & !part.eta) {
      my.dust = my.dust %>%
        sprinkle(cols = c("sumsq", "est", "std.err", "f.val"),
                 round = 2) %>%
        sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
        sprinkle_colnames(
          "Variable",
          "Estimate",
          "Std. Error",
          paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
          "df",
          "F-value",
          "Pr(>F)"
        ) %>%
        sprinkle_align(rows = 1,
                       halign = "center",
                       part = "head") %>%
        sprinkle_pad(rows = 1,
                     pad = 5,
                     part = "head")

    } else if (!VIF) {
      my.dust = my.dust %>%
        sprinkle(
          cols = c("sumsq", "est", "std.err", "f.val", "p.eta"),
          round = 2
        ) %>%
        sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
        sprinkle_colnames(
          "Variable",
          "Estimate",
          "Std. Error",
          paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
          "df",
          "F-value",
          "Pr(>F)",
          "Part <br /> eta"
        ) %>%
        sprinkle_align(rows = 1,
                       halign = "center",
                       part = "head") %>%
        sprinkle_pad(rows = 1,
                     pad = 5,
                     part = "head")

    } else if (!part.eta) {
      my.dust = my.dust %>%
        sprinkle(
          cols = c("sumsq", "est", "std.err", "f.val", "VIF"),
          round = 2
        ) %>%
        sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
        sprinkle_colnames(
          "Variable",
          "Estimate",
          "Std. Error",
          paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
          "df",
          "F-value",
          "Pr(>F)",
          "VIF"
        ) %>%
        sprinkle_align(rows = 1,
                       halign = "center",
                       part = "head") %>%
        sprinkle_pad(rows = 1,
                     pad = 5,
                     part = "head")

    } else{
      my.dust = my.dust %>%
        sprinkle(
          cols = c("sumsq", "est", "std.err", "f.val", "p.eta", "VIF"),
          round = 2
        ) %>%
        sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
        sprinkle_colnames(
          "Variable",
          "Estimate",
          "Std. Error",
          paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
          "df",
          "F-value",
          "Pr(>F)",
          "Part <br /> eta",
          "VIF"
        ) %>%
        sprinkle_align(rows = 1,
                       halign = "center",
                       part = "head") %>%
        sprinkle_pad(rows = 1,
                     pad = 5,
                     part = "head")

    }

  } else if (type == "glm2") {
    my.dust = pixiedust::dust(my.tables.df) %>%
      sprinkle(cols = "p.val", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle_na_string() %>%
      sprinkle(cols = 2:5, round = 2) %>%
      sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
      sprinkle(cols = 2:{7+v.p.rep},
               pad = 5,
               part = "head") %>%
      sprinkle(
        rows = {
          this.temp.var - 2
        },
        cols = 1:v.p.len,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1:the.length,
        cols = 1:v.p.len,
        border = "right",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1:the.length,
        cols = 1,
        border = "left",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1,
        cols = 1:v.p.len,
        border = c("top", "bottom", "left", "right"),
        border_color = "black",
        part = "head"
      ) %>%
      sprinkle(
        rows = this.temp.var - 1,
        cols = 1:v.p.len,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle_align(rows = 1,
                     halign = "center",
                     part = "head")

    if (!VIF) {
      my.dust = my.dust %>% sprinkle_colnames("Variable",
                                              "Odds Ratio",
                                              "Std. Error",
                                              "z-Value",
                                              "Deviance",
                                              "df",
                                              "p-Value")
    } else{
      my.dust = my.dust %>% sprinkle_colnames(
        "Variable",
        "Odds Ratio",
        "Std. Error",
        "z-Value",
        "Deviance",
        "df",
        "p-Value",
        "VIF"
      ) %>%
        sprinkle_round(cols = 8, round = 2)
    }

  } else{
    my.dust = pixiedust::dust(my.tables.df) %>%
      sprinkle(cols = "p.val", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle_na_string() %>%
      sprinkle(cols = 2:5, round = 2) %>%
      sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
      sprinkle(cols = 2:7,
               pad = 5,
               part = "head") %>%
      sprinkle(
        rows = this.temp.var+1,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1,
        cols = 1:2,
        border = c("top", "bottom"),
        border_color = "black",
        part = "head"
      ) %>%
      sprinkle(
        rows = 1,
        cols = 3:4,
        border = c("top", "bottom"),
        border_color = "black",
        part = "head"
      ) %>%
      sprinkle(
        rows = 1,
        cols = 5:{7+v.p.rep},
        border = c("top", "bottom"),
        border_color = "black",
        part = "head"
      ) %>%
      sprinkle(
        rows = this.temp.var - 1,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = this.temp.var - 2,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 1+total.intercepts,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle(
        rows = 2+total.intercepts,
        cols = 1:7,
        border = "bottom",
        border_color = "black"
      ) %>%
      sprinkle_colnames("Variable",
                        "Odds Ratio",
                        "Conf. <br /> 2.5%","Int. <br /> 97.5%",
                        "Deviance",
                        "df",
                        "Pr(>Chi)")%>%
      sprinkle_border(cols=1,border=c("left","right"))%>%
      sprinkle_border(cols={7+v.p.rep},border="right")%>%
      sprinkle_border(cols=1,border=c("left","right"),part="head")%>%
      sprinkle_border(cols={7+v.p.rep},border="right",part="head")
  }
  if(do.glance){
    if (type == "lm") {
      my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
        sprinkle_na_string(part = "foot") %>%
        sprinkle(rows=1,merge=T,halign="center",part="foot")
      # my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
      #   sprinkle(cols = c(2, {
      #     7 + v.p.rep
      #   }),
      #   round = 3,
      #   part = "foot") %>%
      #   sprinkle(cols = 3:{
      #     5 + v.p.rep
      #   },
      #   replace = c(rep("", {
      #     12 + 4 * v.p.rep
      #   })),
      #   part = "foot") %>%
      #   sprinkle(
      #     cols = 1,
      #     replace = c("R-Square", "Adj R-Sq", "F-Statistic", "P-Value"),
      #     part = "foot"
      #   ) %>%
      #   sprinkle(
      #     cols = 2,
      #     rows = 4,
      #     fn = quote(pvalString(
      #       value, digits = 3, format = "default"
      #     )),
      #     part = "foot"
      #   ) %>%
      #   sprinkle(
      #     cols = 1:v.p.len,
      #     rows = 1,
      #     halign = "center",
      #     part = "head"
      #   ) %>%
      #   sprinkle_width(cols = 1,
      #                  width = 90,
      #                  width_units = "pt") %>%
      #   sprinkle_width(cols = 2,
      #                  width = 108,
      #                  width_units = "pt") %>%
      #   sprinkle_width(cols = 4,
      #                  width = 62,
      #                  width_units = "pt") %>%
      #   sprinkle_width(cols = 5,
      #                  width = 68,
      #                  width_units = "pt") %>%
      #   sprinkle_width(cols = 6,
      #                  width = 68,
      #                  width_units = "pt") %>%
      #   sprinkle_width(cols = 7,
      #                  width = 71,
      #                  width_units = "pt") %>%
      #   sprinkle(cols = 2,
      #            halign = "left",
      #            part = "foot")
    } else{
      my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
        sprinkle_na_string(part = "foot") %>%
        sprinkle(rows=1,merge=T,halign="center",part="foot")
    }
  }
  if (pix.int) {
    return(my.dust)
  } else{
    my.dust.print = print(my.dust, quote = F)[1]
    return(my.dust.print)
  }

}

# quick.reg.table.lm = function(my.model,
#                                    part.eta = F,
#                                    VIF = F,
#                                    myDF = my.found.df,
#                                    marginality=T,
#                                    abbrev.length = ab.len,
#                                    pix.int = T,
#                                    pix.method = "html",
#                                    type = my.reg.type,
#                                    my.factor = NULL,
#                                    do.glance=T,
#                                    show.footer=T,
#                                    adjustment = "bonferroni",
#                                    show.contrasts=F,
#                                    show.intercepts=F,
#                                    do.return=T) {
#
#   SS.type = 2
#
#     #### ANOVA TABLES ####
#     my.summary = summary(my.model)
#     my.coefficients = my.summary$coefficients
#     my.coefficients = as.data.frame(my.coefficients)
#     if(marginality){
#       my.III.summary = car::Anova(my.model, type = 2)
#     }else{
#       my.III.summary = car::Anova(my.model, type = 3)
#     }
#
#     the.length = dim(my.III.summary)[1]
#
#     #### Calculate total SS ####
#     my.total = sum(my.III.summary$`Sum Sq`[2:length(my.III.summary$`Sum Sq`)])
#     my.df.total = sum(my.III.summary$Df[2:length(my.III.summary$Df)])
#     total.intercepts = 1
#     my.rownames = c(abbreviate(rownames(my.summary$coefficients), minlength = abbrev.length),
#                     "Residuals",
#                     "Total")
#
#     treat.SS=sum(my.III.summary$`Sum Sq`[2:{length(my.III.summary$`Sum Sq`)-1}])
#     treat.df=sum(my.III.summary$Df[2:{length(my.III.summary$Df)-1}])
#
#
#
#   if (!VIF & !part.eta) {
#     v.p.len = 7
#     v.p.rep = 0
#   } else if (!VIF) {
#     v.p.len = 8
#     v.p.rep = 1
#   } else if (!part.eta) {
#     v.p.len = 8
#     v.p.rep = 1
#   } else{
#     v.p.len = 9
#     v.p.rep = 2
#   }
#
#   #### Make table if not MANOVA ####
#   my.tables.df = as.data.frame(matrix(ncol = v.p.len, nrow = 1))
#
#     if (!VIF & !part.eta) {
#       my.std.error = c(my.coefficients$`Std. Error`, NA)
#       my.estimate = c(my.coefficients$Estimate, NA)
#       names(my.tables.df) = c("rownames",
#                               "sumsq",
#                               "df",
#                               "est",
#                               "std.err",
#                               "f.val",
#                               "p.val")
#     } else if (!part.eta) {
#       my.VIF = car::vif(my.model)
#       names(my.tables.df) = c("rownames",
#                               "sumsq",
#                               "df",
#                               "est",
#                               "std.err",
#                               "f.val",
#                               "p.val",
#                               "VIF")
#     } else if (!VIF) {
#       names(my.tables.df) = c("rownames",
#                               "sumsq",
#                               "df",
#                               "est",
#                               "std.err",
#                               "f.val",
#                               "p.val",
#                               "p.eta")
#     } else{
#       my.VIF = car::vif(my.model)
#       names(my.tables.df) = c("rownames",
#                               "sumsq",
#                               "df",
#                               "est",
#                               "std.err",
#                               "f.val",
#                               "p.val",
#                               "p.eta",
#                               "VIF")
#     }
#
#   #### Make the double table entries ####
#
#   #### Was very annoying...took my frustration out on
#   #### Variable names
#
#   factor.stupid = NULL
#   factor.rownames = NULL
#   num.of.levels = NULL
#   ord.temp = 0
#   if (!is.null(my.factor)) {
#     for (i in 1:length(my.factor)) {
#       factor.stupid = c(factor.stupid, grep(paste("^", my.factor[i], "$", sep =
#                                                     ""), names(myDF)))
#
#         factor.rownames = c(factor.rownames, grep(
#           paste("^", my.factor[i], "$", sep = ""),
#           rownames(my.III.summary)
#         ))
#         num.of.levels = c(num.of.levels, length(levels(myDF[[factor.stupid[i]]])))
#         # }else if(type=="glm"){
#         #   factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+1})
#         #   num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
#
#     }
#
#   } else{
#     factor.rownames = 0L
#   }
#
#
#
#   #### Make phia stuff ####
#   if(show.contrasts){
#       if(!is.null(my.factor)){
#       my.phia.reg=quick.contrast(my.model,skip.me=T,adjustment = adjustment,SS.type = SS.type,abbrev.length = abbrev.length)
#       my.j=0
#       my.big.phia=NULL
#       phia.shift=0
#       real.shift=0
#       for(i in 1:dim(my.phia.reg)[1]){
#         if(!is.na(my.phia.reg[i,1])){
#           my.j=my.j+1
#           real.shift=real.shift+phia.shift
#           my.big.phia[[my.j]]=my.phia.reg[i,2:7]
#         }else{
#           my.big.phia[[my.j]][i-real.shift,]=my.phia.reg[i,2:7]
#         }
#         phia.shift=phia.shift+1
#       }
#
#       my.phia.rownames=NULL
#       my.phia.SS=NULL
#       my.phia.value=NULL
#       my.phia.F=NULL
#       my.phia.p=NULL
#       my.phia.err=NULL
#       for(i in 1:{my.j}){
#         my.phia.rownames[[i]]=my.big.phia[[i]]$names
#         my.phia.SS[[i]]=my.big.phia[[i]]$`Sum of Sq`
#         my.phia.value[[i]]=my.big.phia[[i]]$Value
#         my.phia.F[[i]]=my.big.phia[[i]][,5]
#         my.phia.p[[i]]=my.big.phia[[i]][,6]
#
#       }
#     }
#   }
#
#
#
#   my.factor.var = 1
#   this.temp.var = 1
#   this.shift.temp = 1
#   yet.another.var = 1
#   my.shift = 0
#   other.temp = 2
#   ord.temp = 0
#   phia.temp=1
#
#     dang.length = length(rownames(my.III.summary))
#
#
#
#
#   while (this.shift.temp < dang.length) {
#     if (is.na(factor.rownames[my.factor.var])) {
#       my.factor.rownames = 1
#
#     } else{
#       my.factor.rownames = factor.rownames[my.factor.var]
#
#     }
#     #### LOOP ####
#     if (this.shift.temp == 1) {
#       i = 1
#       if(type=="ord" | type=="glm"){
#         my.tables.df[this.temp.var, ] = c(
#           "Intercept Change",
#           NA,
#           NA,
#           NA,
#           my.int.dev.total,
#           my.int.dev.df,
#           pchisq(my.int.dev.total,my.int.dev.df,lower.tail = F),
#           rep(NA, v.p.rep))
#
#         this.temp.var=this.temp.var+1
#       }
#       while (i <= total.intercepts) {
#         if (type == "lm") {
#           my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
#           my.df = my.III.summary$Df[this.shift.temp]
#           my.est = my.estimate[this.shift.temp]
#           my.std.err = my.std.error[this.shift.temp]
#           my.f.val = my.III.summary$`F value`[this.shift.temp]
#           my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
#           my.tables.df[this.temp.var, ] = c(
#             my.rownames[this.shift.temp],
#             summary(my.model)[[4]][this.shift.temp,1],
#             summary(my.model)[[4]][this.shift.temp,2],
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             rep(NA, v.p.rep)
#           )
#         } else if(type=="glm"){
#           #### HERE IS WHERE I LEFT OFF!! ####
#           # my.or = exp(my.estimate[this.shift.temp])
#           # my.est = my.estimate[this.shift.temp]
#           # my.z.val = my.summary$coefficients[this.shift.temp, 3]
#           # my.std.err = my.std.error[this.shift.temp]
#           #my.dev=my.III.summary$`LR Chisq`[this.shift.temp]
#           #my.df=my.III.summary$Df[this.shift.temp]
#           # my.p.val = my.summary$coefficients[{
#           #   this.shift.temp
#           #}, 4]
#
#           my.tables.df[this.temp.var, ] = c(
#             paste(names(attr(my.model$model[[i]],"labels"))[1],"-",names(attr(my.model$model[[i]],"labels"))[2],sep=""),
#             my.int.dev.or[i],
#             my.int.dev.or.confint[1],
#             my.int.dev.or.confint[2],
#             my.int.dev[i],
#             1,
#             NA,
#             rep(NA, v.p.rep))
#
#           #this.temp.var=this.temp.var+1
#           #my.tables.df[this.temp.var,]=c("Treatment",NA,NA,NA,total.dev,total.df,dchisq(total.dev,total.df),rep(NA,v.p.rep))
#
#         }else{
#           my.tables.df[this.temp.var, ] = c(
#             my.int.names[i],
#             my.int.dev.or[i],
#             my.int.dev.or.confint[i,1],
#             my.int.dev.or.confint[i,2],
#             my.int.dev[i],
#             NA,
#             NA,
#             rep(NA, v.p.rep))
#
#           #this.temp.var=this.temp.var+1
#
#         }
#
#         this.shift.temp = this.shift.temp + 1
#
#         this.temp.var = this.temp.var + 1
#         i = i + 1
#
#       }
#       if(type=="ord" | type=="glm"){
#         my.tables.df[this.temp.var,]=c("Treatment Change",NA,NA,NA,treat.dev,vars.df,pchisq(treat.dev,vars.df,lower.tail = F))
#         this.temp.var=this.temp.var+1
#       }else{
#         my.tables.df[this.temp.var,]=c("Treatment",NA,NA,treat.SS,treat.df,glance(my.model)[4],glance(my.model)[5],rep(NA,v.p.rep))
#         this.temp.var=this.temp.var+1
#       }
#     } else if (this.shift.temp == my.factor.rownames) {
#       if (type == "lm") {
#         my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
#         my.df = my.III.summary$Df[this.shift.temp]
#         my.est = my.estimate[this.shift.temp]
#         my.std.err = my.std.error[this.shift.temp]
#         my.z.val = NA
#         my.f.val = my.III.summary$`F value`[this.shift.temp]
#         my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
#         if(!VIF & !part.eta){
#           my.tables.df[this.temp.var, ] = c(
#             my.factor[yet.another.var],
#             NA,
#             NA,
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             rep(NA, v.p.rep)
#           )
#         }else if(part.eta & !VIF){
#           my.p.eta = my.sumsq / my.total
#           my.tables.df[this.temp.var, ] = c(
#             my.factor[yet.another.var],
#             NA,
#             NA,
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.p.eta
#           )
#         }else if(!part.eta & VIF){
#
#           my.tables.df[this.temp.var, ] = c(
#             my.factor[yet.another.var],
#             NA,
#             NA,
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.VIF[this.shift.temp - 1, 1]
#           )
#         }else{
#           my.p.eta = my.sumsq / my.total
#           my.tables.df[this.temp.var, ] = c(
#             my.factor[yet.another.var],
#             NA,
#             NA,
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.p.eta,
#             my.VIF[this.shift.temp - 1, 1]
#           )
#         }
#       }else{
#
#         # my.or = NA
#         # my.est = NA
#         # my.std.err = NA
#         # my.dev = my.III.summary$`LR Chisq`[ord.temp]
#         # my.df = my.III.summary$Df[ord.temp]
#         # my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
#         my.tables.df[this.temp.var,]=c(my.factor[this.shift.temp-1],
#                                        NA,
#                                        NA,
#                                        NA,
#                                        drop1(new.model,test="Chi")$`LRT`[this.shift.temp],
#                                        drop1(new.model,test="Chi")$`Df`[this.shift.temp],
#                                        drop1(new.model,test="Chi")$`Pr(>Chi)`[this.shift.temp],
#                                        rep(NA,v.p.rep))
#         #
#         # my.tables.df[this.temp.var, ] = c(my.factor[yet.another.var],
#         #                                   my.or,
#         #                                   my.est,
#         #                                   my.std.err,
#         #                                   my.dev,
#         #                                   my.df,
#         #                                   my.p.val)
#         #ord.temp = ord.temp + 1
#       }
#       # else{
#       #   my.or = NA
#       #   my.est = NA
#       #   my.std.err = NA
#       #   my.dev = my.III.summary$Chisq[ord.temp]
#       #   my.df = my.III.summary$Df[ord.temp]
#       #   my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
#       #   my.tables.df[this.temp.var, ] = c(my.factor[yet.another.var],
#       #                                     my.or,
#       #                                     my.est,
#       #                                     my.std.err,
#       #                                     my.dev,
#       #                                     my.df,
#       #                                     my.p.val)
#       #   ord.temp = ord.temp + 1
#       #
#       # }
#       yet.another.var = yet.another.var + 1
#       this.temp.var = this.temp.var + 1
#       this.shift.temp = this.shift.temp + 1
#       if(type=="lm"){
#         other.other.temp = 2
#       }else{
#         other.other.temp=1
#       }
#
#       #### INTERACTION EFFECTS? NOT WORRIED YET ####
#       if ({length(grepl(":", my.summary$coefficients[other.temp, 1])) >
#           0}) {
#
#         #### I think the while should not be +1
#         if(type=="glm" | type=="ord"){
#           num.of.levels[my.factor.var]=num.of.levels[my.factor.var]-1
#         }
#         while (other.other.temp < {
#           num.of.levels[my.factor.var] + 1
#         }) {
#           if (type == "lm") {
#             if(show.contrasts){
#               my.sumsq = NA
#               my.df = NA
#               my.est = my.estimate[other.temp]
#               my.std.err = my.std.error[other.temp]
#               #### NEED TO FIX ####
#               my.f.val = {
#                 my.summary$coefficients[other.temp, 3] ^ 2
#               }
#               my.p.val = my.summary$coefficients[other.temp, 4]
#               my.tables.df[this.temp.var, ] = c(
#                 my.phia.rownames[[phia.temp]][other.temp-1],
#                 my.phia.value[[phia.temp]][other.temp-1],
#                 my.est,
#                 my.phia.SS[[phia.temp]][other.temp-1],
#                 1,
#                 my.phia.F[[phia.temp]][other.temp-1],
#                 my.phia.p[[phia.temp]][other.temp-1],
#                 rep(NA, v.p.rep)
#               )
#             }
#             ord.temp=ord.temp + 1
#           } else if (type == "glm") {
#             # my.or = exp(my.estimate[other.temp])
#             # my.est = my.estimate[other.temp]
#             # my.std.err = my.std.error[other.temp]
#             # my.z.val = my.summary$coefficients[other.temp, 3]
#             # my.dev = my.III.summary$`LR Chisq`[other.temp]
#             # my.df = my.III.summary$Df[other.temp]
#             # my.p.val = my.summary$coefficients[other.temp, 4]
#             # my.tables.df[this.temp.var, ] = c(my.rownames[other.temp],
#             #                                   my.or,
#             #                                   my.std.err,
#             #                                   my.z.val,
#             #                                   NA,
#             #                                   NA,
#             #                                   my.p.val)
#             if(show.contrasts){
#               my.tables.df[this.temp.var, ] = c(my.phia.rownames[other.temp-1],
#                                                 vars.or[other.temp-1],
#                                                 vars.or.confint[other.temp-1,1],
#                                                 vars.or.confint[other.temp-1,2],
#                                                 my.phia[other.temp-1,3],
#                                                 my.phia[other.temp-1,4],
#                                                 my.phia[other.temp-1,6],rep(NA,v.p.rep))
#               #this.shift.temp = this.shift.temp + 1
#             }
#             ord.temp=ord.temp + 1
#           } else{
#             my.or = exp(my.estimate[other.temp + total.intercepts - 1])
#             my.est = my.estimate[other.temp + total.intercepts - 1]
#             my.std.err = my.std.error[other.temp + total.intercepts -
#                                         1]
#             my.z.val = my.summary$coefficients[{
#               other.temp + total.intercepts - 1
#             }, 3]
#             #my.dev=my.III.summary$Chisq[other.temp]
#             #my.df=my.III.summary$Df[other.temp]
#             my.p.val = my.summary$coefficients[{
#               other.temp + total.intercepts - 1
#             }, 4]
#             if (!VIF) {
#               my.tables.df[this.temp.var, ] = c(my.rownames[other.temp +
#                                                               total.intercepts - 1],
#                                                 my.or,
#                                                 my.std.err,
#                                                 my.z.val,
#                                                 NA,
#                                                 NA,
#                                                 my.p.val)
#             } else{
#               my.tables.df[this.temp.var, ] = c(my.rownames[other.temp +
#                                                               total.intercepts - 1],
#                                                 my.or,
#                                                 my.std.err,
#                                                 my.z.val,
#                                                 NA,
#                                                 NA,
#                                                 my.p.val,
#                                                 my.VIF[this.shift.temp - 1, 1])
#
#             }
#             this.shift.temp = this.shift.temp + 1
#           }
#           this.temp.var = this.temp.var + 1
#           other.temp = other.temp + 1
#           other.other.temp = other.other.temp + 1
#           the.length = the.length + 1
#
#         }
#         phia.temp=phia.temp+1
#         if(!show.contrasts){this.temp.var=this.temp.var-ord.temp}
#       } else{
#
#       }
#
#       if (my.factor.var == 1) {
#         my.shift = {
#           my.shift + other.temp - 2
#         }
#
#       } else{
#         my.shift = my.shift + other.temp
#
#       }
#
#       my.factor.var = my.factor.var + 1
#
#     } else{
#       if (type == "lm") {
#         my.sumsq = my.III.summary$`Sum Sq`[this.shift.temp]
#         my.df = my.III.summary$Df[this.shift.temp]
#         my.est = my.estimate[this.shift.temp]
#         my.std.err = my.std.error[this.shift.temp]
#         my.f.val = my.III.summary$`F value`[this.shift.temp]
#         my.p.val = my.III.summary$`Pr(>F)`[this.shift.temp]
#         if (!VIF & !part.eta) {
#           my.tables.df[this.temp.var, ] = c(
#             rownames(my.III.summary)[this.shift.temp],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val
#           )
#         } else if (!part.eta) {
#           my.tables.df[this.temp.var, ] = c(
#             rownames(my.III.summary)[this.shift.temp],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.VIF[this.shift.temp - 1, 1]
#           )
#         } else if (!VIF) {
#           my.p.eta = my.sumsq / my.total
#           my.tables.df[this.temp.var, ] = c(
#             rownames(my.III.summary)[this.shift.temp],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.p.eta
#           )
#         } else{
#           my.p.eta = my.sumsq / my.total
#           my.tables.df[this.temp.var, ] = c(
#             rownames(my.III.summary)[this.shift.temp],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,1],
#             summary(my.model)[[4]][this.shift.temp+ord.temp-1,2],
#             my.sumsq,
#             my.df,
#             my.f.val,
#             my.p.val,
#             my.p.eta,
#             my.VIF[this.shift.temp - 1, 1]
#           )
#         }
#         ord.temp=ord.temp+1
#       } else if (type == "glm") {
#         # if (!is.null(my.factor)) {
#         #   this.shift.temp = this.shift.temp - ord.temp-1
#         # }
#         # my.or = exp(my.estimate[this.shift.temp])
#         # my.est = my.estimate[this.shift.temp]
#         # my.std.err = my.std.error[this.shift.temp]
#         # my.z.val = my.summary$coefficients[other.temp, 3]
#         # my.dev = my.III.summary$`LR Chisq`[ord.temp]
#         # my.df = my.III.summary$Df[ord.temp]
#         # my.p.val = my.III.summary$`Pr(>Chisq)`[ord.temp]
#         if (!VIF) {
#           # my.tables.df[this.temp.var, ] = c(my.rownames[this.shift.temp],
#           #                                   my.or,
#           #                                   my.std.err,
#           #                                   NA,
#           #                                   my.dev,
#           #                                   my.df,
#           #                                   my.p.val)
#           my.tables.df[this.temp.var,]=c(names(vars.or)[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
#                                          vars.or[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
#                                          vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},1],
#                                          vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},2],
#                                          vars.dev[this.shift.temp],
#                                          vars.dev.df[this.shift.temp],
#                                          vars.dev.p[this.shift.temp],rep(NA,v.p.len))
#         } else{
#           my.tables.df[this.temp.var, ] = c(names(vars.or)[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
#                                             vars.or[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1}],
#                                             vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},1],
#                                             vars.or.confint[{this.shift.temp+sum(vars.dev.df[1:{this.shift.temp-1}],na.rm = T)-yet.another.var+1},2],
#                                             vars.dev[this.shift.temp],
#                                             vars.dev.df[this.shift.temp],
#                                             vars.dev.p[this.shift.temp],
#                                             rep(NA,v.p.len))
#         }
#         # if (!is.null(my.factor)) {
#         #   this.shift.temp = this.shift.temp + ord.temp+1
#         # }
#         yet.another.var=yet.another.var+1
#         ord.temp = ord.temp + 1
#       } else{
#         if (!is.null(my.factor)) {
#           this.shift.temp = this.shift.temp - ord.temp
#         }else{
#           this.shift.temp=this.shift.temp-total.intercepts
#         }
#         my.tables.df[this.temp.var, ] = c(names(vars.or)[this.shift.temp],
#                                           vars.or[{this.shift.temp}],
#                                           if(!is.null(dim(vars.or.confint))){
#                                             vars.or.confint[{this.shift.temp},1]
#                                           }else{
#                                             vars.or.confint[1]
#                                           },
#                                           if(!is.null(dim(vars.or.confint))){
#                                             vars.or.confint[{this.shift.temp},2]
#                                           }else{
#                                             vars.or.confint[2]
#                                           },
#                                           vars.dev[this.shift.temp-ord.temp],
#                                           vars.dev.df[this.shift.temp-ord.temp],
#                                           vars.dev.p[this.shift.temp-ord.temp],
#                                           rep(NA,v.p.len))
#         if (!is.null(my.factor)) {
#           this.shift.temp = this.shift.temp + ord.temp
#         }else{
#           this.shift.temp=this.shift.temp+total.intercepts
#         }
#         #ord.temp = ord.temp + 1
#       }
#       this.shift.temp = this.shift.temp + 1
#       this.temp.var = this.temp.var + 1
#
#     }
#   }
#
#   if (type == "lm") {
#     my.tables.df[this.temp.var, ] = c(
#       "Residuals",
#       NA,
#       NA,
#       my.III.summary$`Sum Sq`[this.shift.temp],
#       my.III.summary$Df[this.shift.temp],
#       NA,
#       NA,
#       rep(NA, v.p.rep)
#     )
#     my.tables.df[this.temp.var + 1, ] = c("Total Change",
#                                           NA,
#                                           NA,
#                                           my.total,
#                                           my.df.total,
#                                           rep(NA, v.p.rep + 2))
#     my.tables.df[this.temp.var+2,]=c("Total SS",NA,NA,my.total+my.III.summary$`Sum Sq`[1],my.df.total+1,rep(NA,v.p.rep+2))
#     if (!VIF & !part.eta) {
#
#     } else if (!VIF) {
#       my.tables.df$p.eta = as.numeric(my.tables.df$p.eta)
#     } else if (!part.eta) {
#       my.tables.df$VIF = as.numeric(my.tables.df$VIF)
#     } else{
#       my.tables.df$p.eta = as.numeric(my.tables.df$p.eta)
#       my.tables.df$VIF = as.numeric(my.tables.df$VIF)
#     }
#     my.tables.df$f.val = as.numeric(my.tables.df$f.val)
#     my.tables.df$est = as.numeric(my.tables.df$est)
#     my.tables.df$sumsq = as.numeric(my.tables.df$sumsq)
#     my.tables.df$df = as.numeric(my.tables.df$df)
#     my.tables.df$std.err = as.numeric(my.tables.df$std.err)
#     my.tables.df$p.val = as.numeric(my.tables.df$p.val)
#   } else{
#
#     my.tables.df[this.temp.var,]=c("Total Change",NA,NA,NA,total.dev.change,total.dev.change.df,pchisq(total.dev.change,total.dev.change.df,lower.tail = F),rep(NA,v.p.rep))
#     my.tables.df[this.temp.var+1,]=c("Residuals",NA,NA,NA,resid.dev,resid.df,NA,rep(NA,v.p.rep))
#     my.tables.df[this.temp.var+2,]=c("Total",NA,NA,NA,total.dev,total.df,NA,rep(NA,v.p.rep))
#
#
#     #my.tables.df[this.temp.var,] = c("Change from Null", NA, NA, NA, ddeviance2, ddf2, fit2)
#     my.tables.df$p.odd = as.numeric(my.tables.df$p.odd)
#     my.tables.df$p.odd.2.5 = as.numeric(my.tables.df$p.odd.2.5)
#     my.tables.df$p.odd.97.5 = as.numeric(my.tables.df$p.odd.97.5)
#     my.tables.df$deviance=as.numeric(my.tables.df$deviance)
#     my.tables.df$p.val=as.numeric(my.tables.df$p.val)
#     if (VIF) {
#       my.tables.df$VIF = as.numeric(my.tables.df$VIF)
#     }
#   }
#
#
#
#   #### Make custom glance stats ####
#   if(do.glance){
#     #### Can eventually make it options
#     if (type == "lm") {
#       glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
#       glance_stats[1,]=c(paste("Method: ","QR Decomposition",if(show.contrasts){paste("<br />Adjustment Method: ",adjustment,sep="")},sep=""),rep(NA,6),rep(NA,v.p.rep))
#
#       # glance_stats = broom::glance(my.model)
#       # glance_stats = tidyr::gather(glance_stats)
#       #
#       #
#       # glance_stats[3:{
#       #   5 + v.p.rep
#       # }] = NA
#       # glance_stats[6 + v.p.rep] = c(glance_stats$key[7:8],
#       #                               NA,
#       #                               glance_stats$key[9:11],
#       #                               NA,
#       #                               NA,
#       #                               NA,
#       #                               NA,
#       #                               NA)
#       # glance_stats[7 + v.p.rep] = c(glance_stats$value[7:8],
#       #                               NA,
#       #                               glance_stats$value[9:11],
#       #                               NA,
#       #                               NA,
#       #                               NA,
#       #                               NA,
#       #                               NA)
#       # glance_stats = glance_stats[-7:-11, ]
#       # glance_stats = glance_stats[-3, ]
#       # glance_stats = glance_stats[-5, ]
#
#     } else if (type == "glm") {
#       glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
#       glance_stats[1,]=c(paste("Family: ",new.model$family$family," <br /> Link: ",new.model$family$link,if(show.contrasts){paste(" <br />Adjustment: ",adjustment,sep="")},sep=""),rep(NA,6),rep(NA,v.p.rep))
#
#       # glance_stats = broom::glance(my.model)
#       # glance_stats = tidyr::gather(glance_stats)
#       #
#       #
#       # glance_stats[[3]] = c("Family: ", "Link: ", rep(NA, 5))
#       # glance_stats[[4]] = c(my.model$family$family,
#       #                       my.model$family$link,
#       #                       rep(NA, 5))
#       # glance_stats[5:{
#       #   5 + v.p.rep
#       # }] = NA
#       # glance_stats[6 + v.p.rep] = c(glance_stats$key[3:6], NA, NA, NA)
#       # glance_stats[7 + v.p.rep] = c(glance_stats$value[3:6], NA, NA, NA)
#       # glance_stats[[1]] = c("Null.dev", "Chi-Sq", "dF", "Pr(>Chisq)", NA, NA, NA)
#       # glance_stats[[2]] = c(
#       #   glance_stats$value[1],
#       #   ddeviance2,
#       #   ddf2,
#       #   pvalString(fit2, digits = 3, format = "default"),
#       #   NA,
#       #   NA,
#       #   NA
#       # )
#       # glance_stats = glance_stats[-5:-7, ]
#       # glance_stats2=as.data.frame(matrix(ncol=7,nrow=1))
#       # glance_stats2[1,]=c(glance_stats$key[1],glance_stats$value[1],NA,NA,NA,glance_stats$key[6],glance_stats$value[6])
#       # glance_stats2[2,]=c(glance_stats$key[2],glance_stats$value[2],NA,NA,NA,glance_stats$key[7],glance_stats$value[7])
#       # glance_stats2[3,]=c(glance_stats$key[3],glance_stats$value[3],NA,NA,NA,glance_stats$key[4],glance_stats$value[4])
#       # glance_stats=glance_stats2
#     } else if(type=="ord"){
#       glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
#       glance_stats[1,]=c(paste("Family: Ordinal <br /> Link: ",new.model$info$link,sep=""),rep(NA,6),rep(NA,v.p.rep))
#       # glance_stats[[1]]=c("Null.dev","Chi-Sq","dF","Pr(>Chisq)")
#       # glance_stats[[2]]=c({{-2*my.model$logLik}-ddeviance2},ddeviance2,ddf2,pvalString(fit2, digits = 3, format = "default"))
#       # glance_stats[[3]]=c("Family: ","Link: ",NA,NA)
#       # glance_stats[[4]]=c("ordinal",levels(my.model$info$link),NA,NA)
#       # glance_stats[[5:{5+v.p.rep}]]=c(NA,NA,NA,NA)
#       # glance_stats[[{6+v.p.rep}]]=c("logLik","AIC","BIC","deviance")
#       # glance_stats[[{7+v.p.rep}]]=c(my.model$logLik,as.numeric(levels(my.model$info$AIC)),BIC(my.model),{-2*my.model$logLik})
#     }else{
#
#     }
#   }
#   #### For total
#   the.length = the.length + 1
#   this.temp.var = this.temp.var + 1
#   #### Make table ####
#
#   options(pixie_interactive = pix.int,
#           pixie_na_string = "")
#
#
#   if (type == "lm") {
#     the.length=the.length+2
#     my.dust = pixiedust::dust(my.tables.df) %>%
#       sprinkle(cols = "p.val", fn = quote(pvalString(
#         value, digits = 3, format = "default"
#       ))) %>%
#       sprinkle_print_method(pix.method) %>%
#       sprinkle_na_string() %>%
#       sprinkle(
#         rows = 1:the.length,
#         cols = v.p.len,
#         border = "right",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1:the.length,
#         cols = 1,
#         border = "left",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows=2,
#         cols=1:v.p.len,
#         border=c("top","bottom")
#       )%>%
#       sprinkle(rows=this.temp.var-1,
#                border="top")%>%
#       sprinkle(
#         rows = 1,
#         cols = 1:v.p.len,
#         border = c("top", "bottom"),
#         border_color = "black",
#         part = "head"
#       ) %>%
#       sprinkle(rows=1,cols=1,border="left",part="head")%>%
#       sprinkle(rows=1,cols=v.p.len,border="right",part="head")%>%
#       sprinkle(
#         rows = this.temp.var+1,
#         cols = 1:v.p.len,
#         border = "bottom",
#         border_color = "black"
#       )
#
#     if (!VIF & !part.eta) {
#       my.dust = my.dust %>%
#         sprinkle(cols = c("sumsq", "est", "std.err", "f.val"),
#                  round = 2) %>%
#         sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
#         sprinkle_colnames(
#           "Variable",
#           "Estimate",
#           "Std. Error",
#           paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
#           "df",
#           "F-value",
#           "Pr(>F)"
#         ) %>%
#         sprinkle_align(rows = 1,
#                        halign = "center",
#                        part = "head") %>%
#         sprinkle_pad(rows = 1,
#                      pad = 5,
#                      part = "head")
#
#     } else if (!VIF) {
#       my.dust = my.dust %>%
#         sprinkle(
#           cols = c("sumsq", "est", "std.err", "f.val", "p.eta"),
#           round = 2
#         ) %>%
#         sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
#         sprinkle_colnames(
#           "Variable",
#           "Estimate",
#           "Std. Error",
#           paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
#           "df",
#           "F-value",
#           "Pr(>F)",
#           "Part <br /> eta"
#         ) %>%
#         sprinkle_align(rows = 1,
#                        halign = "center",
#                        part = "head") %>%
#         sprinkle_pad(rows = 1,
#                      pad = 5,
#                      part = "head")
#
#     } else if (!part.eta) {
#       my.dust = my.dust %>%
#         sprinkle(
#           cols = c("sumsq", "est", "std.err", "f.val", "VIF"),
#           round = 2
#         ) %>%
#         sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
#         sprinkle_colnames(
#           "Variable",
#           "Estimate",
#           "Std. Error",
#           paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
#           "df",
#           "F-value",
#           "Pr(>F)",
#           "VIF"
#         ) %>%
#         sprinkle_align(rows = 1,
#                        halign = "center",
#                        part = "head") %>%
#         sprinkle_pad(rows = 1,
#                      pad = 5,
#                      part = "head")
#
#     } else{
#       my.dust = my.dust %>%
#         sprinkle(
#           cols = c("sumsq", "est", "std.err", "f.val", "p.eta", "VIF"),
#           round = 2
#         ) %>%
#         sprinkle(cols = c(2, 3, 4, 5, 6, 7, 8), pad = 5) %>%
#         sprinkle_colnames(
#           "Variable",
#           "Estimate",
#           "Std. Error",
#           paste("Type ", SS.type, "<br /> Sums of Sq", sep = ""),
#           "df",
#           "F-value",
#           "Pr(>F)",
#           "Part <br /> eta",
#           "VIF"
#         ) %>%
#         sprinkle_align(rows = 1,
#                        halign = "center",
#                        part = "head") %>%
#         sprinkle_pad(rows = 1,
#                      pad = 5,
#                      part = "head")
#
#     }
#
#   } else if (type == "glm2") {
#     my.dust = pixiedust::dust(my.tables.df) %>%
#       sprinkle(cols = "p.val", fn = quote(pvalString(
#         value, digits = 3, format = "default"
#       ))) %>%
#       sprinkle_print_method(pix.method) %>%
#       sprinkle_na_string() %>%
#       sprinkle(cols = 2:5, round = 2) %>%
#       sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
#       sprinkle(cols = 2:{7+v.p.rep},
#                pad = 5,
#                part = "head") %>%
#       sprinkle(
#         rows = {
#           this.temp.var - 2
#         },
#         cols = 1:v.p.len,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1:the.length,
#         cols = 1:v.p.len,
#         border = "right",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1:the.length,
#         cols = 1,
#         border = "left",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1,
#         cols = 1:v.p.len,
#         border = c("top", "bottom", "left", "right"),
#         border_color = "black",
#         part = "head"
#       ) %>%
#       sprinkle(
#         rows = this.temp.var - 1,
#         cols = 1:v.p.len,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle_align(rows = 1,
#                      halign = "center",
#                      part = "head")
#
#     if (!VIF) {
#       my.dust = my.dust %>% sprinkle_colnames("Variable",
#                                               "Odds Ratio",
#                                               "Std. Error",
#                                               "z-Value",
#                                               "Deviance",
#                                               "df",
#                                               "p-Value")
#     } else{
#       my.dust = my.dust %>% sprinkle_colnames(
#         "Variable",
#         "Odds Ratio",
#         "Std. Error",
#         "z-Value",
#         "Deviance",
#         "df",
#         "p-Value",
#         "VIF"
#       ) %>%
#         sprinkle_round(cols = 8, round = 2)
#     }
#
#   } else{
#     my.dust = pixiedust::dust(my.tables.df) %>%
#       sprinkle(cols = "p.val", fn = quote(pvalString(
#         value, digits = 3, format = "default"
#       ))) %>%
#       sprinkle_print_method(pix.method) %>%
#       sprinkle_na_string() %>%
#       sprinkle(cols = 2:5, round = 2) %>%
#       sprinkle(cols = c(2, 3, 4, 5, 6, 7), pad = 5) %>%
#       sprinkle(cols = 2:7,
#                pad = 5,
#                part = "head") %>%
#       sprinkle(
#         rows = this.temp.var+1,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1,
#         cols = 1:2,
#         border = c("top", "bottom"),
#         border_color = "black",
#         part = "head"
#       ) %>%
#       sprinkle(
#         rows = 1,
#         cols = 3:4,
#         border = c("top", "bottom"),
#         border_color = "black",
#         part = "head"
#       ) %>%
#       sprinkle(
#         rows = 1,
#         cols = 5:{7+v.p.rep},
#         border = c("top", "bottom"),
#         border_color = "black",
#         part = "head"
#       ) %>%
#       sprinkle(
#         rows = this.temp.var - 1,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = this.temp.var - 2,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 1+total.intercepts,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle(
#         rows = 2+total.intercepts,
#         cols = 1:7,
#         border = "bottom",
#         border_color = "black"
#       ) %>%
#       sprinkle_colnames("Variable",
#                         "Odds Ratio",
#                         "Conf. <br /> 2.5%","Int. <br /> 97.5%",
#                         "Deviance",
#                         "df",
#                         "Pr(>Chi)")%>%
#       sprinkle_border(cols=1,border=c("left","right"))%>%
#       sprinkle_border(cols={7+v.p.rep},border="right")%>%
#       sprinkle_border(cols=1,border=c("left","right"),part="head")%>%
#       sprinkle_border(cols={7+v.p.rep},border="right",part="head")
#   }
#   if(do.glance){
#     if (type == "lm") {
#       my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
#         sprinkle_na_string(part = "foot") %>%
#         sprinkle(rows=1,merge=T,halign="center",part="foot")
#       # my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
#       #   sprinkle(cols = c(2, {
#       #     7 + v.p.rep
#       #   }),
#       #   round = 3,
#       #   part = "foot") %>%
#       #   sprinkle(cols = 3:{
#       #     5 + v.p.rep
#       #   },
#       #   replace = c(rep("", {
#       #     12 + 4 * v.p.rep
#       #   })),
#       #   part = "foot") %>%
#       #   sprinkle(
#       #     cols = 1,
#       #     replace = c("R-Square", "Adj R-Sq", "F-Statistic", "P-Value"),
#       #     part = "foot"
#       #   ) %>%
#       #   sprinkle(
#       #     cols = 2,
#       #     rows = 4,
#       #     fn = quote(pvalString(
#       #       value, digits = 3, format = "default"
#       #     )),
#       #     part = "foot"
#       #   ) %>%
#       #   sprinkle(
#       #     cols = 1:v.p.len,
#       #     rows = 1,
#       #     halign = "center",
#       #     part = "head"
#       #   ) %>%
#       #   sprinkle_width(cols = 1,
#       #                  width = 90,
#       #                  width_units = "pt") %>%
#       #   sprinkle_width(cols = 2,
#       #                  width = 108,
#       #                  width_units = "pt") %>%
#       #   sprinkle_width(cols = 4,
#       #                  width = 62,
#       #                  width_units = "pt") %>%
#       #   sprinkle_width(cols = 5,
#       #                  width = 68,
#       #                  width_units = "pt") %>%
#       #   sprinkle_width(cols = 6,
#       #                  width = 68,
#       #                  width_units = "pt") %>%
#       #   sprinkle_width(cols = 7,
#       #                  width = 71,
#       #                  width_units = "pt") %>%
#       #   sprinkle(cols = 2,
#       #            halign = "left",
#       #            part = "foot")
#     } else{
#       my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
#         sprinkle_na_string(part = "foot") %>%
#         sprinkle(rows=1,merge=T,halign="center",part="foot")
#     }
#   }
#   if (pix.int) {
#     return(my.dust)
#   } else{
#     my.dust.print = print(my.dust, quote = F)[1]
#     return(my.dust.print)
#   }
#
# }

# quick.reg.table.default=function(my.model, myDF=my.found.df, my.factor=NULL, SS.type=3, pix.int=T,pix.method="html",type=my.reg.type,test.stat="Wilks"){
#   library(pixiedust)
#   library(broom)
#
#   #### Find type
#   # my.call=as.character(my.model$call)
#   # my.split.call=strsplit(my.call,"\\\\(")
#   # my.reg.type2=my.split.call[[1]][1]
#   # if(my.reg.type2=="lm" | my.reg.type2 == "stats::lm"){
#   #   my.reg.type="lm"
#   # }else if(my.reg.type2=="glm" | my.reg.type2== "stats::glm"){
#   #   my.reg.type="glm"
#   # }else if(my.reg.type2=="manova" | my.reg.type2=="stats::manova"){
#   #   my.reg.type="manova"
#   # }else if(my.reg.type2=="clm" | my.reg.type2=="ordinal::clm"){
#   #   my.reg.type="ord"
#   # }else{
#   #   stop("Type not supported")
#   # }
#   my.found.df=my.model$call$data
#   if(is.null(my.found.df)){
#     stop(paste("No data frame found"))
#   }
#
#   #### Make factor list
#   if(type=="manova" | type=="stats::manova"){
#     x3=capture.output(car::Anova(my.model,type=SS.type,test=test.stat))
#     my.manova.test=data.frame(matrix(ncol=7,nrow=1))
#     my.var.temp=4
#
#     while(my.var.temp<{length(x3)-1}){
#
#       test=strsplit(x3[my.var.temp],"\\\\s+")
#
#       if(length(test[[1]])==9){
#
#         test2=test[[1]][-9]
#         test2=test2[-7]
#
#       }else if(length(test[[1]])==8){
#
#         test2=test[[1]][-8]
#
#       }else{
#
#         test2=test[[1]]
#       }
#
#       my.manova.test[{my.var.temp-3},]=test2
#       my.var.temp=my.var.temp+1
#
#     }
#
#     my.manova.test[[2]]=as.numeric(my.manova.test[[2]])
#     my.manova.test[[3]]=as.numeric(my.manova.test[[3]])
#     my.manova.test[[4]]=as.numeric(my.manova.test[[4]])
#     my.manova.test[[5]]=as.numeric(my.manova.test[[5]])
#     my.manova.test[[6]]=as.numeric(my.manova.test[[6]])
#     my.manova.test[[7]]=as.numeric(my.manova.test[[7]])
#
#     options(pixie_interactive = pix.int)
#     my.manova.pixie=pixiedust::dust(my.manova.test)%>%
#       sprinkle_print_method(pix.method)%>%
#       sprinkle(cols="X7",fn=quote(pvalString(value,digits=3,format="default")))%>%
#       sprinkle(cols="X3",round=3)%>%
#       sprinkle_colnames("","df",paste(test.stat," <br /> Statistic"),"approx <br /> F-value","num df","den df","Pr(>F)")%>%
#       sprinkle(cols=1:7,rows={length(x3)-5},border=c("bottom","left","right"))%>%
#       sprinkle(cols=1:7,pad=10)%>%
#       sprinkle(cols=1:7,rows=1:{length(x3)-5},border=c("left","right"))%>%
#       sprinkle(cols=1:7,rows=1,border=c("top","bottom","left","right"),part="head")%>%
#       sprinkle(cols=1,rows=1,border="left",part="head")%>%
#       sprinkle(cols=7,rows=1,border="right",part="head")%>%
#       sprinkle_width(cols=1,rows=1:2,width=90,width_units="pt")%>%
#       sprinkle_width(cols=2,rows=1:2,width=30,width_units="pt")%>%
#       sprinkle_width(cols=3,width=60,width_units="pt")%>%
#       sprinkle_width(cols=4,width=60,width_units="pt")%>%
#       sprinkle_width(cols=5,width=50,width_units="pt")%>%
#       sprinkle_width(cols=6,width=50,width_units="pt")%>%
#       sprinkle_width(cols=7,width=70,width_units="pt")%>%
#       sprinkle(rows=1,halign="center",part="head")
#
#
#     if(pix.int){
#       return(my.manova.pixie)
#     }else{
#       my.manova.pixie=print(my.manova.pixie,quote=F)[1]
#       return(my.manova.pixie)
#     }
#   }else{
#     #### Use car::Anova to get SS Type 3
#
#     my.summary=summary(my.model)
#     my.coefficients=my.summary$coefficients
#     my.coefficients=as.data.frame(my.coefficients)
#     my.III.summary=car::Anova(my.model,type=SS.type)
#     if(type=="glm" & is.null(my.factor)){
#       the.length=dim(my.III.summary)[1]+1
#     }else{
#       the.length=dim(my.III.summary)[1]
#     }
#
#
#     if(type=="lm"){
#       #### Calculate total SS
#       my.total=sum(my.III.summary$`Sum Sq`[2:length(my.III.summary$`Sum Sq`)])
#       my.df.total=sum(my.III.summary$Df[2:length(my.III.summary$Df)])
#       total.intercepts=1
#       my.rownames=c(rownames(my.summary$coefficients),"Residuals","Total")
#     }else if(type=="glm"){
#       #### Calculate model deviance stats
#       ddeviance2=my.model$null.deviance-my.model$deviance
#       ddf2=my.model$df.null-my.model$df.residual
#       fit2=1-pchisq(ddeviance2,ddf2)
#       total.intercepts=1
#       my.rownames=c(rownames(my.summary$coefficients),"Total")
#     }else if(type=="ord"){
#       my.temp.ord=update(my.model,~1)
#       ddeviance2=my.model$logLik-my.temp.ord$logLik
#       ddf2=my.model$edf-my.temp.ord$edf
#       fit2=1-pchisq(ddeviance2,ddf2)
#       total.intercepts=my.temp.ord$edf
#       my.rownames=c(rownames(my.summary$coefficients),"Total")
#     }else{
#       print("Error")
#       return()
#     }
#
#
#
#     my.std.error=c(my.coefficients$`Std. Error`,NA)
#     my.estimate=c(my.coefficients$Estimate,NA)
#     my.tables.df=as.data.frame(matrix(ncol=7,nrow=1))
#
#     if(type=="lm"){
#       names(my.tables.df)=c("rownames","sumsq","df","est","std.err","f.val","p.val")
#     }else{
#       names(my.tables.df)=c("var","od.rat","std.err","z.val","dev","df","p.val")
#     }
#
#     #### Make the double table entries ####
#
#     #### Was very annoying...took my frustration out on
#     #### Variable names
#
#     factor.stupid=NULL
#     factor.rownames=NULL
#     num.of.levels=NULL
#     ordinal.temp=0
#     if(!is.null(my.factor)){
#
#       for(i in 1:length(my.factor)){
#         factor.stupid=c(factor.stupid,grep(paste("^",my.factor[i],"$",sep=""),names(myDF)))
#
#         if(type=="lm"){
#           factor.rownames=c(factor.rownames,grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary)))
#           num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
#           # }else if(type=="glm"){
#           #   factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+1})
#           #   num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
#         }else{
#           factor.rownames=c(factor.rownames,{grep(paste("^",my.factor[i],"$",sep=""),rownames(my.III.summary))+total.intercepts+sum(num.of.levels)-ordinal.temp})
#           num.of.levels=c(num.of.levels,length(levels(myDF[[factor.stupid[i]]])))
#           ordinal.temp=ordinal.temp+1
#         }
#       }
#
#     }else{
#
#       factor.rownames=0L
#
#     }
#
#
#     my.factor.var=1
#     this.temp.var=1
#     this.shift.temp=1
#     yet.another.var=1
#     my.shift=0
#     other.temp=2
#     ord.temp=1
#     if(type=="lm"){
#       dang.length=length(rownames(my.III.summary))
#     }else if (type=="glm"){
#       dang.length=total.intercepts+length(rownames(my.III.summary))+max({sum(num.of.levels)-length(my.factor)},0)+1
#     }else{
#       dang.length=total.intercepts+length(rownames(my.III.summary))+max({sum(num.of.levels)-length(my.factor)},0)+1
#       #dang.length=7
#     }
#
#     while(this.shift.temp<dang.length){
#
#       if(is.na(factor.rownames[my.factor.var])){
#
#         my.factor.rownames=1
#
#       }else{
#
#         my.factor.rownames=factor.rownames[my.factor.var]
#
#       }
#       if(this.shift.temp==1){
#         i=1
#         while(i<=total.intercepts){
#           if(type=="lm"){
#             my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
#             my.df=my.III.summary$Df[this.shift.temp]
#             my.est=my.estimate[this.shift.temp]
#             my.std.err=my.std.error[this.shift.temp]
#             my.f.val=my.III.summary$`F value`[this.shift.temp]
#             my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
#             my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val)
#           }else{
#             my.or=exp(my.estimate[this.shift.temp])
#             my.est=my.estimate[this.shift.temp]
#             my.z.val=my.summary$coefficients[this.shift.temp,3]
#             my.std.err=my.std.error[this.shift.temp]
#             #my.dev=my.III.summary$`LR Chisq`[this.shift.temp]
#             #my.df=my.III.summary$Df[this.shift.temp]
#             my.p.val=my.summary$coefficients[{this.shift.temp},4]
#             my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
#           }
#           this.shift.temp=this.shift.temp+1
#           this.temp.var=this.temp.var+1
#           i=i+1
#         }
#       }else if(this.shift.temp==my.factor.rownames){
#         if(type=="lm"){
#           my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
#           my.df=my.III.summary$Df[this.shift.temp]
#           my.est=NA
#           my.std.err=NA
#           my.z.val=NA
#           my.f.val=my.III.summary$`F value`[this.shift.temp]
#           my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
#           my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.sumsq,my.df,my.std.err,my.z.val,my.f.val,my.p.val)
#         }else if(type=="glm"){
#           my.or=NA
#           my.est=NA
#           my.std.err=NA
#           my.dev=my.III.summary$`LR Chisq`[ord.temp]
#           my.df=my.III.summary$Df[ord.temp]
#           my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
#           my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.or,my.est,my.std.err,my.dev,my.df,my.p.val)
#           ord.temp=ord.temp+1
#         }else{
#           my.or=NA
#           my.est=NA
#           my.std.err=NA
#           my.dev=my.III.summary$Chisq[ord.temp]
#           my.df=my.III.summary$Df[ord.temp]
#           my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
#           my.tables.df[this.temp.var,]=c(my.factor[yet.another.var],my.or,my.est,my.std.err,my.dev,my.df,my.p.val)
#           ord.temp=ord.temp+1
#
#         }
#         yet.another.var=yet.another.var+1
#         this.temp.var=this.temp.var+1
#         this.shift.temp=this.shift.temp+1
#         other.other.temp=2
#
#
#         if(length(grepl(":",my.summary$coefficients[other.temp,1]))>0){
#
#           while(other.other.temp<{num.of.levels[my.factor.var]+1}){
#             if(type=="lm"){
#               my.sumsq=NA
#               my.df=NA
#               my.est=my.estimate[other.temp]
#               my.std.err=my.std.error[other.temp]
#               #### NEED TO FIX ####
#               my.f.val={my.summary$coefficients[other.temp,3]^2}
#               my.p.val=my.summary$coefficients[other.temp,4]
#               my.tables.df[this.temp.var,]=c(my.rownames[other.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val)
#             }else if(type=="glm"){
#               my.or=exp(my.estimate[other.temp])
#               my.est=my.estimate[other.temp]
#               my.std.err=my.std.error[other.temp]
#               my.z.val=my.summary$coefficients[other.temp,3]
#               my.dev=my.III.summary$`LR Chisq`[other.temp]
#               my.df=my.III.summary$Df[other.temp]
#               my.p.val=my.summary$coefficients[other.temp,4]
#               my.tables.df[this.temp.var,]=c(my.rownames[other.temp],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
#               this.shift.temp=this.shift.temp+1
#             }else{
#               my.or=exp(my.estimate[other.temp+total.intercepts-1])
#               my.est=my.estimate[other.temp+total.intercepts-1]
#               my.std.err=my.std.error[other.temp+total.intercepts-1]
#               my.z.val=my.summary$coefficients[{other.temp+total.intercepts-1},3]
#               #my.dev=my.III.summary$Chisq[other.temp]
#               #my.df=my.III.summary$Df[other.temp]
#               my.p.val=my.summary$coefficients[{other.temp+total.intercepts-1},4]
#               my.tables.df[this.temp.var,]=c(my.rownames[other.temp+total.intercepts-1],my.or,my.std.err,my.z.val,NA,NA,my.p.val)
#               this.shift.temp=this.shift.temp+1
#             }
#             this.temp.var=this.temp.var+1
#             other.temp=other.temp+1
#             other.other.temp=other.other.temp+1
#             the.length=the.length+1
#
#           }
#
#         }else{
#
#         }
#
#         if(my.factor.var==1){
#
#           my.shift={my.shift+other.temp-2}
#
#         }else{
#
#           my.shift=my.shift+other.temp
#
#         }
#
#         my.factor.var=my.factor.var+1
#
#       }else{
#         if(type=="lm"){
#           my.sumsq=my.III.summary$`Sum Sq`[this.shift.temp]
#           my.df=my.III.summary$Df[this.shift.temp]
#           my.est=my.estimate[this.shift.temp]
#           my.std.err=my.std.error[this.shift.temp]
#           my.f.val=my.III.summary$`F value`[this.shift.temp]
#           my.p.val=my.III.summary$`Pr(>F)`[this.shift.temp]
#           my.tables.df[this.temp.var,]=c(rownames(my.III.summary)[this.shift.temp],my.sumsq,my.df,my.est,my.std.err,my.f.val,my.p.val)
#         }else if(type=="glm"){
#           if(!is.null(my.factor)){
#             this.shift.temp=this.shift.temp-ordinal.temp
#           }
#           my.or=exp(my.estimate[this.shift.temp])
#           my.est=my.estimate[this.shift.temp]
#           my.std.err=my.std.error[this.shift.temp]
#           my.z.val=my.summary$coefficients[other.temp,3]
#           my.dev=my.III.summary$`LR Chisq`[ord.temp]
#           my.df=my.III.summary$Df[ord.temp]
#           my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
#           my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,NA,my.dev,my.df,my.p.val)
#           if(!is.null(my.factor)){
#             this.shift.temp=this.shift.temp+ordinal.temp
#           }
#           ord.temp=ord.temp+1
#         }else{
#           if(!is.null(my.factor)){
#             this.shift.temp=this.shift.temp-ordinal.temp
#           }
#           my.or=exp(my.estimate[this.shift.temp])
#           my.est=my.estimate[this.shift.temp]
#           my.std.err=my.std.error[this.shift.temp]
#           my.z.val=my.summary$coefficients[this.shift.temp,3]
#           my.dev=my.III.summary$Chisq[ord.temp]
#           my.df=my.III.summary$Df[ord.temp]
#           my.p.val=my.III.summary$`Pr(>Chisq)`[ord.temp]
#           my.tables.df[this.temp.var,]=c(my.rownames[this.shift.temp],my.or,my.std.err,NA,my.dev,my.df,my.p.val)
#           if(!is.null(my.factor)){
#             this.shift.temp=this.shift.temp+ordinal.temp
#           }
#           ord.temp=ord.temp+1
#         }
#         this.shift.temp=this.shift.temp+1
#         this.temp.var=this.temp.var+1
#
#       }
#     }
#
#     if(type=="lm"){
#       my.tables.df[this.temp.var,]=c("Residuals",my.III.summary$`Sum Sq`[this.shift.temp],my.III.summary$Df[this.shift.temp],NA,NA,NA,NA)
#       my.tables.df[this.temp.var+1,]=c("Total",my.total,my.df.total,NA,NA,NA,NA)
#       my.tables.df$f.val=as.numeric(my.tables.df$f.val)
#       my.tables.df$est=as.numeric(my.tables.df$est)
#       my.tables.df$sumsq=as.numeric(my.tables.df$sumsq)
#     }else{
#       my.tables.df[this.temp.var,]=c("Change from Null",NA,NA,NA,ddeviance2,ddf2,fit2)
#       my.tables.df$od.rat=as.numeric(my.tables.df$od.rat)
#       my.tables.df$z.val=as.numeric(my.tables.df$z.val)
#       my.tables.df$dev=as.numeric(my.tables.df$dev)
#     }
#
#     my.tables.df$df=as.numeric(my.tables.df$df)
#     my.tables.df$std.err=as.numeric(my.tables.df$std.err)
#     my.tables.df$p.val=as.numeric(my.tables.df$p.val)
#
#     #### Make custom glance stats ####
#
#     #### Can eventually make it options
#     if(type=="lm"){
#       glance_stats=broom::glance(my.model)
#       glance_stats=tidyr::gather(glance_stats)
#
#
#       glance_stats[3:5]=NA
#       glance_stats[6]=c(glance_stats$key[7:8],NA,glance_stats$key[9:11],NA,NA,NA,NA,NA)
#       glance_stats[7]=c(glance_stats$value[7:8],NA,glance_stats$value[9:11],NA,NA,NA,NA,NA)
#       glance_stats=glance_stats[-7:-11,]
#       glance_stats=glance_stats[-3,]
#       glance_stats=glance_stats[-5,]
#     }else{
#       # glance_stats2=as.data.frame(matrix(ncol=7,nrow=1))
#       # glance_stats2[1,]=c(glance_stats$key[1],glance_stats$value[1],NA,NA,NA,glance_stats$key[6],glance_stats$value[6])
#       # glance_stats2[2,]=c(glance_stats$key[2],glance_stats$value[2],NA,NA,NA,glance_stats$key[7],glance_stats$value[7])
#       # glance_stats2[3,]=c(glance_stats$key[3],glance_stats$value[3],NA,NA,NA,glance_stats$key[4],glance_stats$value[4])
#       # glance_stats=glance_stats2
#     }
#     #### For total
#     the.length=the.length+1
#     this.temp.var=this.temp.var+1
#     #### Make table ####
#
#     options(pixie_interactive = pix.int,pixie_na_string="")
#
#
#     if(type=="lm"){
#       my.dust=pixiedust::dust(my.tables.df)%>%
#         sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
#         sprinkle_print_method(pix.method)%>%
#         sprinkle_na_string()%>%
#         sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
#         sprinkle(rows=1:the.length,cols=1:7,
#                  border="right",border_color="black")%>%
#         sprinkle(rows=1:the.length,cols=1,
#                  border="left",border_color="black")%>%
#         sprinkle(rows=1,cols=1:7,
#                  border=c("top","bottom","left","right"),border_color="black",part="head")%>%
#         sprinkle(rows=this.temp.var,cols=1:7,
#                  border="bottom",border_color="black")%>%
#         sprinkle(cols=c("sumsq","est","std.err","f.val"),round=2)%>%
#         sprinkle_colnames("Variable","Sums of Squares","df","Estimate","Std. Error","F-value","Pr(>F)")
#
#     }else if(type=="glm"){
#       my.dust=pixiedust::dust(my.tables.df)%>%
#         sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
#         sprinkle_print_method(pix.method)%>%
#         sprinkle_na_string()%>%
#         sprinkle(cols=2:5,round=2)%>%
#         sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
#         sprinkle(cols=2:7,pad=5,part="head")%>%
#         sprinkle(rows={this.temp.var-2},cols=1:7,border="bottom",border_color="black")%>%
#         sprinkle(rows=1:the.length,cols=1:7,
#                  border="right",border_color="black")%>%
#         sprinkle(rows=1:the.length,cols=1,
#                  border="left",border_color="black")%>%
#         sprinkle(rows=1,cols=1:7,
#                  border=c("top","bottom","left","right"),border_color="black",part="head")%>%
#         sprinkle(rows=this.temp.var-1,cols=1:7,
#                  border="bottom",border_color="black")%>%
#         sprinkle_colnames("Variable","Odds Ratio","Std. Error","z-Value","Deviance","df","p-Value")
#
#     }else{
#       my.dust=pixiedust::dust(my.tables.df)%>%
#         sprinkle(cols="p.val",fn=quote(pvalString(value,digits=3,format="default")))%>%
#         sprinkle_print_method(pix.method)%>%
#         sprinkle_na_string()%>%
#         sprinkle(cols=2:5,round=2)%>%
#         sprinkle(cols=c(2,3,4,5,6,7),pad=5)%>%
#         sprinkle(cols=2:7,pad=5,part="head")%>%
#         sprinkle(rows=total.intercepts,cols=1:7,border="bottom",border_color="black")%>%
#         sprinkle(rows={this.temp.var-2},cols=1:7,border="bottom",border_color="black")%>%
#         sprinkle(rows=1:{the.length+total.intercepts},cols=1:7,
#                  border="right",border_color="black")%>%
#         sprinkle(rows=1:{the.length+total.intercepts},cols=1,
#                  border="left",border_color="black")%>%
#         sprinkle(rows=1,cols=1:7,
#                  border=c("top","bottom","left","right"),border_color="black",part="head")%>%
#         sprinkle(rows=this.temp.var-1,cols=1:7,
#                  border="bottom",border_color="black")%>%
#         sprinkle_colnames("Variable","Odds Ratio","Std. Error","z-Value","Deviance","df","p-Value")
#     }
#
#     if(type=="lm"){
#       my.dust=pixiedust::redust(my.dust,glance_stats,part="foot")%>%
#         sprinkle(cols=c(2,7),round=3,part="foot")%>%
#         sprinkle(cols=3:5,replace=c("","","","","","","","","","","",""),part="foot")%>%
#         sprinkle(cols=1, replace=c("R-Square","Adj R-Sq","F-Statistic","P-Value"),part="foot")%>%
#         sprinkle(cols=2,rows=4,fn=quote(pvalString(value,digits=3,format="default")),part="foot")%>%
#         sprinkle(cols=1:7,rows=1,halign="center",part="head")%>%
#         sprinkle_width(cols=1,width=90,width_units="pt")%>%
#         sprinkle_width(cols=2,width=108,width_units="pt")%>%
#         sprinkle_width(cols=4,width=62,width_units="pt")%>%
#         sprinkle_width(cols=5,width=68,width_units="pt")%>%
#         sprinkle_width(cols=6,width=68,width_units="pt")%>%
#         sprinkle_width(cols=7,width=71,width_units="pt")%>%
#         sprinkle(cols=2,halign="left",part="foot")
#     }else{
#       # my.dust=pixiedust::redust(my.dust,glance_stats,part="foot")%>%
#       #   sprinkle(cols=c(2,7),round=2,part="foot")%>%
#       #   sprinkle(cols=3:5,replace=c("","","","","","","","",""),part="foot")
#     }
#     if(pix.int){
#       return(my.dust)
#     }else{
#       my.dust.print=print(my.dust,quote=F)[1]
#       return(my.dust.print)
#     }
#
#   }
# }

#### Formula dissecter ####

library(gtools)
#### Split and get pieces
for(i in 1:length(my.model$call)){
  my.formula.grep=grep("~",my.model$call[[i]])
  if(length(my.formula.grep)>0){

    #### Get vars
    my.vars=strsplit(as.character(my.model$call[[i]]),"~")
    my.dep.var=my.vars[[2]]
    my.ind.vars=strsplit(my.vars[[3]],"\\+")

    #### Null formula
    my.null.formula=paste(my.dep.var,"~1",sep="")

    #### Other formulas
    my.ind.vars.perm=NULL
    if(length(my.ind.vars[[1]])>2){

      my.ind.vars.perm[[1]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]])-1,my.ind.vars[[1]])
      my.ind.vars.perm[[2]]=gtools::permutations(length(my.ind.vars[[1]]),length(my.ind.vars[[1]]),my.ind.vars[[1]])
    }else{
      my.ind.vars.perm[[1]]=as.data.frame(matrix(ncol=1,nrow=2))
      my.ind.vars.perm[[1]]=rbind(names(my.model$model)[2],names(my.model$model)[1])
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
    my.pre.formula.list=c(my.pre.formula.list,my.temp.form)
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
    my.formula.list=c(my.formula.list,my.temp.form)
  }
}


my.return.list=list(my.null.formula,my.pre.formula.list,my.formula.list)

return(my.return.list)

my.vars=strsplit(as.character(my.formula),"~")





my.vars=strsplit(as.character(new.model$formula),"~")
my.dep.var=my.vars[[2]]
my.other.vars=names(new.model$model)
my.var.grep=grep(paste("^",my.dep.var,"$",sep=""),my.other.vars)
my.other.vars=my.other.vars[-my.var.grep]
my.formula=paste("~",my.other.vars[1])
for(i in 2:length(my.other.vars)){
  my.formula=paste(my.formula,"*",my.other.vars[i],sep="")}





my.new.model=NULL
my.new.model[[1]]=my.null.model


#### Type I & II Dataframes ####
j=dim(my.new.df)[2]
line=2
for(i in 2:{dim(my.new.df)[2]}){
  j=i
  my.rh=1
  k=1
  while(k <{dim(my.new.df)[2]}){
    if(j>{dim(my.new.df)[2]}){
      my.j=abs(j-dim(my.new.df)[2])+1
    }else{
      my.j=j
    }
    my.rh=paste(my.rh,"+",names(my.new.df)[my.j])
    my.new.model[[line]]=manova(eval(parse(text=paste("my.new.df[[1]]~",my.rh))),data=my.new.df)

    j=j+1
    line=line+1
    k=k+1
  }
}
my.SS.type.1=NULL
my.SS.type.1.change=NULL


#### Type I SS
#### Also treatment change
#1+n
my.SS.type.1.dfs=NULL
the.var.length=dim(my.new.df)[2]
for(i in 1:the.var.length){
  my.SS.type.1.dfs[[i]]=my.new.model[[i]]
}




#### Null change
my.SS.type.1.change[[1]]=summary(my.SS.type.1.dfs[[2]])$SS[[1]]
for(z in 2:the.var.length-1){
  my.SS.type.1.change[[z]]=summary(my.SS.type.1.dfs[[z+1]])$SS[[z]]
}


#### Make treatment
my.SS.type.1.change.total=0
for(y in 1:length(my.SS.type.1.change)){
  my.SS.type.1.change.total=my.SS.type.1.change.total+my.SS.type.1.change[[y]]
}


#### Type II SS
#they are type 2, here is proof.
my.SS.type.2=NULL
my.SS.type.2.change=NULL
my.SSP.type.2.change.fa=NULL
my.latent.SSP.type.2.change=NULL
my.SSP.type.2.treat.fa=NULL
my.latent.SSP.type.2.treat=NULL


#### Type II SS
the.var.length=dim(my.new.df)[2]-1
kappa=the.var.length
#3n+1
my.SS.type.2.change[[1]]=summary(my.new.model[[kappa+1]])$SS[the.var.length]








#### Make latent ####
#### Make latent variables (~equivalent to ANOVAs) of Full model (i.e. type I)
my.SSP.type.2.change.fa.eigen=eigen(my.SS.type.2.change[[1]][[1]])
my.SSP.type.2.change.fa.values=1/sqrt(my.SSP.type.2.change.fa.eigen$values)
for(n in 2:length(my.SSP.type.2.change.fa.values)){
  my.SSP.type.2.change.fa.values=rbind(my.SSP.type.2.change.fa.values,my.SSP.type.2.change.fa.values)
}
my.SSP.type.2.change.fa.values.t=t(my.SSP.total.fa.values)

my.SSP.type.2.change.fa[[1]]=my.SSP.type.2.change.fa.eigen$vectors%*%my.SSP.type.2.change.fa.values.t
my.latent.SSP.type.2.change[[1]]=my.SS.type.2.change[[1]][[1]]%*%my.SSP.type.2.change.fa[[1]]

#### Make residual totals ####
for(l in 2:{dim(my.new.df)[2]-1}){
  my.gamma=kappa*l+1
  my.SS.type.2.change[[l]]=summary(my.new.model[[my.gamma]])$SS[the.var.length]
  my.SSP.type.2.change.fa.eigen=eigen(my.SS.type.2.change[[l]][[1]])
  my.SSP.type.2.change.fa.values=1/sqrt(my.SSP.type.2.change.fa.eigen$values)
  for(n in 2:length(my.SSP.type.2.change.fa.values)){
    my.SSP.type.2.change.fa.values=rbind(my.SSP.type.2.change.fa.values,my.SSP.type.2.change.fa.values)
  }
  my.SSP.type.2.change.fa.values.t=t(my.SSP.total.fa.values)

  my.SSP.type.2.change.fa[[l]]=my.SSP.type.2.change.fa.eigen$vectors%*%my.SSP.type.2.change.fa.values.t
  my.latent.SSP.type.2.change[[l]]=my.SS.type.2.change[[l]][[1]]%*%my.SSP.type.2.change.fa[[l]]
}


#### Get residual stuff
the.resid.SS=sum(diag(my.SSP.err))
the.resid.df=my.model$df.residual



#### Make SSCPs  of full matrix####
#### Get SSP treatment and error

my.SSP.treat=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSP
my.SSP.err=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSPE

### Get total SSP
my.SSP.treat.total=my.SSP.treat[[1]]
for(i in 2:{length(my.SSP.treat)}){
  my.SSP.treat.total=my.SSP.treat.total+my.SSP.treat[[i]]
}
my.SSP.total=my.SSP.treat.total+my.SSP.err


#### Get all df ####
my.SSP.total.df=my.y.levels*dim(my.model$model)[1]-1
my.SSP.err.df=my.model$df.residual
my.SSP.treat.df.total=my.SSP.total.df-my.SSP.err.df
my.SSP.treat.df=1
for(i in 2:length(my.SSP.treat)){
  manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
  if(length(manova.grep)>0){
    my.SSP.treat.df=c(my.SSP.treat.df,{length(my.model$xlevels[[manova.grep]])-1})
  }else{
    my.SSP.treat.df=c(my.SSP.treat.df,1)
  }
}


# #### MAKE REGULAR AND LATENT CONTRASTS ####
# for(i in 1:{my.y.levels}){
#   #### Pick right model
#   if(i==1){
#     my.contrast.model=my.new.model[[my.y.levels^2+1]]
#   }else{
#     my.contrast.model=my.new.model[[{i-1}*my.y.levels+1]]
#   }
#
#
#   #### Mean Square Error
#   my.MSE=mean(my.contrast.model$residuals^2)
#
#
#   #### Latent Mean Square Error
#   my.latent.MSE=NULL
#   for(j in 1:my.y.levels){
#     if(i==1){
#       my.latent.MSE=as.numeric(mean(my.contrast.model$residuals[i]^2))
#     }else{
#       my.latent.MSE=c(my.latent.MSE,as.numeric(mean(my.contrast.model$residuals[i]^2)))
#     }
#   }
#
#
#
#   #### Get mean responses for variables longer than 2
#   #### WARN! LEVEL NAMES IS GOING TO BREAK IT!!!! ####
#   my.count.means=NULL
#   my.count.n=NULL
#   level.names=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       if(q==2){
#         level.names=as.vector(levels(my.new.df[[q]]))
#       }else{
#         level.names=rbind(level.names,as.vector(levels(my.new.df[[q]])))
#       }
#       for(j in 1:{my.SSP.treat.df[q]+1})
#         if(j==1){
#           my.count.means[[q-1]]=as.list(mean(my.new.df[as.character(my.new.df[[q]])==level.names[j],1]))
#           my.count.n[[q-1]]=as.list(dim(my.new.df[as.character(my.new.df[[q]])==level.names[j],1])[1])
#         }else{
#           my.count.means[[q-1]]=c(my.count.means[[q-1]],mean(my.new.df[as.character(my.new.df[[q]])==level.names[j],1]))
#           my.count.n[[q-1]]=c(my.count.n[[q-1]],dim(my.new.df[as.character(my.new.df[[q]])==level.names[j],1])[1])
#         }
#     }
#   }
#
#   #### Make contrasts
#   my.contrasts=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       level.names=levels(my.new.df[[q]])
#       for(j in 1:{my.SSP.treat.df[q]})
#         if(j==1){
#           my.contrasts[[q-1]]=c(1,-1,rep(0,my.SSP.treat.df[q]-1))
#         }else{
#           my.contrasts[[q-1]]=rbind(my.contrasts[[q-1]],c(1,rep(0,j-1),-1,rep(0,my.SSP.treat.df[q]-j)))
#         }
#     }
#   }
#
#   #### Compute F values
#   my.contrasts.F=NULL
#   my.latent.contrasts.F=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       for(j in 1:{my.SSP.treat.df[q]})
#         if(j==1){
#           my.contrasts.I=as.integer(t(as.matrix(my.count.means[[q-1]])))%*%as.integer(as.matrix(my.contrasts[[q-1]][j,]))
#           my.contrasts.denom=sum(my.contrasts[[q-1]][j,]^2/as.integer(my.count.n[[q-1]]))
#           my.contrasts.F[[q-1]]=as.list({my.contrasts.I^2}/{my.latent.MSE[j]*my.contrasts.denom})
#           for(l in 1:{my.y.levels}){
#             if(l==1){
#               my.latent.contrasts.F[[q-1]]=as.list({my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
#             }else{
#               my.latent.contrasts.F[[q-1]]=c(my.latent.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
#             }
#           }
#         }else{
#           my.contrasts.I=as.integer(t(as.matrix(my.count.means[[q-1]])))%*%as.integer(as.matrix(my.contrasts[[q-1]][j,]))
#           my.contrasts.denom=sum(my.contrasts[[q-1]][j,]^2/as.integer(my.count.n[[q-1]]))
#           my.contrasts.F[[q-1]]=c(my.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
#           for(l in 1:{my.y.levels}){
#             my.latent.contrasts.F[[q-1]]=c(my.latent.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
#           }
#         }
#     }
#   }
#
#   #### Compute SS from F values
#   #### NEED TO FIX
#   #### Fixed
#   #### change to partial resid.SS
#
#   my.contrasts.SS=NULL
#   my.latent.contrasts.SS=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       for(j in 1:{my.SSP.treat.df[q]}){
#         if(j==1){
#           my.contrasts.SS[[q-1]]=as.list({as.numeric(my.contrasts.F[[q-1]][j])*{the.resid.SS}*my.y.levels}/{the.resid.df})
#           for(l in 1:{my.y.levels}){
#             if(l==1){
#               my.latent.contrasts.SS[[q-1]]=as.list({as.numeric(my.latent.contrasts.F[[q-1]][1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
#             }else{
#               my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][my.y.levels*l-1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
#             }
#           }
#         }else{
#           my.contrasts.SS[[q-1]]=c(my.contrasts.SS[[q-1]],{as.numeric(my.contrasts.F[[q-1]][j])*{the.resid.SS}*my.y.levels}/{the.resid.df})
#           for(l in 1:{my.y.levels}){
#             if(l==1){
#               my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][j])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
#             }else{
#               my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][my.y.levels*j-1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
#             }
#           }
#         }
#       }
#     }
#   }
#
#
#   #### Make rownames
#   my.contrasts.names=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       for(j in 1:{my.SSP.treat.df[q]}){
#         if(j==1){
#           my.contrasts.names[[q-1]]=as.list(paste(level.names[[1]],"-",level.names[[j+1]]))
#         }else{
#           my.contrasts.names[[q-1]]=c(my.contrasts.names[[q-1]],paste(level.names[[1]],"-",level.names[[j+1]]))
#         }
#       }
#     }
#   }
#
#
#   #### Make table
#   #### NEED TO ADD PVAL SO CAN ADJUST ####
#   #### p.adjust works on set of p-vals!
#   my.contrasts.table=NULL
#   for(q in 2:length(my.SSP.treat.df)){
#     if(my.SSP.treat.df[q]>1){
#       num.of.contrasts=my.SSP.treat.df[q]
#       my.contrasts.table[[q-1]]=cbind(my.contrasts.names[[1]],my.contrasts.F[[1]],my.contrasts.SS[[1]])
#       colnames(my.contrasts.table[[q-1]])=c("name","F.val","SS")
#     }
#   }
# }
#
#
# #### Latent tables
# my.latent.contrasts.F.R=NULL
# my.latent.contrasts.SS.R=NULL
# for(q in 2:length(my.SSP.treat.df)){
#   if(my.SSP.treat.df[q]>1){
#     if(q==2){
#       the.latent.levels=1
#       num.of.contrasts=my.SSP.treat.df[q]
#       for(s in 2:my.SSP.treat.df[q]){
#         the.latent.levels=c(the.latent.levels,my.y.levels*s-1)}
#       my.latent.contrasts.F.R[[q-1]]=as.list(my.latent.contrasts.F[[q-1]][[1]])
#       my.latent.contrasts.SS.R[[q-1]]=as.list(my.latent.contrasts.SS[[q-1]][[1]])
#       for(r in 2:length(the.latent.levels)){
#         my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
#         my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
#       }
#       for(v in 2:my.y.levels){
#         the.latent.levels=NULL
#         for(s in 1:my.SSP.treat.df[q]){
#           the.latent.levels=c(the.latent.levels,my.y.levels*s)}
#         for(r in 1:length(the.latent.levels)){
#           my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
#           my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
#         }
#       }
#     }else{
#       for(v in 1:my.y.levels){
#         the.latent.levels=NULL
#         #the.latent.levels=q-1
#         for(s in 1:my.SSP.treat.df[q]){
#           the.latent.levels=c(the.latent.levels,my.y.levels*s)}
#         for(r in 1:length(the.latent.levels)){
#           my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
#           my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
#         }
#       }
#     }
#   }
# }
# my.latent.contrasts.table=NULL
# for(q in 2:length(my.SSP.treat.df)){
#   num.of.contrasts=my.SSP.treat.df[q]
#   if(my.SSP.treat.df[q]>1){
#     for(g in 1:my.y.levels){
#       my.latent.contrasts.table[[q-1]][[g]]=cbind(my.contrasts.names[[1]],my.latent.contrasts.F.R[[q-1]][{g+ifelse(g>1,{g-1}*my.SSP.treat.df[q]-1,0)}:{g*my.SSP.treat.df[q]}],my.latent.contrasts.SS.R[[q-1]][{g+ifelse(g>1,{g-1}*my.SSP.treat.df[q]-1,0)}:{g*my.SSP.treat.df[q]}])
#     }
#   }
# }


#### LATENT MODEL ####
#### Make latent variables (~equivalent to ANOVAs) of Full matrix
my.SSP.total.fa.eigen=eigen(my.SSP.total)
my.SSP.total.fa.values=1/sqrt(my.SSP.total.fa.eigen$values)
for(i in 2:length(my.SSP.total.fa.values)){
  my.SSP.total.fa.values=rbind(my.SSP.total.fa.values,my.SSP.total.fa.values)
}
my.SSP.total.fa.values.t=t(my.SSP.total.fa.values)
my.SSP.total.fa=my.SSP.total.fa.eigen$vectors%*%my.SSP.total.fa.values.t
my.latent.SSP.total=my.SSP.total%*%my.SSP.total.fa

my.SSP.treat.fa.eigen=NULL
my.SSP.treat.fa.values=NULL
my.SSP.treat.fa.values.t=NULL
my.SSP.treat.fa=NULL
my.latent.SSP.treat=NULL
for(i in 1:length(my.SSP.treat)){
  my.SSP.treat.fa.eigen[[i]]=eigen(my.SSP.treat[[i]])


  my.SSP.treat.fa.values[[i]]=tryCatch(as.matrix({1/sqrt(my.SSP.treat.fa.eigen[[i]]$values)}),
                                       warning=function(w){
                                         my.SSP.treat.fa.values[[i]]=matrix(ncol=1,nrow=length(my.SSP.treat.fa.eigen[[i]]$values))
                                         for(j in 1:length(my.SSP.treat.fa.eigen[[i]]$values)){
                                           my.SSP.treat.fa.values[[i]][j]={1/sqrt(max(my.SSP.treat.fa.eigen[[i]]$values[j],0))}
                                           if(my.SSP.treat.fa.values[[i]][j]==Inf){
                                             my.SSP.treat.fa.values[[i]][j]=0
                                           }
                                         }
                                         return(my.SSP.treat.fa.values[[i]])
                                       }
  )

  for(j in 2:length(my.SSP.treat.fa.values[[i]])){
    my.SSP.treat.fa.values[[i]]=cbind(my.SSP.treat.fa.values[[i]],my.SSP.treat.fa.values[[i]])
  }
  my.SSP.treat.fa[[i]]=my.SSP.treat.fa.eigen[[i]]$vectors%*%my.SSP.treat.fa.values[[i]]
  my.latent.SSP.treat[[i]]=my.SSP.treat[[i]]%*%my.SSP.treat.fa[[i]]
}





#### Get SSP treatment and error for Type II error SS ####

my.SSP.treat=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSP
my.SSP.err=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSPE




#### Get total SSP ####
my.SSP.treat.total=my.SSP.treat[[1]]
for(i in 2:{length(my.SSP.treat)}){
  my.SSP.treat.total=my.SSP.treat.total+my.SSP.treat[[i]]
}
my.SSP.total=my.SSP.treat.total+my.SSP.err







#### Make error matrices ####
my.SSP.err.fa.eigen=eigen(my.SSP.err)
my.SSP.err.fa.values=1/sqrt(my.SSP.err.fa.eigen$values)
for(i in 2:length(my.SSP.err.fa.values)){
  my.SSP.err.fa.values=rbind(my.SSP.err.fa.values,my.SSP.err.fa.values)
}
my.SSP.err.fa.values.t=t(my.SSP.err.fa.values)
my.SSP.err.fa=my.SSP.err.fa.eigen$vectors%*%my.SSP.err.fa.values.t

my.latent.SSP.err=my.SSP.err%*%my.SSP.err.fa











#### Get treatment change totals ####
my.SSP.treat.change.total=0
my.SSP.treat.change=as.data.frame(matrix(ncol=my.y.levels,nrow=my.y.levels))
my.latent.SSP.treat.change.total=0
my.latent.SSP.treat.change=as.data.frame(matrix(ncol=my.y.levels,nrow=my.y.levels))
for(i in 2:length(my.SSP.treat)){
  if(i==2){
    my.SSP.treat.change=my.SSP.treat[[2]]
    my.latent.SSP.treat.change=my.latent.SSP.treat[[2]]
  }else{
    my.SSP.treat.change=my.SSP.treat.change+my.SSP.treat[[i]]
    my.latent.SSP.treat.change=my.latent.SSP.treat.change+my.SSP.treat[[i]]
  }
}
my.SSP.treat.change.df=sum(my.SSP.treat.df[-1])*my.y.levels
my.latent.SSP.treat.change.df=my.SSP.treat.change.df/my.y.levels
my.SSP.treat.change.total=sum(diag(my.SSP.treat.change))
my.latent.SSP.treat.change.total=sum(diag(my.latent.SSP.treat.change))








#### Get total change stuff ####
my.total.change=my.SSP.err+my.SSP.treat.change
the.total.change.SS=the.resid.SS+sum(diag(my.SS.type.1.change.total))
the.total.change.df=the.resid.df+{my.SSP.treat.change.df/my.y.levels}

#### Make basic table ####
my.table.names=c("var","test.stat","f.val","SS","df","mult.df","resid df","p.val")


my.manova.table=as.data.frame(matrix(ncol=v.p.len,nrow=1))
names(my.manova.table)=my.table.names
my.line.var=1
for(i in 1:length(my.SSP.treat)){

  #### Put in basic line ####
  my.treat.err=solve(my.SSP.err)%*%my.SSP.treat[[i]]
  my.test.stat=quick.m.test(my.treat.err,test.stat)
  my.SS=quick.tr(my.SSP.treat[[i]])
  my.df=my.y.levels*my.SSP.treat.df[i]
  my.resid.df=ifelse(i==1,{my.y.levels*the.resid.df},min({my.SSP.err.df*my.SSP.treat.df[i]-my.SSP.treat.df[i]},my.y.levels*my.SSP.err.df-my.SSP.treat.df[i]))
  my.f.val={my.SS/my.df}/{the.resid.SS/my.resid.df}
  my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

  my.manova.table[my.line.var,]=c(names(my.SSP.treat)[i],my.test.stat,my.f.val,my.SS,ifelse(i==1,NA,my.df/my.y.levels),ifelse(i==1,my.y.levels,my.df),my.resid.df,my.p.val)
  my.line.var=my.line.var+1


  #### Put in decomposed factors (y constrasts) ####
  #### Intercept
  if({show.intercepts & i==1}){
    for(y in 1:my.y.levels){
      my.name=rownames(my.SSP.total)[y]
      my.test.stat=NA
      my.SS=my.SSP.treat[[i]][y,y]
      my.df=my.SSP.treat.df[i]
      my.resid.df=my.SSP.err.df
      my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
      my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

      my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
      my.line.var=my.line.var+1
    }
  }
  #### Other
  if({show.y.contrasts & i!=1}){
    for(y in 1:my.y.levels){
      my.name=paste(ifelse(real.names,rownames(my.SSP.total)[y],y),"|",names(my.SSP.treat)[i],sep="")
      my.test.stat=NA
      my.SS=my.SSP.treat[[i]][y,y]
      my.df=my.SSP.treat.df[i]
      my.resid.df=min(abs(my.y.levels*my.SSP.err.df-my.df*my.SSP.err.df-my.df),abs(my.y.levels*my.SSP.err.df-my.df))
      my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
      my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

      my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
      my.line.var=my.line.var+1


      #### Put in latents (ANOVA) ####
      my.counter=NULL
      for(r in 2:my.y.levels){
        my.counter=c(my.counter,r)
      }
      my.counter=c(my.counter,1)
      if(show.latent &i!=1){
        my.i.temp=my.counter[i-1]
        my.y=y
        my.name=NA
        my.df=my.SSP.treat.df[my.i.temp]
        my.test.stat=NA
        my.resid.df=my.SSP.err.df-my.df
        my.SS=my.latent.SSP.type.2.change[[my.i.temp]][my.y,1]
        my.f.val={my.SS/my.df}/{my.latent.SSP.err[my.y,my.y]/my.resid.df}
        my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

        my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
        my.line.var=my.line.var+1
      }

      #### Put in latent contrasts ####
      if(show.contrasts & show.latent & i!=1){
        #### Check length
        if({my.SSP.treat.df[i]>1}){
          other.manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
          for(k in 1:{my.SSP.treat.df[i]}){
            my.name=my.contrasts.table[[other.manova.grep]][k,1]
            my.f.val=my.latent.contrasts.table[[other.manova.grep]][[i-1]][k,2]
            my.SS=my.latent.contrasts.table[[other.manova.grep]][[i-1]][k,3]
            my.test.stat=NA
            my.df=1
            my.mult.df=NA
            my.resid.df=the.resid.df-my.y.levels+1
            my.p.val=NA

            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,p.adjust(my.p.val,adjustment))
            my.line.var=my.line.var+1
          }
        }
      }

      #### Put in contrasts ####
      #### Not right. Don't have it decomposed this way.
      if(show.contrasts & !show.latent & F){
        #### Check length
        if({my.SSP.treat.df[i]>1}){
          other.manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
          for(k in 1:{my.SSP.treat.df[i]}){
            my.name=my.contrasts.table[[other.manova.grep]][k,1]
            my.f.val=my.contrasts.table[[other.manova.grep]][k,2]
            my.SS=my.contrasts.table[[other.manova.grep]][k,3]
            my.test.stat=NA
            my.df=NA
            my.mult.df=my.y.levels
            my.resid.df=the.resid.df-my.y.levels+1
            my.p.val=pf(as.numeric(my.f.val),as.numeric(my.mult.df),as.numeric(my.resid.df),lower.tail = F)

            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,p.adjust(my.p.val,adjustment))
            my.line.var=my.line.var+1
          }
        }
      }
    }
  }

  #### Put in contrasts ####
  if(show.contrasts & !show.latent & !show.y.contrasts){
    #### Check length
    if({my.SSP.treat.df[i]>1}){
      other.manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
      for(k in 1:{my.SSP.treat.df[i]}){
        my.name=my.contrasts.table[[other.manova.grep]][k,1]
        my.f.val=my.contrasts.table[[other.manova.grep]][k,2]
        my.SS=my.contrasts.table[[other.manova.grep]][k,3]
        my.test.stat=NA
        my.df=NA
        my.mult.df=my.y.levels
        my.resid.df=the.resid.df-my.y.levels+1
        my.p.val=pf(as.numeric(my.f.val),as.numeric(my.mult.df),as.numeric(my.resid.df),lower.tail = F)

        my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,p.adjust(my.p.val,adjustment))
        my.line.var=my.line.var+1
      }
    }
  }


  #### Put in latents (ANOVA) ####
  my.counter=NULL
  for(r in 2:my.y.levels){
    my.counter=c(my.counter,r)
  }
  my.counter=c(my.counter,1)
  if(show.latent &i!=1 & !show.y.contrasts){
    my.i.temp=my.counter[i-1]
    for(y in 1:{my.y.levels}){
      my.y=y
      my.name=paste(ifelse(real.names,rownames(my.SSP.total)[y],y),"|",names(my.SSP.treat)[i],sep="")
      my.df=my.SSP.treat.df[my.i.temp]
      my.test.stat=NA
      my.resid.df=my.SSP.err.df-my.df
      my.SS=my.latent.SSP.type.2.change[[my.i.temp]][my.y,1]
      my.f.val={my.SS/my.df}/{my.latent.SSP.err[my.y,my.y]/my.resid.df}
      my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

      my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
      my.line.var=my.line.var+1

      #### Put in latent contrasts ####
      if(show.contrasts & !{show.contrasts & show.latent}){
        #### Check length
        if({my.SSP.treat.df[i]>1}){
          other.manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
          for(k in 1:{my.SSP.treat.df[i]}){
            my.name=my.contrasts.table[[other.manova.grep]][k,1]
            my.f.val=my.latent.contrasts.table[[other.manova.grep]][[y]][k,2]
            my.SS=my.latent.contrasts.table[[other.manova.grep]][[y]][k,3]
            my.test.stat=NA
            my.df=1
            my.mult.df=NA
            my.resid.df=the.resid.df-my.y.levels+1
            my.p.val=NA

            my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,my.df,my.mult.df,my.resid.df,p.adjust(my.p.val,adjustment))
            my.line.var=my.line.var+1
          }
        }
      }
    }
  }


  #### Put in treatment ####
  #### From Type I statistics
  if(i==1){
    my.test.stat=quick.m.test(my.SS.type.1.change.total,test.stat)
    my.SS=sum(diag(my.SS.type.1.change.total))
    my.df=my.SSP.treat.change.df
    #### Should really be a min statement, but for later...
    my.resid.df=my.y.levels*my.SSP.err.df
    my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
    my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

    my.manova.table[my.line.var,]=c("Treatment",my.test.stat,my.f.val,my.SS,{my.df/my.y.levels},my.df,my.resid.df,my.p.val)
    my.line.var=my.line.var+1

    #### Treatment Y Contrasts ####
    if(show.y.contrasts){
      for(b in 1:my.y.levels){
        my.test.stat=NA
        my.name=paste(ifelse(real.names,rownames(my.SSP.total)[b],b),"|Treatment",sep="")
        my.SS=my.SS.type.1.change.total[b,b]
        my.df=my.SSP.treat.change.df
        #### Should really be a min statement, but for later...
        my.resid.df=my.y.levels*my.SSP.err.df
        my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
        my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

        my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),my.test.stat,my.f.val,my.SS,{my.df/my.y.levels},NA,{my.resid.df/my.y.levels},my.p.val)
        my.line.var=my.line.var+1

        #### Treatment Latents ####
        if(show.latent){
          #### Show latent treatments (ANOVAs)
          #### Need to change latents to type II
          my.SS=sum(weighted.residuals(my.null.model)[,b]^2)-sum(weighted.residuals(my.model)[,b]^2)
          my.df=my.SSP.treat.change.df/my.y.levels
          #### Should really be a min statement, but for later...
          my.resid.df=my.SSP.err.df-my.df
          my.f.val={my.SS/my.df}/{my.latent.SSP.err[y,y]/my.resid.df}
          my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
          my.manova.table[my.line.var,]=c(NA,NA,my.f.val,my.SS,my.df,NA,{my.resid.df},my.p.val)
          my.line.var=my.line.var+1
        }
      }
    }

    #### Latents ####
    my.counter=NULL
    for(r in 2:my.y.levels){
      my.counter=c(my.counter,r)
    }
    my.counter=c(my.counter,1)
    if(show.latent & !show.y.contrasts){
      for(b in 1:{my.y.levels}){
        my.name=paste(ifelse(real.names,rownames(my.SSP.total)[b],b),"|Treatment",sep="")
        my.SS=sum(weighted.residuals(my.null.model)[,b]^2)-sum(weighted.residuals(my.model)[,b]^2)
        my.df=my.SSP.treat.change.df/my.y.levels
        #### Should really be a min statement, but for later...
        my.resid.df=my.SSP.err.df-my.df
        my.f.val={my.SS/my.df}/{my.latent.SSP.err[y,y]/my.resid.df}
        my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
        my.manova.table[my.line.var,]=c(abbreviate(my.name,abbrev.length),NA,my.f.val,my.SS,my.df,NA,{my.resid.df},my.p.val)
        my.line.var=my.line.var+1
      }
    }
  }
}

