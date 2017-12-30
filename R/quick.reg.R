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

quick.reg.table = function(my.model,
                           part.eta = F,
                           VIF = F,
                           myDF = my.found.df,
                           SS.type = 2,
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



  UseMethod("quick.reg.table",my.model)
}

quick.reg.table.manova = function(my.model,
                                  part.eta = F,
                                  VIF = F,
                                  myDF = my.found.df,
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
                                  do.return=T,marginality=T) {

  #### Begin MANOVA ####
  #### Inits ####
  my.y.levels=dim(my.model$model[[1]])[2]
  my.envir=environment()
  SS.type = 2
  my.nested.table=quick.SSCP(my.model, myDF, marginality, show.contrasts, show.latent,my.envir)
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
  if(!marginality){
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
    if(!marginality==3){
      latent.part.resid.total=latent.treat.model[[1]][[2]]
    }else{
      latent.part.resid.total=latent.treat.model[[1]][[2]][[1]]
    }
  }

  #### Get totals ####
  the.total=total.resid+treat.total+if(!marginality){sum(diag(treat.model[[1]][[1]][[1]]))}else{0}
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
    my.test.stat=ifelse({i==1 & marginality},NA,quick.m.test(my.treat.err,test.stat))
    my.SS=ifelse({i==1 & marginality},0,sum(diag(treat.model[[1]][[1]][[i]])))
    my.df=my.y.levels*ifelse(my.i!=1,as.numeric(my.nested.table[my.i,7]),1)
    my.resid.df=ifelse(my.i==1,{my.y.levels*the.resid.df},min({the.resid.df*as.numeric(my.nested.table[my.i,5])-as.numeric(my.nested.table[my.i,5])},my.y.levels*the.resid.df-as.numeric(my.nested.table[my.i,5])))
    my.f.val=ifelse({i==1 & marginality},NA,{my.SS/my.df}/{the.resid/the.resid.df})
    my.p.val=ifelse({i==1 & marginality},NA,pf(my.f.val,my.df,my.resid.df,lower.tail = F))

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
  my.html.table=quick.table(my.manova.table,"manova",test=test.stat,marginality=marginality, abbrev.length = abbrev.length,the.footer = the.footer)
  if(do.return){
    return(my.html.table)
  }else{
    return()
  }
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
}

quick.reg.table.lm = function(my.model,
                              part.eta = F,
                              VIF = F,
                              myDF = my.found.df,
                              marginality=T,
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

  SS.type = 2
  my.envir=environment()

  my.nested.table=quick.SSCP(my.model, myDF, marginality, show.contrasts, show.latent,my.envir)
  #### ANOVA TABLES ####
  my.summary = summary(my.model)
  my.coefficients = my.summary$coefficients
  my.coefficients = as.data.frame(my.coefficients)

    if(!marginality){
    my.III.summary = car::Anova(my.model, type = 3)
    }else{
      my.III.summary = car::Anova(my.model, type = 2)
    }

    if (type == "glm" & is.null(my.factor)) {
      the.length = dim(my.III.summary)[1] + 1
    } else{
      the.length = dim(my.III.summary)[1]
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


  my.tables.df = as.data.frame(matrix(ncol = v.p.len, nrow = 1))

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
      factor.rownames = c(factor.rownames, grep(
        paste("^", my.factor[i], "$", sep = ""),
        rownames(my.III.summary)
      ))
      num.of.levels = c(num.of.levels, length(levels(myDF[[factor.stupid[i]]])))
    }

  } else{
    factor.rownames = 0L

  }



  #### Make phia stuff ####
  if(show.contrasts){
    if(type=="lm" & !is.null(my.factor)){
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

    dang.length = length(rownames(my.III.summary))




    while (this.shift.temp < dang.length) {
      if (is.na(factor.rownames[my.factor.var])) {
        my.factor.rownames = 1

      } else{
        my.factor.rownames = factor.rownames[my.factor.var]

      }
      #### LOOP ####
      if (this.shift.temp == 1) {
        i = 1
        while (i <= total.intercepts) {
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

      }

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


      #### Make custom glance stats ####
      if(do.glance){
        #### Can eventually make it options
        glance_stats=as.data.frame(matrix(ncol={7+v.p.rep},nrow=1))
        glance_stats[1,]=c(paste("Method: ","QR Decomposition",if(show.contrasts){paste("<br />Adjustment Method: ",adjustment,sep="")},sep=""),rep(NA,6),rep(NA,v.p.rep))
      }
      #### For total
      the.length = the.length + 1
      this.temp.var = this.temp.var + 1
      #### Make table ####

      options(pixie_interactive = pix.int,
              pixie_na_string = "")


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

      if(do.glance){
        my.dust = pixiedust::redust(my.dust, glance_stats, part = "foot") %>%
          sprinkle_na_string(part = "foot") %>%
          sprinkle(rows=1,merge=T,halign="center",part="foot")
      }
      if (pix.int) {
        return(my.dust)
      } else{
        my.dust.print = print(my.dust, quote = F)[1]
        return(my.dust.print)
      }

    }
}
