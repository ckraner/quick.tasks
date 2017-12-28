#### Contrasts function ####
##### Things needed #####
#my.y.levels - can get
#my.new.model - NEED - FROM THE TABLE
#my.SSP.treat.df - NEED - FROM THE LINE
#SS.type - NEED
#the.resid.SS
#the.resid.df
#my.latent.SSP.err - NEED - FROM THE LINE
#my.new.df

quick.part.cont=function(my.nested.table,SS.type,latent.cont=F){
  #### Get y levels
  my.y.levels=dim(my.nested.table[1,4][[1]][[1]])[1]

  #### Put in NA rows
  temp.colnames=colnames(my.nested.table)
  my.nested.table2=cbind(my.nested.table,as.matrix(c(rep(NA,dim(my.nested.table)[1]))))
  if(latent.cont){
    my.nested.table2=cbind(my.nested.table2,as.matrix(c(rep(NA,dim(my.nested.table)[1]))))
    colnames(my.nested.table2)=c(temp.colnames,"Contrasts","Latent Contrasts")
  }else{
    colnames(my.nested.table2)=c(temp.colnames,"Contrasts")
  }
  #### For each row besides the null model
  for(p in 2:dim(my.nested.table)[1]){
    #### Check if should have contrast
    if(!is.na(my.nested.table[p,7][[1]]) & my.nested.table[p,7][[1]]>1){
      #### If should, make contrast
      my.new.df=my.nested.table[p,3][[1]]$model
      my.MSE=mean(my.nested.table[p,3][[1]]$residuals^2)
      the.resid.SS=sum(diag(my.nested.table[p,4][[1]][[2]]))
      the.resid.df=my.nested.table[p,3][[1]]$df.residual

      if(latent.cont){
        #### Latent Mean Square Error
        my.latent.MSE=NULL
        for(j in 1:my.y.levels){
          if(i==1){
            my.latent.MSE=as.numeric(mean(my.contrast.model$residuals[i]^2))
          }else{
            my.latent.MSE=c(my.latent.MSE,as.numeric(mean(my.contrast.model$residuals[i]^2)))
          }
        }
      }

      #### Get mean responses for variables longer than 2
      #### WARN! LEVEL NAMES IS GOING TO BREAK IT!!!! ####
      level.name.grep=grep(my.nested.table[p,1],names(my.nested.table[p,3][[1]]$xlevels))
      level.names=as.vector(my.nested.table[p,3][[1]]$xlevels[[level.name.grep]])
      num.of.contrasts=my.nested.table[p,7][[1]]

      count.grep=grep(my.nested.table[p,1],names(my.new.df))
      my.count.means=NULL
      my.latent.count.means=NULL
      my.count.n=NULL
      for(j in 1:{num.of.contrasts+1}){
        if(j==1){
          my.count.means=as.vector(mean(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1]))
          my.count.n=as.vector(dim(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1])[1])
          if(latent.cont){
            for(E in 1:{my.y.levels}){
              my.latent.count.means[[E]]=as.list(mean(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1][,E]))
            }
          }
        }else{
          my.count.means=c(my.count.means,mean(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1]))
          my.count.n=c(my.count.n,dim(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1])[1])
          if(latent.cont){
            for(E in 1:{my.y.levels}){
              my.latent.count.means[[E]]=c(my.latent.count.means[[E]],mean(my.new.df[which(my.new.df[[count.grep]]==level.names[j]),1][,E]))
            }
          }
        }
      }


      #### Make contrasts
      my.contrasts=NULL
      for(j in 1:{num.of.contrasts}){
        if(j==1){
          my.contrasts=c(1,-1,rep(0,num.of.contrasts-1))
        }else{
          my.contrasts=rbind(my.contrasts,c(1,rep(0,j-1),-1,rep(0,num.of.contrasts-j)))
        }
      }


      #### Compute F values & SS
      my.contrasts.F=NULL
      my.contrasts.SSC=NULL
      if(latent.cont){
        my.latent.contrasts.F=NULL
        my.latent.contrasts.SSC=NULL
      }
      for(j in 1:{num.of.contrasts}){
        if(j==1){
          my.contrasts.I=as.integer(t(as.matrix(my.count.means)))%*%as.integer(as.matrix(my.contrasts[j,]))
          my.contrasts.denom=sum(my.contrasts[j,]^2/as.integer(my.count.n))
          my.contrasts.SSC=as.list({my.contrasts.I^2}/{my.contrasts.denom})
          my.contrasts.F=as.list({my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          if(latent.cont){
            for(l in 1:{my.y.levels}){
              my.latent.contrasts.I=as.integer(t(as.matrix(my.latent.count.means[[l]])))%*%as.integer(as.matrix(my.contrasts[j,]))
              my.latent.contrasts.SSC[[l]]=as.list({my.latent.contrasts.I^2}/{my.contrasts.denom})
              my.latent.contrasts.F[[l]]=as.list({my.latent.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
            }
          }
        }else{
          my.contrasts.I=as.integer(t(as.matrix(my.count.means)))%*%as.integer(as.matrix(my.contrasts[j,]))
          my.contrasts.denom=sum(my.contrasts[j,]^2/as.integer(my.count.n))
          my.contrasts.SSC=c(my.contrasts.SSC,{my.contrasts.I^2}/{my.contrasts.denom})
          my.contrasts.F=c(my.contrasts.F,{my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          if(latent.cont){
            for(l in 1:{my.y.levels}){
              my.latent.contrasts.I=as.integer(t(as.matrix(my.latent.count.means[[l]])))%*%as.integer(as.matrix(my.contrasts[j,]))
              my.latent.contrasts.SSC[[l]]=c(my.latent.contrasts.SSC[[l]],{my.latent.contrasts.I^2}/{my.contrasts.denom})
              my.latent.contrasts.F[[l]]=c(my.latent.contrasts.F[[l]],{my.latent.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
            }
          }
        }
      }


      #### Make rownames
      my.contrasts.names=NULL
      for(j in 1:{num.of.contrasts}){
        if(j==1){
          my.contrasts.names=as.list(paste(trimws(level.names[[1]]),"-",trimws(level.names[[j+1]])))
        }else{
          my.contrasts.names=c(my.contrasts.names,paste(trimws(level.names[[1]]),"-",trimws(level.names[[j+1]])))
        }
      }


      #### Add to table
      my.contrasts.4.table=cbind(as.matrix(unlist(my.contrasts.names)),as.numeric(as.matrix(unlist(my.contrasts.F))),as.matrix(as.numeric(unlist(my.contrasts.SSC))))
      contr.grep=grep("^Contrasts$",colnames(my.nested.table2))
      my.nested.table2[p,contr.grep]=list(my.contrasts.4.table)

      if(latent.cont){
        my.latent.contrasts.4.table=NULL
        for(V in 1:my.y.levels){
          if(V==1){
          my.latent.contrasts.4.table=cbind(as.matrix(unlist(my.contrasts.names)),as.numeric(as.matrix(unlist(my.latent.contrasts.F[[V]]))),as.matrix(as.numeric(unlist(my.latent.contrasts.SSC[[V]]))))
          }else{
            my.latent.contrasts.4.table=cbind(my.latent.contrasts.4.table,as.matrix(unlist(my.contrasts.names)),as.numeric(as.matrix(unlist(my.latent.contrasts.F[[V]]))),as.matrix(as.numeric(unlist(my.latent.contrasts.SSC[[V]]))))
          }
        }
        latent.contr.grep=grep("^Latent Contrasts$",colnames(my.nested.table2))
        my.nested.table2[p,latent.contr.grep]=list(my.latent.contrasts.4.table)

      }
    }
  }
  return(my.nested.table2)
}





#### MAKE REGULAR AND LATENT CONTRASTS ####
for(i in 1:{my.y.levels}){
  #### Pick right model
  if(i==1){
    my.contrast.model=my.new.model[[my.y.levels^2+1]]
  }else{
    my.contrast.model=my.new.model[[{i-1}*my.y.levels+1]]
  }

  #### my.model
  my.model=my.new.model[[{i-1}*my.y.levels+1]]

  my.SSP.treat=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSP

  my.SSP.err=car::Anova(my.model[[my.y.levels+1]],type=SS.type)$SSPE
  the.resid.SS=sum(diag(my.SSP.err))
  the.resid.df=my.model$df.residual

  #### Mean Square Error
  my.MSE=mean(my.contrast.model$residuals^2)


  #### Latent Mean Square Error
  my.latent.MSE=NULL
  for(j in 1:my.y.levels){
    if(i==1){
      my.latent.MSE=as.numeric(mean(my.contrast.model$residuals[i]^2))
    }else{
      my.latent.MSE=c(my.latent.MSE,as.numeric(mean(my.contrast.model$residuals[i]^2)))
    }
  }



  #### Get mean responses for variables longer than 2
  #### WARN! LEVEL NAMES IS GOING TO BREAK IT!!!! ####
  my.count.means=NULL
  my.count.n=NULL
  level.names=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      if(q==2){
        level.names=as.vector(levels(my.contrast.model[[q]]))
      }else{
        level.names=rbind(level.names,as.vector(levels(my.contrast.model[[q]])))
      }
      for(j in 1:{my.SSP.treat.df[q]+1})
        if(j==1){
          my.count.means[[q-1]]=as.list(mean(my.new.df[as.character(my.contrast.model[[q]])==level.names[j],1]))
          my.count.n[[q-1]]=as.list(dim(my.new.df[as.character(my.contrast.model[[q]])==level.names[j],1])[1])
        }else{
          my.count.means[[q-1]]=c(my.count.means[[q-1]],mean(my.new.df[as.character(my.contrast.model[[q]])==level.names[j],1]))
          my.count.n[[q-1]]=c(my.count.n[[q-1]],dim(my.new.df[as.character(my.contrast.model[[q]])==level.names[j],1])[1])
        }
    }
  }

  #### Make contrasts
  my.contrasts=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      level.names=levels(my.contrast.model[[q]])
      for(j in 1:{my.SSP.treat.df[q]})
        if(j==1){
          my.contrasts[[q-1]]=c(1,-1,rep(0,my.SSP.treat.df[q]-1))
        }else{
          my.contrasts[[q-1]]=rbind(my.contrasts[[q-1]],c(1,rep(0,j-1),-1,rep(0,my.SSP.treat.df[q]-j)))
        }
    }
  }

  #### Compute F values
  my.contrasts.F=NULL
  my.latent.contrasts.F=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      for(j in 1:{my.SSP.treat.df[q]})
        if(j==1){
          my.contrasts.I=as.integer(t(as.matrix(my.count.means[[q-1]])))%*%as.integer(as.matrix(my.contrasts[[q-1]][j,]))
          my.contrasts.denom=sum(my.contrasts[[q-1]][j,]^2/as.integer(my.count.n[[q-1]]))
          my.contrasts.F[[q-1]]=as.list({my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          for(l in 1:{my.y.levels}){
            if(l==1){
              my.latent.contrasts.F[[q-1]]=as.list({my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
            }else{
              my.latent.contrasts.F[[q-1]]=c(my.latent.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
            }
          }
        }else{
          my.contrasts.I=as.integer(t(as.matrix(my.count.means[[q-1]])))%*%as.integer(as.matrix(my.contrasts[[q-1]][j,]))
          my.contrasts.denom=sum(my.contrasts[[q-1]][j,]^2/as.integer(my.count.n[[q-1]]))
          my.contrasts.F[[q-1]]=c(my.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          for(l in 1:{my.y.levels}){
            my.latent.contrasts.F[[q-1]]=c(my.latent.contrasts.F[[q-1]],{my.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
          }
        }
    }
  }

  #### Compute SS from F values
  #### NEED TO FIX
  #### Fixed
  #### change to partial resid.SS

  my.contrasts.SS=NULL
  my.latent.contrasts.SS=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      for(j in 1:{my.SSP.treat.df[q]}){
        if(j==1){
          my.contrasts.SS[[q-1]]=as.list({as.numeric(my.contrasts.F[[q-1]][j])*{the.resid.SS}*my.y.levels}/{the.resid.df})
          for(l in 1:{my.y.levels}){
            if(l==1){
              my.latent.contrasts.SS[[q-1]]=as.list({as.numeric(my.latent.contrasts.F[[q-1]][1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
            }else{
              my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][my.y.levels*l-1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
            }
          }
        }else{
          my.contrasts.SS[[q-1]]=c(my.contrasts.SS[[q-1]],{as.numeric(my.contrasts.F[[q-1]][j])*{the.resid.SS}*my.y.levels}/{the.resid.df})
          for(l in 1:{my.y.levels}){
            if(l==1){
              my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][j])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
            }else{
              my.latent.contrasts.SS[[q-1]]=c(my.latent.contrasts.SS[[q-1]],{as.numeric(my.latent.contrasts.F[[q-1]][my.y.levels*j-1])*{my.SSP.err[q-1,q-1]}*1}/{the.resid.df-my.y.levels+1})
            }
          }
        }
      }
    }
  }


  #### Make rownames
  my.contrasts.names=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      for(j in 1:{my.SSP.treat.df[q]}){
        if(j==1){
          my.contrasts.names[[q-1]]=as.list(paste(level.names[[1]],"-",level.names[[j+1]]))
        }else{
          my.contrasts.names[[q-1]]=c(my.contrasts.names[[q-1]],paste(level.names[[1]],"-",level.names[[j+1]]))
        }
      }
    }
  }


  #### Make table ####
  #### NEED TO ADD PVAL SO CAN ADJUST ####
  #### p.adjust works on set of p-vals!
  my.contrasts.table=NULL
  for(q in 2:length(my.SSP.treat.df)){
    if(my.SSP.treat.df[q]>1){
      num.of.contrasts=my.SSP.treat.df[q]
      my.contrasts.table[[q-1]]=cbind(my.contrasts.names[[1]],my.contrasts.F[[1]],my.contrasts.SS[[1]])
      colnames(my.contrasts.table[[q-1]])=c("name","F.val","SS")
    }
  }
}


#### Latent tables ####
my.latent.contrasts.F.R=NULL
my.latent.contrasts.SS.R=NULL
for(q in 2:length(my.SSP.treat.df)){
  if(my.SSP.treat.df[q]>1){
    if(q==2){
      the.latent.levels=1
      num.of.contrasts=my.SSP.treat.df[q]
      for(s in 2:my.SSP.treat.df[q]){
        the.latent.levels=c(the.latent.levels,my.y.levels*s-1)}
      my.latent.contrasts.F.R[[q-1]]=as.list(my.latent.contrasts.F[[q-1]][[1]])
      my.latent.contrasts.SS.R[[q-1]]=as.list(my.latent.contrasts.SS[[q-1]][[1]])
      for(r in 2:length(the.latent.levels)){
        my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
        my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
      }
      for(v in 2:my.y.levels){
        the.latent.levels=NULL
        for(s in 1:my.SSP.treat.df[q]){
          the.latent.levels=c(the.latent.levels,my.y.levels*s)}
        for(r in 1:length(the.latent.levels)){
          my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
          my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
        }
      }
    }else{
      for(v in 1:my.y.levels){
        the.latent.levels=NULL
        #the.latent.levels=q-1
        for(s in 1:my.SSP.treat.df[q]){
          the.latent.levels=c(the.latent.levels,my.y.levels*s)}
        for(r in 1:length(the.latent.levels)){
          my.latent.contrasts.F.R[[q-1]]=c(my.latent.contrasts.F.R[[q-1]],my.latent.contrasts.F[[q-1]][[the.latent.levels[r]]])
          my.latent.contrasts.SS.R[[q-1]]=c(my.latent.contrasts.SS.R[[q-1]],my.latent.contrasts.SS[[q-1]][[the.latent.levels[r]]])
        }
      }
    }
  }
}
my.latent.contrasts.table=NULL
for(q in 2:length(my.SSP.treat.df)){
  num.of.contrasts=my.SSP.treat.df[q]
  if(my.SSP.treat.df[q]>1){
    for(g in 1:my.y.levels){
      my.latent.contrasts.table[[q-1]][[g]]=cbind(my.contrasts.names[[1]],my.latent.contrasts.F.R[[q-1]][{g+ifelse(g>1,{g-1}*my.SSP.treat.df[q]-1,0)}:{g*my.SSP.treat.df[q]}],my.latent.contrasts.SS.R[[q-1]][{g+ifelse(g>1,{g-1}*my.SSP.treat.df[q]-1,0)}:{g*my.SSP.treat.df[q]}])
    }
  }
}
}





























#### ORDINAL #####
#### Analysis of Devience table
#### Based on ordinal package primer by
#### Rune Hauabo B Christensen


#### Check for null model

#### Pairwise deletion of variables in data frame performed

new.df=my.model$model
new.model=update(my.model,data=new.df)
null.model=update(new.model,~1)

if(type=="ord"){
  vars.df=new.model$edf-null.model$edf
}else{
  vars.df=new.model$df.null-new.model$df.residual
}
if(vars.df==0){
  stop("This is the null model.")
}

if(type=="ord"){
  treat.dev=abs(as.numeric(levels(new.model$info$AIC)[1])-as.numeric(levels(null.model$info$AIC)[1]))
  vars.dev=drop1(new.model, test="Chi")$LRT[-1]
  vars.dev.df=drop1(new.model,test="Chi")$Df[-1]
  vars.dev.p=drop1(new.model,test="Chi")$`Pr(>Chi)`[-1]
  #total.dev=-2*{as.integer(levels(null.model$info$logLik)[1])-as.integer(levels(new.model$info$logLik)[1])}
  total.dev=anova(null.model,new.model)$LR.stat[2]
  resid.dev=total.dev-treat.dev
  resid.df=null.model$edf
  total.df=new.model$edf

  vars.or=exp(coef(new.model))[{null.model$edf+1}:length(coef(new.model))]
  vars.or.confint=exp(confint(new.model,type="Wald"))[{null.model$edf+1}:length(coef(new.model)),]
}else{
  resid.dev=new.model$deviance
  resid.df=new.model$df.residual
  vars.dev=anova(new.model,test="Chi")$Deviance[-1]
  vars.dev.df=anova(new.model,test="Chi")$Df[-1]
  vars.dev.p=anova(new.model,test="Chi")$`Pr(>Chi)`[-1]
  total.dev=new.model$null.deviance-new.model$deviance
  vars.dev.total=sum(vars.dev)

  vars.or=exp(coef(new.model))[-1]
  vars.or.confint=exp(confint(new.model,type="Wald"))[-1,]
}
my.names=names(vars.or)


dev.table=as.data.frame(matrix(ncol=7,nrow=1))
names(dev.table)=c("var","p.odd","p.odd.2.5","p.odd.97.5","deviance","df","p-val")
if(type=="ord"){
  dev.table[1,]=c("Treatment",NA,NA,NA,treat.dev,vars.df,dchisq(treat.dev,vars.df))

  dev.var2=2
}else{
  dev.var2=1
}
dev.var=1
while(dev.var < {vars.df+1}){
  dev.table[dev.var2,]=c(my.names[dev.var],vars.or[dev.var],vars.or.confint[dev.var,1],vars.or.confint[dev.var,2],vars.dev[dev.var],vars.dev.df[dev.var],vars.dev.p[dev.var])
  dev.var2=dev.var2+1
  dev.var=dev.var+1
}

dev.table[dev.var2+1,]=c("Residual",NA,NA,NA,resid.dev,resid.df,dchisq(resid.dev,resid.df))
dev.table[dev.var2+2,]=c("Total",NA,NA,NA,total.dev,total.df,dchisq(total.dev,total.df))

for(i in 2:7){
  dev.table[,i]=as.numeric(dev.table[,i])
}

options(pixie_interactive = pix.int,
        pixie_na_string = "")

dev.table.dusted=pixiedust::dust(dev.table)%>%
  sprinkle(cols = "p-val", fn = quote(pvalString(
    value, digits = 3, format = "default"
  ))) %>%
  sprinkle_print_method(pix.method) %>%
  sprinkle_na_string()%>%
  sprinkle_border(rows=1:dim(dev.table)[1],border="all")%>%
  sprinkle_round(cols=2,round=2)%>%
  sprinkle_round(cols=5:6,round=2)%>%
  sprinkle_round(cols=3:4,round=4)%>%
  sprinkle_colnames("Variable","Odds <br /> Ratio","Confidence <br /> 2.5%","Interval <br /> 97.5%","Deviance","dF","P(>Chisq)")%>%
  sprinkle_pad(cols=1:7,pad=5,part="head")%>%
  sprinkle_border(cols=1:2,border="all",part="head")%>%
  sprinkle_border(cols=5:7,border="all",part="head")%>%
  sprinkle_border(cols=3:4,border=c("top","bottom"),part="head")%>%
  sprinkle_align(rows=1,halign="center",part="head")

family.glance=as.data.frame(matrix(ncol=7,nrow=1))
family.glance[1,]=c(paste("Family: Ordinal <br />","Link method: ",levels(new.model$info$link),sep=""),rep(NA,6))
#family.glance[2,]=c(paste("Link method: ",levels(new.model$info$link),sep=""),rep(NA,6))

dev.table.dusted=redust(dev.table.dusted,family.glance,part="foot")%>%
  sprinkle_na_string()%>%
  sprinkle(rows=1,merge=T,halign="center",part="foot")

if (pix.int) {
  return(dev.table.dusted)
} else{
  my.manova.pixie = print(dev.table.dusted, quote = F)[1]
  return(my.manova.pixie)
}
