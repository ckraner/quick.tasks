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








#### Contrasts function ####
##### Things needed #####
#my.y.levels - can get

my.new.model - NEED
my.SSP.treat.df - NEED

SS.type - NEED
#the.resid.SS
#the.resid.df
my.latent.SSP.err - NEED
#my.new.df

quick.cont=function(my.model,my.ssp)
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
