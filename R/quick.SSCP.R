quick.SSCP=function(my.model, myDF, SS.type, show.contrasts, my.envir, ...){
  UseMethod("quick.SSCP", my.model)
}
quick.SSCP.lm=function(my.model,myDF,SS.type,show.contrasts,my.envir){
  my.formula.lists=quick.formula(my.model,my.envir)

  #### New DF and Number (levels) of dependent
  my.new.df=my.model$model
  my.y.levels=1
}
quick.SSCP.manova=function(my.model, myDF, SS.type, show.contrasts, show.latent, my.envir){
  #### Split and get formulas ####

  my.formula.lists=quick.tasks::quick.formula(my.model,my.envir)



  #### New DF and Number (levels) of dependent
  my.new.df=my.model$model
  my.y.levels=dim(my.model$model[[1]])[2]


  #### Run all models ####
  #### Only running Type II
  #### Since it is now lapply, can make it parallel easily later
  attach(my.new.df)
  my.null.model=manova(my.formula.lists[[1]])
  if(SS.type==3){
  my.null.model.SSCP=c(car::Anova(my.null.model,type=3)$SSP,car::Anova(my.null.model,type=3)$SSPE)
  }else{
  SSCP.temp=summary(my.null.model,intercept = T)$SS
  my.null.model.SSCP=c(SSCP.temp[-length(SSCP.temp)],SSCP.temp[length(SSCP.temp)])
  }
  my.pre.models=lapply(my.formula.lists[[2]],manova)
  my.pre.models.SSCP=lapply(my.pre.models,quick.SSP.matrix,SS.type)
  my.full.models=lapply(my.formula.lists[[3]],manova)
  my.full.models.SSCP=lapply(my.full.models,quick.SSP.matrix,SS.type)
  detach(my.new.df)


  #### Make nested data frame ####
  ## variable (str) | formula (str) | model
  my.nested.table=NULL
  my.nested.table=list("null",my.formula.lists[[1]],my.null.model,my.null.model.SSCP,NA,NA,NA)

  model.names=names(my.model$model)[-1]

  for(v in 1:length(model.names)){
    my.nested.table=rbind(my.nested.table,list(model.names[v],my.formula.lists[[2]][[v]],
                                               my.pre.models[[v]],my.pre.models.SSCP[[v]],
                                               {my.null.model$df.residual-my.pre.models[[v]]$df.residual},
                                               NA,NA))
    my.nested.table=rbind(my.nested.table,list(model.names[v],my.formula.lists[[3]][[v]],
                                               my.full.models[[v]],my.full.models.SSCP[[v]],
                                               {my.null.model$df.residual-my.full.models[[v]]$df.residual},
                                               my.full.models.SSCP[[v]][[1]][{length(my.full.models.SSCP[[v]][[1]])}],
                                               {my.null.model$df.residual-my.full.models[[v]]$df.residual}-{my.null.model$df.residual-my.pre.models[[v]]$df.residual}))

  }

  colnames(my.nested.table)=c("Variable","Formula","Model","SSCP","df","Change","df Change")

  #### Create latent variables ####
  if(show.latent){


    latent.SSCP=lapply(my.nested.table[-1,4],quick.latent,SS.type)
    latent.change=NA
    for(S in 2:{length(latent.SSCP)+1}){
      if(S %% 2 == 0){
        latent.change=c(latent.change,NA)
      }else{
        latent.change=c(latent.change,as.list(latent.SSCP[[S-1]][[1]][{length(latent.SSCP[[S-1]][[1]])}]))
      }
    }

    #### Add to table
    my.nested.table=cbind(my.nested.table,c(NA,latent.SSCP),latent.change)

    colnames(my.nested.table)=c("Variable","Formula","Model","SSCP","df","Change","df Change","Latent SSCP","Latent Change")
  }

  #### Get Contrasts ####
  if(show.contrasts){
    my.nested.table=quick.part.cont(my.nested.table,latent.cont=ifelse(show.latent,T,F))
  }

  return(my.nested.table)
}
