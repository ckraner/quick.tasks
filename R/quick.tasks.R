###################### QUICK.TABLES ########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################## 12/2017 ###########################
############################################################

#' Partial contrasts for ANOVA and MANOVA tables
#'
#' Adds Contrasts and Latent Contrasts to nested table
#' for creating overall MANOVA tables.
#'
#' This package takes the SSCP matrices from [quick.reg()] and
#' calculates partial contrasts directly. Also can create
#' latent contrasts.
#'
#' @param my.nested.table Table of models and SSCP matrices from [quick.reg()]
#' @param adjustment P-value adjustment sent to p.adjust (default = "bonferroni").
#' @param latent.cont Should latent contrasts also be computed? (default = F).
#' @param abbrev.length Integer telling how long of a label is too long. Longer than this and will be condensed (default=15)
#' @return Nested table with added rows and information (Unfortunately as character at the moment)
#' @keywords Explore
#' @examples
#' quick.contrast(my.nested.table)

quick.part.cont=function(my.nested.table,
                         adjustment="bonferroni",
                         latent.cont=F,
                         abbrev.length=30){
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
        the.resid.SS=sum(diag(my.nested.table[p,4][[1]][[2]][[1]]))
      the.resid.df=my.nested.table[p,3][[1]]$df.residual

      if(latent.cont){
        #### Latent Mean Square Error
        my.latent.MSE=NULL
        for(j in 1:my.y.levels){
          if(j==1){
            my.latent.MSE=as.numeric(mean(my.nested.table[p,3][[1]]$residuals[j,]^2))
          }else{
            my.latent.MSE=c(my.latent.MSE,as.numeric(mean(my.nested.table[p,3][[1]]$residuals[j,]^2)))
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

      #pf(my.f.val,my.df,my.resid.df,lower.tail = F)
      #### Compute F values & SS & P-val
      my.contrasts.F=NULL
      my.contrasts.SSC=NULL
      my.contrasts.p=NULL
      if(latent.cont){
        my.latent.contrasts.F=NULL
        my.latent.contrasts.SSC=NULL
        my.latent.contrasts.p=NULL
      }
      for(j in 1:{num.of.contrasts}){
        if(j==1){
          my.contrasts.I=as.integer(t(as.matrix(my.count.means)))%*%as.integer(as.matrix(my.contrasts[j,]))
          my.contrasts.denom=sum(my.contrasts[j,]^2/as.integer(my.count.n))
          my.contrasts.SSC=as.list({my.contrasts.I^2}/{my.contrasts.denom})
          my.contrasts.F=as.list({my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          my.contrasts.p=as.list(pf(my.contrasts.F[[j]],1,the.resid.df,lower.tail = F))
          if(latent.cont){
            for(l in 1:{my.y.levels}){
              my.latent.contrasts.I=as.integer(t(as.matrix(my.latent.count.means[[l]])))%*%as.integer(as.matrix(my.contrasts[j,]))
              my.latent.contrasts.SSC[[l]]=as.list({my.latent.contrasts.I^2}/{my.contrasts.denom})
              my.latent.contrasts.F[[l]]=as.list({my.latent.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
              my.latent.contrasts.p[[l]]=as.list(pf(my.latent.contrasts.F[[l]][[j]],1,the.resid.df,lower.tail = F))
            }
          }
        }else{
          my.contrasts.I=as.integer(t(as.matrix(my.count.means)))%*%as.integer(as.matrix(my.contrasts[j,]))
          my.contrasts.denom=sum(my.contrasts[j,]^2/as.integer(my.count.n))
          my.contrasts.SSC=c(my.contrasts.SSC,{my.contrasts.I^2}/{my.contrasts.denom})
          my.contrasts.F=c(my.contrasts.F,{my.contrasts.I^2}/{my.MSE*my.contrasts.denom})
          my.contrasts.p=c(my.contrasts.p,pf(my.contrasts.F[[j]],1,the.resid.df,lower.tail = F))
          if(latent.cont){
            for(l in 1:{my.y.levels}){
              my.latent.contrasts.I=as.integer(t(as.matrix(my.latent.count.means[[l]])))%*%as.integer(as.matrix(my.contrasts[j,]))
              my.latent.contrasts.SSC[[l]]=c(my.latent.contrasts.SSC[[l]],{my.latent.contrasts.I^2}/{my.contrasts.denom})
              my.latent.contrasts.F[[l]]=c(my.latent.contrasts.F[[l]],{my.latent.contrasts.I^2}/{my.latent.MSE[l]*my.contrasts.denom})
              my.latent.contrasts.p[[l]]=c(my.latent.contrasts.p[[l]],pf(my.latent.contrasts.F[[l]][[j]],1,the.resid.df,lower.tail = F))
            }
          }
        }
      }

      #### Make p-val Adjustments
      my.contrasts.p=p.adjust(my.contrasts.p,method=adjustment)
      if(latent.cont){
        for(l in 1:{my.y.levels}){
          my.latent.contrasts.p[[l]]=p.adjust(my.latent.contrasts.p[[l]],method=adjustment)
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
      my.contrasts.4.table=cbind(as.matrix(unlist(my.contrasts.names)),as.matrix(unlist(my.contrasts.F)),as.matrix(unlist(my.contrasts.SSC)),as.matrix(my.contrasts.p))
      contr.grep=grep("^Contrasts$",colnames(my.nested.table2))
      my.nested.table2[p,contr.grep]=list(my.contrasts.4.table)

      if(latent.cont){
        my.latent.contrasts.4.table=NULL
        for(V in 1:my.y.levels){
          if(V==1){
            my.latent.contrasts.4.table=cbind(as.matrix(unlist(my.contrasts.names)),as.numeric(as.matrix(unlist(my.latent.contrasts.F[[V]]))),as.matrix(as.numeric(unlist(my.latent.contrasts.SSC[[V]]))),as.matrix(my.latent.contrasts.p[[V]]))
          }else{
            my.latent.contrasts.4.table=cbind(my.latent.contrasts.4.table,as.matrix(unlist(my.contrasts.names)),as.numeric(as.matrix(unlist(my.latent.contrasts.F[[V]]))),as.matrix(as.numeric(unlist(my.latent.contrasts.SSC[[V]]))),as.matrix(my.latent.contrasts.p[[V]]))
          }
        }
        latent.contr.grep=grep("^Latent Contrasts$",colnames(my.nested.table2))
        my.nested.table2[p,latent.contr.grep]=list(my.latent.contrasts.4.table)

      }
    }
  }
  return(my.nested.table2)
}

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

quick.contrast = function(my.model,
                          SS.type = 3,
                          adjustment = "bonferroni",
                          test.stat = "Wilks",
                          abbrev.length = 15,
                          pix.int = T,
                          pix.method = "html",
                          my.factors = my.contrasts,
                          my.type = my.reg.type,
                          skip.me=F) {
  #### Find type
  my.call = as.character(my.model$call)
  my.split.call = strsplit(my.call, "\\(")
  my.reg.type = my.split.call[[1]][1]


  #### Find factors
  my.contrasts = names(my.model$contrasts)

  if (my.type != "manova" | my.type != "stats::manova") {
    #### Find levels
    my.x.levels = NULL
    for (i in 1:length(my.factors)) {
      my.x.levels = c(my.x.levels, my.model$xlevels[[i]])
    }

    #### Find num non var
    my.non.var = length(my.model$coefficients) - length(my.x.levels) - 1 +
      length(my.contrasts)
  }

  library(pixiedust)
  library(phia)
  if (my.type == "manova" | my.type == "stats::manova") {
    my.phia.print = as.data.frame(matrix(ncol = 8, nrow = 1))
    my.lengths = NULL
    this.table.var = 1

    while (this.table.var < {
      length(my.factors) + 1
    }) {
      my.phia = phia::testInteractions(
        my.model,
        fixed = my.factors[this.table.var],
        adjustment = adjustment,
        abbrev.levels = abbrev.length
      )
      my.phia$names = attr(my.phia, "row.names")
      my.phia = my.phia[c("names",
                          "Df",
                          "test stat",
                          "approx F",
                          "num Df",
                          "den Df",
                          "Pr(>F)")]
      attr(my.phia, "class") = attr(my.phia, "class")[-1]
      my.lengths = c(my.lengths, nrow(my.phia))
      this.stuff = c(my.factors[this.table.var], NA)

      if (my.lengths[this.table.var] > 2) {
        for (i in 1:{
          my.lengths[this.table.var] - 2
        }) {
          this.stuff = c(this.stuff, NA)

        }

      }

      my.phia = cbind(this.stuff, my.phia)

      if (this.table.var == 1) {
        my.phia.print = my.phia

      } else{
        my.phia.print = rbind(my.phia.print, my.phia)

      }

      this.table.var = this.table.var + 1
      #print(my.phia.print)

    }

    rownames(my.phia.print) = NULL

    if(skip.me){
      return(my.phia.print)
    }
    phia.length = dim(my.phia.print)[1]

    options(pixie_interactive = pix.int)

    my.phia.pixie = pixiedust::dust(my.phia.print) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle(cols = "Pr(>F)", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle(cols = "test stat", round = 3) %>%
      sprinkle(cols = "approx F", round = 3) %>%
      sprinkle_colnames(
        "",
        "Levels",
        "df",
        "Pillai <br /> Statistic",
        "approx <br /> F-value",
        "num <br /> df",
        "den <br /> df",
        "Pr(>F)"
      ) %>%
      sprinkle(cols = 1:8,
               rows = {
                 sum(my.lengths)
               },
               border = "bottom") %>%
      sprinkle(cols = 1:8, pad = 10) %>%
      sprinkle(cols = 1,
               rows = 1:{
                 sum(my.lengths)
               },
               border = "left") %>%
      sprinkle(
        cols = 3:8,
        rows = 1:{
          sum(my.lengths)
        },
        border = c("right", "left")
      ) %>%
      sprinkle(
        cols = 1:8,
        rows = 1,
        border = c("top", "bottom"),
        part = "head"
      ) %>%
      sprinkle(
        cols = 1,
        rows = 1,
        border = "left",
        part = "head"
      ) %>%
      sprinkle(
        cols = 3:8,
        rows = 1,
        border = c("right", "left"),
        part = "head"
      ) %>%
      sprinkle_na_string(na_string = "") %>%
      sprinkle_width(
        cols = 1,
        rows = 1:2,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(
        cols = 2,
        rows = 1:2,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(cols = 3,
                     width = 30,
                     width_units = "pt") %>%
      sprinkle_width(cols = 4,
                     width = 60,
                     width_units = "pt") %>%
      sprinkle_width(cols = 5,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 6,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 7,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 8,
                     width = 70,
                     width_units = "pt") %>%
      sprinkle(rows = 1,
               halign = "center",
               part = "head")

    adj.method = as.data.frame(matrix(ncol = 8, nrow = 1))
    adj.method[1, ] = c(paste("p adjustment method: ", adjustment, sep =
                                ""),
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        NA)
    my.phia.pixie = redust(my.phia.pixie, adj.method, part = "foot") %>%
      sprinkle(merge = T,
               halign = "center",
               part = "foot")

    if (pix.int & !skip.me) {
      return(my.phia.pixie)
    } else if(!pix.int & !skip.me){
      my.phia.pixie = print(my.phia.pixie, quote = F)[1]
      return(my.phia.pixie)
    }else{
      return(my.phia.print)
    }
  } else if (my.type == "lm" | my.type == "stats::lm") {
    my.phia.print = as.data.frame(matrix(ncol = 7, nrow = 1))
    my.lengths = NULL
    this.table.var = 1

    while (this.table.var < {
      length(my.factors) + 1
    }) {
      my.phia = phia::testInteractions(
        my.model,
        pairwise = my.factors[this.table.var],
        adjustment = adjustment,
        abbrev.levels = abbrev.length
      )
      my.phia = my.phia[-{
        dim(my.phia)[1]
      }, ]
      my.phia$names = attr(my.phia, "row.names")
      my.phia = my.phia[c("names", "Value", "Df", "Sum of Sq", "F", "Pr(>F)")]
      attr(my.phia, "class") = attr(my.phia, "class")[-1]
      my.lengths = c(my.lengths, {
        nrow(my.phia)
      })
      this.stuff = c(my.factors[this.table.var], NA)

      if (my.lengths[this.table.var] > 2) {
        for (i in 1:{
          my.lengths[this.table.var] - 2
        }) {
          this.stuff = c(this.stuff, NA)

        }

      }

      my.phia = cbind(this.stuff, my.phia)

      if (my.lengths[this.table.var] == 1) {
        my.phia = my.phia[-2, ]
      }
      if (this.table.var == 1) {
        my.phia.print = my.phia

      } else{
        my.phia.print = rbind(my.phia.print, my.phia)

      }

      this.table.var = this.table.var + 1
      #print(my.phia.print)

    }
    #print(my.phia.print)
    #my.phia.print=my.phia.print[-{dim(my.phia.print)[1]-1},]
    rownames(my.phia.print) = NULL
    phia.length = dim(my.phia.print)[1]

    options(pixie_interactive = pix.int)

    my.phia.pixie = pixiedust::dust(my.phia.print) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle(cols = "Pr(>F)", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle(cols = "test stat", round = 3) %>%
      sprinkle(cols = "approx F", round = 3) %>%
      sprinkle_colnames("",
                        "Levels",
                        "Value",
                        "df",
                        "Sums of <br /> Squares",
                        "F-value",
                        "Pr(>F)") %>%
      sprinkle(cols = 1:7,
               rows = {
                 sum(my.lengths)
               },
               border = "bottom") %>%
      sprinkle(cols = 1:7, pad = 10) %>%
      sprinkle(cols = 1,
               rows = 1:{
                 sum(my.lengths)
               },
               border = "left") %>%
      sprinkle(
        cols = 3:7,
        rows = 1:{
          sum(my.lengths)
        },
        border = c("right", "left")
      ) %>%
      sprinkle(
        cols = 1:7,
        rows = 1,
        border = c("top", "bottom"),
        part = "head"
      ) %>%
      sprinkle(
        cols = 1,
        rows = 1,
        border = "left",
        part = "head"
      ) %>%
      sprinkle(
        cols = 3:7,
        rows = 1,
        border = c("right", "left"),
        part = "head"
      ) %>%
      sprinkle_na_string(na_string = "") %>%
      sprinkle_width(
        cols = 1,
        rows = 1,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(
        cols = 2,
        rows = 1,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(cols = 3,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 4,
                     width = 30,
                     width_units = "pt") %>%
      sprinkle_width(cols = 5,
                     width = 60,
                     width_units = "pt") %>%
      sprinkle_width(cols = 6,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 7,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle(rows = 1,
               halign = "center",
               part = "head")

    adj.method = as.data.frame(matrix(ncol = 7, nrow = 1))
    adj.method[1, ] = c(paste("p adjustment method: ", adjustment, sep =
                                ""),
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        NA)
    my.phia.pixie = redust(my.phia.pixie, adj.method, part = "foot") %>%
      sprinkle(merge = T,
               halign = "center",
               part = "foot")

    if (pix.int & !skip.me) {
      return(my.phia.pixie)
    } else if(!pix.int & !skip.me){
      my.phia.pixie = print(my.phia.pixie, quote = F)[1]
      return(my.phia.pixie)
    }else{
      return(my.phia.print)
    }
  } else{
    my.phia.print = as.data.frame(matrix(ncol = 6, nrow = 1))
    my.lengths = NULL
    this.table.var = 1

    while (this.table.var < {
      length(my.factors) + 1
    }) {
      my.phia = phia::testInteractions(
        my.model,
        pairwise = my.factors[this.table.var],
        adjustment = adjustment,
        abbrev.levels = abbrev.length
      )
      my.phia = my.phia[-{
        dim(my.phia)[1]
      }, ]
      my.phia$names = attr(my.phia, "row.names")
      my.phia = my.phia[c("names", "Value", "Df", "Chisq", "Pr(>Chisq)")]
      attr(my.phia, "class") = attr(my.phia, "class")[-1]
      my.lengths = c(my.lengths, {
        nrow(my.phia)
      })
      this.stuff = c(my.factors[this.table.var], NA)

      if (my.lengths[this.table.var] > 2) {
        for (i in 1:{
          my.lengths[this.table.var] - 2
        }) {
          this.stuff = c(this.stuff, NA)

        }

      }

      my.phia = cbind(this.stuff, my.phia)

      if (my.lengths[this.table.var] == 1) {
        my.phia = my.phia[-2, ]
      }
      if (this.table.var == 1) {
        my.phia.print = my.phia

      } else{
        my.phia.print = rbind(my.phia.print, my.phia)

      }

      this.table.var = this.table.var + 1
      #print(my.phia.print)

    }
    #print(my.phia.print)
    #my.phia.print=my.phia.print[-{dim(my.phia.print)[1]-1},]
    rownames(my.phia.print) = NULL
    phia.length = dim(my.phia.print)[1]
    my.phia.print[[6]] = as.numeric(my.phia.print[[6]])


    options(pixie_interactive = pix.int)
    my.phia.pixie = pixiedust::dust(my.phia.print) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle(cols = "Pr(>Chisq)", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle(cols = "value", round = 3) %>%
      sprinkle(cols = "Chisq", round = 3) %>%
      sprinkle_colnames("", "Levels", "Value", "df", "Chi-Sq", "Pr(>Chisq)") %>%
      sprinkle(cols = 1:6,
               rows = {
                 sum(my.lengths)
               },
               border = "bottom") %>%
      sprinkle(cols = 1:6, pad = 10) %>%
      sprinkle(cols = 1,
               rows = 1:{
                 sum(my.lengths)
               },
               border = "left") %>%
      sprinkle(
        cols = 3:6,
        rows = 1:{
          sum(my.lengths)
        },
        border = c("right", "left")
      ) %>%
      sprinkle(
        cols = 1:6,
        rows = 1,
        border = c("top", "bottom"),
        part = "head"
      ) %>%
      sprinkle(
        cols = 1,
        rows = 1,
        border = "left",
        part = "head"
      ) %>%
      sprinkle(
        cols = 3:6,
        rows = 1,
        border = c("right", "left"),
        part = "head"
      ) %>%
      sprinkle_na_string(na_string = "") %>%
      sprinkle_width(
        cols = 1,
        rows = 1,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(
        cols = 2,
        rows = 1,
        width = 70,
        width_units = "pt"
      ) %>%
      sprinkle_width(cols = 3,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle_width(cols = 4,
                     width = 30,
                     width_units = "pt") %>%
      sprinkle_width(cols = 5,
                     width = 60,
                     width_units = "pt") %>%
      sprinkle_width(cols = 6,
                     width = 50,
                     width_units = "pt") %>%
      sprinkle(rows = 1,
               halign = "center",
               part = "head")

    adj.method = as.data.frame(matrix(ncol = 6, nrow = 1))
    adj.method[1, ] = c(paste("p adjustment method: ", adjustment, sep =
                                ""),
                        NA,
                        NA,
                        NA,
                        NA,
                        NA)
    my.phia.pixie = redust(my.phia.pixie, adj.method, part = "foot") %>%
      sprinkle(merge = T,
               halign = "center",
               part = "foot")
    if (pix.int & !skip.me) {
      return(my.phia.pixie)
    } else if(!pix.int & !skip.me){
      my.phia.pixie = print(my.phia.pixie, quote = F)[1]
      return(my.phia.pixie)
    }else{
      return(my.phia.print)
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

quick.lavaan = function(myfit) {
  require(lavaan)
  require(pixiedust)
  require(tibble)

  prev.width = getOption("width")
  options(width = 80)

  mysummary = capture.output(lavaan::summary(myfit,
                                             standardized = TRUE, rsq = T))
  #### Fit Table ####
  my.fit.table = as.data.frame(matrix(nrow = 7, ncol = 4))
  my.fit.table[1, ] = c("Number of Iterations",
                        lavInspect(myfit, what = "iterations"),
                        NA,
                        NA)
  my.fit.table[2, ] = c("Total Observations",
                        lavInspect(myfit, what = "ntotal"),
                        NA,
                        NA)
  my.fit.table[3, ] = c(
    "Chi-Sq Test of Fit",
    round(fitMeasures(myfit)[3], 2),
    fitMeasures(myfit)[4],
    pixiedust::pval_string(fitMeasures(myfit)[5])
  )
  my.fit.table[4, ] = c("Comparitive Fit Index", round(fitMeasures(myfit)[9], 3), NA, NA)
  my.fit.table[5, ] = c("Tucker-Lewis Index", round(fitMeasures(myfit)[10], 3), NA, NA)
  my.fit.table[6, ] = c("RMSEA",
                        round(fitMeasures(myfit)[23], 3),
                        NA,
                        pixiedust::pval_string(fitMeasures(myfit)[26]))
  my.fit.table[7, ] = c("SRMR", round(fitMeasures(myfit)[29], 3), NA, NA)
  colnames(my.fit.table) = c("Name", "Value", "df", "p-val")
  #my.fit.table

  #### R Sq ####
  my.r2 = lavaan::lavInspect(myfit, what = "r2")
  my.r2 = as.data.frame(my.r2)
  my.r2 = tibble::rownames_to_column(my.r2)
  my.r2$my.r2 = round(my.r2$my.r2, 3)
  #my.r2

  #### Other tables ####
  lavaan.latent = NULL
  lavaan.covar = NULL
  lavaan.vari = NULL
  lavaan.reg = NULL
  lavaan.latent.temp = 1
  lavaan.temp = 1
  while (lavaan.temp < {
    length(mysummary) + 1
  }) {
    ### Latent Variables matrix
    if (length(grep("Latent Variables:", mysummary[lavaan.temp])) > 0) {
      lavaan.temp = lavaan.temp + 2
      while ((length(grep("Covariances:", mysummary[lavaan.temp])) == 0) &&
             (length(grep("Regressions:", mysummary[lavaan.temp])) == 0) &&
             (length(grep("Variances:", mysummary[lavaan.temp])) == 0) &&
             (length(grep("R-Square:", mysummary[lavaan.temp])) == 0)) {
        str.temp = strsplit(mysummary[lavaan.temp], "  ")
        if (length(str.temp[[1]]) > 0) {
          lavaan.latent[[lavaan.latent.temp]] = list()
          for (i in 1:(length(str.temp[[1]]))) {
            if (nchar(str.temp[[1]][i]) > 0) {
              lavaan.latent[[lavaan.latent.temp]] = c(lavaan.latent[[lavaan.latent.temp]], str.temp[[1]][i])
            }
          }
          lavaan.latent.temp = lavaan.latent.temp + 1
        }
        lavaan.temp = lavaan.temp + 1
      }

      if (length(grep("Regressions:", mysummary[lavaan.temp])) != 0) {
        #### Regressions Matrix
        lavaan.temp = lavaan.temp + 2
        lavaan.latent.temp = 1
        while ((length(grep("Covariances:", mysummary[lavaan.temp])) == 0) &&
               (length(grep("Variances:", mysummary[lavaan.temp])) == 0) &&
               (length(grep("R-Square:", mysummary[lavaan.temp])) == 0)) {
          str.temp = strsplit(mysummary[lavaan.temp], "  ")
          if (length(str.temp[[1]]) > 0) {
            lavaan.reg[[lavaan.latent.temp]] = list()
            for (i in 1:(length(str.temp[[1]]))) {
              if (nchar(str.temp[[1]][i]) > 0) {
                lavaan.reg[[lavaan.latent.temp]] = c(lavaan.reg[[lavaan.latent.temp]], str.temp[[1]][i])
              }
            }
            lavaan.latent.temp = lavaan.latent.temp + 1
          }
          lavaan.temp = lavaan.temp + 1
        }
      }

      #### Covariance matrix
      if ((length(grep("Covariances:", mysummary[lavaan.temp])) != 0)) {
        lavaan.temp = lavaan.temp + 2
        lavaan.latent.temp = 1
        while ((length(grep("Variances:", mysummary[lavaan.temp])) == 0) &&
               (length(grep("R-Square:", mysummary[lavaan.temp])) == 0)) {
          str.temp = strsplit(mysummary[lavaan.temp], "  ")
          if (length(str.temp[[1]]) > 0) {
            lavaan.covar[[lavaan.latent.temp]] = list()
            for (i in 1:(length(str.temp[[1]]))) {
              if (nchar(str.temp[[1]][i]) > 0) {
                lavaan.covar[[lavaan.latent.temp]] = c(lavaan.covar[[lavaan.latent.temp]], str.temp[[1]][i])
              }
            }
            lavaan.latent.temp = lavaan.latent.temp + 1
          }
          lavaan.temp = lavaan.temp + 1
        }
      }

      #### Variance matrix
      if ((length(grep("Variances:", mysummary[lavaan.temp])) != 0)) {
        lavaan.temp = lavaan.temp + 2
        lavaan.latent.temp = 1
        while (length(grep("R-Square:", mysummary[lavaan.temp])) == 0) {
          str.temp = strsplit(mysummary[lavaan.temp], "  ")
          if (length(str.temp[[1]]) > 0) {
            lavaan.vari[[lavaan.latent.temp]] = list()
            for (i in 1:(length(str.temp[[1]]))) {
              if (nchar(str.temp[[1]][i]) > 0) {
                lavaan.vari[[lavaan.latent.temp]] = c(lavaan.vari[[lavaan.latent.temp]], str.temp[[1]][i])
              }
            }
            lavaan.latent.temp = lavaan.latent.temp + 1
          }
          lavaan.temp = lavaan.temp + 1
        }
      }

    }
    lavaan.temp = lavaan.temp + 1
  }


  #### lavaan.covar
  if (length(lavaan.covar) > 0) {
    covar.table = as.data.frame(matrix(nrow = length(lavaan.covar), ncol = 7))
    table.count = 1
    for (i in 1:length(lavaan.covar)) {
      if (length(lavaan.covar[[i]]) < 7) {
        covar.table[table.count, ] = c(lavaan.covar[[i]][1], NA, NA, NA, NA, NA, NA)
      }
      else{
        covar.table[table.count, ] = c(
          lavaan.covar[[i]][1],
          lavaan.covar[[i]][2],
          lavaan.covar[[i]][3],
          lavaan.covar[[i]][6],
          lavaan.covar[[i]][7],
          lavaan.covar[[i]][4],
          lavaan.covar[[i]][5]
        )
      }
      table.count = table.count + 1
    }
    colnames(covar.table) = c("Name",
                              "Estimate",
                              "Std.Err",
                              "Std.lv",
                              "Std.all",
                              "z-value",
                              "P(>|z|)")
    #covar.table
  } else{
    covar.table = NULL
  }

  #### lavaan.latent
  if (length(lavaan.latent) > 0) {
    latent.table = as.data.frame(matrix(nrow = length(lavaan.latent), ncol =
                                          7))
    latent.table.count = 1
    for (i in 1:length(lavaan.latent)) {
      if (length(lavaan.latent[[i]]) < 3) {
        latent.table[latent.table.count, ] = c(lavaan.latent[[i]][[1]], NA, NA, NA, NA, NA, NA)
      } else if (length(lavaan.latent[[i]]) == 4) {
        latent.table[latent.table.count, ] = c(
          lavaan.latent[[i]][[1]],
          lavaan.latent[[i]][[2]],
          NA,
          lavaan.latent[[i]][[3]],
          lavaan.latent[[i]][[4]],
          NA,
          NA
        )
      } else{
        latent.table[latent.table.count, ] = c(
          lavaan.latent[[i]][[1]],
          lavaan.latent[[i]][[2]],
          lavaan.latent[[i]][[3]],
          lavaan.latent[[i]][[6]],
          lavaan.latent[[i]][[7]],
          lavaan.latent[[i]][[4]],
          lavaan.latent[[i]][[5]]
        )
      }
      latent.table.count = latent.table.count + 1
    }
    colnames(latent.table) = c("Name",
                               "Estimate",
                               "Std.Err",
                               "Std.lv",
                               "Std.all",
                               "z-value",
                               "P(>|z|)")
    #latent.table
  } else{
    latent.table = NULL
  }

  #### lavaan.vari
  if (length(lavaan.vari) > 0) {
    vari.table = as.data.frame(matrix(nrow = length(lavaan.vari), ncol = 7))
    vari.table.count = 1
    for (i in 1:length(lavaan.vari)) {
      vari.table[vari.table.count, ] = c(
        lavaan.vari[[i]][[1]],
        lavaan.vari[[i]][[2]],
        lavaan.vari[[i]][[3]],
        lavaan.vari[[i]][[6]],
        lavaan.vari[[i]][[7]],
        lavaan.vari[[i]][[4]],
        lavaan.vari[[i]][[5]]
      )
      vari.table.count = vari.table.count + 1
    }
    colnames(vari.table) = c("Name",
                             "Estimate",
                             "Std.Err",
                             "Std.lv",
                             "Std.all",
                             "z-value",
                             "P(>|z|)")
  } else{
    vari.table = NULL
  }

  #### lavaan.reg
  if (length(lavaan.reg) > 0) {
    reg.table = as.data.frame(matrix(nrow = length(lavaan.reg), ncol = 7))
    reg.table.count = 1
    for (i in 1:length(lavaan.reg)) {
      if (length(lavaan.reg[[i]]) < 7) {
        reg.table[reg.table.count, ] = c(lavaan.reg[[i]][1], NA, NA, NA, NA, NA, NA)
      }
      else{
        reg.table[reg.table.count, ] = c(
          lavaan.reg[[i]][1],
          lavaan.reg[[i]][2],
          lavaan.reg[[i]][3],
          lavaan.reg[[i]][6],
          lavaan.reg[[i]][7],
          lavaan.reg[[i]][4],
          lavaan.reg[[i]][5]
        )
      }
      reg.table.count = reg.table.count + 1
    }
    colnames(reg.table) = c("Name",
                            "Estimate",
                            "Std.Err",
                            "Std.lv",
                            "Std.all",
                            "z-value",
                            "P(>|z|)")
  } else{
    reg.table = NULL
  }

  dusted.fit.table = pixiedust::dust(my.fit.table) %>%
    pixiedust::sprinkle_na_string(na_string = "") %>%
    pixiedust::sprinkle_print_method(print_method = "html") %>%
    pixiedust::sprinkle_colnames("", "Value", "df", "P-Val") %>%
    pixiedust::sprinkle_border(border = "all") %>%
    pixiedust::sprinkle_border(border = "all", part = "head") %>%
    pixiedust::sprinkle_pad(pad = 7) %>%
    pixiedust::sprinkle_align(halign = "center", part = "head")
  #dusted.fit.table

  latent.table$`P(>|z|)` = as.numeric(latent.table$`P(>|z|)`)
  dusted.latent.table = pixiedust::dust(latent.table) %>%
    pixiedust::sprinkle_na_string(na_string = "") %>%
    pixiedust::sprinkle_print_method(print_method = "html") %>%
    pixiedust::sprinkle_colnames("Latent",
                                 "Estimate",
                                 "Std.Err",
                                 "Std.lv",
                                 "Std.all",
                                 "z-value",
                                 "P(>|z|)") %>%
    pixiedust::sprinkle_border(border = "all") %>%
    pixiedust::sprinkle_border(border = "all", part = "head") %>%
    pixiedust::sprinkle_pad(pad = 7) %>%
    pixiedust::sprinkle_pad(pad = 7, part = "head") %>%
    pixiedust::sprinkle(cols = 7, fn = quote(pvalString(value))) %>%
    pixiedust::sprinkle_align(halign = "center", part = "head")
  #dusted.latent.table

  if (length(reg.table) > 0) {
    reg.table$`P(>|z|)` = as.numeric(reg.table$`P(>|z|)`)
    dusted.reg.table = pixiedust::dust(reg.table) %>%
      pixiedust::sprinkle_na_string(na_string = "") %>%
      pixiedust::sprinkle_print_method(print_method = "html") %>%
      pixiedust::sprinkle_colnames("Regression",
                                   "Estimate",
                                   "Std.Err",
                                   "Std.lv",
                                   "Std.all",
                                   "z-value",
                                   "P(>|z|)") %>%
      pixiedust::sprinkle_border(border = "all") %>%
      pixiedust::sprinkle_border(border = "all", part = "head") %>%
      pixiedust::sprinkle_pad(pad = 7) %>%
      pixiedust::sprinkle_pad(pad = 7, part = "head") %>%
      pixiedust::sprinkle(cols = 7, fn = quote(pvalString(value))) %>%
      pixiedust::sprinkle_align(halign = "center", part = "head")
    #dusted.reg.table
  } else{
    dusted.reg.table = NULL
  }

  if (length(covar.table) > 0) {
    covar.table$`P(>|z|)` = as.numeric(covar.table$`P(>|z|)`)
    dusted.covar.table = pixiedust::dust(covar.table) %>%
      pixiedust::sprinkle_na_string(na_string = "") %>%
      pixiedust::sprinkle_print_method(print_method = "html") %>%
      pixiedust::sprinkle_colnames("Covariance",
                                   "Estimate",
                                   "Std.Err",
                                   "Std.lv",
                                   "Std.all",
                                   "z-value",
                                   "P(>|z|)") %>%
      pixiedust::sprinkle_border(border = "all") %>%
      pixiedust::sprinkle_border(border = "all", part = "head") %>%
      pixiedust::sprinkle_pad(pad = 7) %>%
      pixiedust::sprinkle_pad(pad = 7, part = "head") %>%
      pixiedust::sprinkle(cols = 7, fn = quote(pvalString(value))) %>%
      pixiedust::sprinkle_align(halign = "center", part = "head")
    #dusted.covar.table
  } else{
    dusted.covar.table = NULL
  }

  if (length(vari.table)) {
    vari.table$`P(>|z|)` = as.numeric(vari.table$`P(>|z|)`)
    dusted.vari.table = pixiedust::dust(vari.table) %>%
      pixiedust::sprinkle_na_string(na_string = "") %>%
      pixiedust::sprinkle_print_method(print_method = "html") %>%
      pixiedust::sprinkle_colnames("Variance",
                                   "Estimate",
                                   "Std.Err",
                                   "Std.lv",
                                   "Std.all",
                                   "z-value",
                                   "P(>|z|)") %>%
      pixiedust::sprinkle_border(border = "all") %>%
      pixiedust::sprinkle_border(border = "all", part = "head") %>%
      pixiedust::sprinkle_pad(pad = 7) %>%
      pixiedust::sprinkle_pad(pad = 7, part = "head") %>%
      pixiedust::sprinkle(cols = 7, fn = quote(pvalString(value))) %>%
      pixiedust::sprinkle_align(halign = "center", part = "head")
    #dusted.vari.table
  } else{
    dusted.vari.table = NULL
  }

  dusted.r2 = pixiedust::dust(my.r2) %>%
    pixiedust::sprinkle_na_string(na_string = "") %>%
    pixiedust::sprinkle_print_method(print_method = "html") %>%
    pixiedust::sprinkle_border(border = "all") %>%
    pixiedust::sprinkle_border(border = "all", part = "head") %>%
    pixiedust::sprinkle_pad(pad = 7) %>%
    pixiedust::sprinkle_pad(pad = 7, part = "head") %>%
    pixiedust::sprinkle_align(halign = "center", part = "head") %>%
    pixiedust::sprinkle_colnames("Variable", "R^2")
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
  mydustlist = list(
    dusted.fit.table,
    dusted.latent.table,
    dusted.reg.table,
    dusted.covar.table,
    dusted.vari.table,
    dusted.r2
  )
  #print(mydustlist[[1]])
  dust.num = 1
  while (dust.num < {
    length(mydustlist) + 1
  }) {
    print(mydustlist[[dust.num]])
    x = readline(prompt = "Press z and enter to go back, or enter to go forward\n")
    if (x != "z") {
      dust.num = dust.num + 1
    } else{
      if (dust.num > 1) {
        dust.num = dust.num - 1
      }
    }
  }
  options(width = prev.width)
}



#' Survey designs for multinomial logistic functions.
#'
#' The "survey" package does not directly compute designs with multinomial functions.
#' Instead, replicate weights must be computed and applied tediously to every function.
#' This wrapper function will compute the major pieces of information using your design.
#'
#' @param my.formula Formula to be used
#' @param design Design to be used
#' @param my.df Dataframe
#' @param type Type of repetition function. see svrepdesign for more info (default: bootstrap)
#' @param replicates Number of replications to use (default: 10)
#'
#' @return List of: table of beta, z, probability; ANODE table; omnibus test; Pseudo R-2 values
#' @export
#'

quick.multinom.survey=function(my.formula,design,my.df,type="bootstrap",replicates=10){
  ### Replicates use based on https://rpubs.com/corey_sparks/58926 who cites survey package author.
  library(nnet)
  library(AER)
  library(survey)
  library(pscl)

  #### Make intercept formula
  my.formula=as.character(my.formula)
  my.base.form=paste(my.formula[2],my.formula[1],"1")
  my.formula=paste(my.formula[2],my.formula[1],my.formula[3])

  #### Make normal formulas to take df
  temp.mult=multinom(as.formula(my.formula),trace=FALSE,data=my.df)
  temp.mult.null=multinom(as.formula(my.base.form),trace=FALSE,data=my.df)
  mult.df=temp.mult$edf
  mult.null.df=temp.mult.null$edf
  mult.df.change=mult.df-mult.null.df

  #### Make "replicate weights"
  des.rep=as.svrepdesign(design,type="bootstrap",replicates=10)

  #### Get null deviance
  mdev.null=as.table(withReplicates(des.rep,substitute(deviance(multinom(eval(as.formula(my.base.form),envir = .GlobalEnv),weights=.weights,trace=F)))))

  #### Get everything from new model
  mfitcoef=as.table(withReplicates(des.rep,substitute(lmtest::coeftest(multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=.weights,trace=F)))))
  mfitcoef=as.data.frame(cbind(mfitcoef,exp(mfitcoef[,1]),{1/exp(mfitcoef[,1])}))
  names(mfitcoef)=c(names(mfitcoef)[1:4],"OR","1/OR")
  mfitcoef=mfitcoef[,c(1,2,5,6,3,4)]
  mdev=as.table(withReplicates(des.rep,substitute(deviance(multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=.weights,trace=F)))))
  mcar=as.table(withReplicates(des.rep,substitute(as.matrix(car::Anova(multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=.weights,trace=F),type=3)))))
  #mconfint=as.table(withReplicates(des.rep,substitute(confint(multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=.weights,trace=F)))))
  mpr2=as.table(withReplicates(des.rep,substitute(pscl::pR2(multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=.weights,trace=F)))))

  #### Get omnibus chi square change
  #### Still using central distribution
  dev.change=as.numeric(mdev-mdev.null)
  dev.p=pchisq(dev.change,mult.df.change)

  return(list(mfitcoef,mcar,list(as.numeric(mdev),as.numeric(mdev.null),dev.change,dev.p),mpr2))
}




#' Make ANODE table of pooled values
#'
#' Helper for pooled multinomial analyses analyses
#'
#' @param my.mult multinomial mira dataframe
#' @param task ANODE for ANODE table, coef for z-table, R2 for pR2
#' @param type Type for car::Anova (default=3)
#'
#' @return Pooled values for task
#' @export
#'

quick.multi.anode=function(my.mult,task,type=3){


  #my.mult=with(imputed.df,substitute(nnet::multinom(eval(as.formula(my.formula),envir = .GlobalEnv),weights=eval(weight))))
  if(task=="ANODE"){
  x.man=lapply(my.mult$analyses,car::Anova,type=type)
  ncols=3
  }else if(task=="coef"){
    library(AER)
    x.man=lapply(my.mult$analyses,lmtest::coeftest)
    ncols=4
  }else if(task=="R2"){
    library(pscl)
    x.man=invisible(lapply(my.mult$analyses,pscl::pR2))
    x.man=x.man
    ncols=6
    }else{
    warning("Not supported...yet!")
  }
  my.Anode=NULL
  for(j in 1:dim(x.man[[1]])[1]){
    my.val=0
    for(i in 1:length(x.man)){
      my.val=my.val+x.man[[i]][j,1:ncols]
    }
    my.val=my.val/length(x.man)
    if(is.null(my.Anode)){
      my.Anode=my.val
    }else{
      my.Anode=rbind(my.Anode,my.val)
    }

  }
  if(task=="coef"){
    rownames(my.Anode)=rownames(x.man[[1]])
  }
  return(my.Anode)
}


#' Preform weighted multnomial regression with imputed data frames
#'
#' This will perform multinomial logistic regression with the option of including weights
#' using an imputed data frame from the MICE package. Returns ANODE table, z-table, deviance,
#' partial R2 statistics, and the analyses. Can return null model as well. Can be combined with
#' quick.table.multi.model for omnibus table of multiple models.
#'
#' @param my.formula Formula to test. Cannot test null model, but can have it returned with any other model
#' @param weights Weights to use
#' @param imputed.df Imputed df from MICE
#' @param include.null Include the null model as well?
#' @param round If want values rounded. Default: T
#'
#' @return object of class "pooled", "multinom"
#' @export
#'
#' @examples
quick.pooled.multinom=function(my.formula,weights,imputed.df,include.null=F,round=T){
  library(mice)
  library(nnet)
  library(AER)
  library(pscl)

  #### Make intercept formula
  imputed.df.complete=complete(imputed.df)
  my.formula=as.character(my.formula)
  my.base.form=paste(my.formula[2],my.formula[1],"1")
  my.formula=paste(my.formula[2],my.formula[1],my.formula[3])


  my.mult1=with(imputed.df,nnet::multinom(eval(parse(text=as.character(my.base.form))),weights=weight,trace=F))
  null.dev=mean(sapply(my.mult1$analyses,deviance))
  null.df=nnet::multinom(rschl_F~1,trace=F,data=imputed.df.complete)$edf


  my.mult=with(imputed.df,nnet::multinom(eval(parse(text=my.formula)),weights=weight,trace=F))
  my.mult.anode=quick.tasks::quick.multi.anode(my.mult,task = "ANODE")
  my.mult.coef=quick.tasks::quick.multi.anode(my.mult,task = "coef")
  #quick.tasks::quick.multi.anode(my.mult2,task = "R2")
  my.mult.r2=lapply(my.mult$analyses,pscl::pR2)
  my.mult.r2t=t(as.data.frame(my.mult.r2))
  my.mult.r21=apply(my.mult.r2t,2,mean)

  my.mult.dev=mean(sapply(my.mult$analyses,deviance))
  df.orig=sum(my.mult.anode[,2])
  my.mult.dev.change=my.mult.dev-null.dev
  df.change=df.orig-null.df
  my.p=pchisq(my.mult.dev.change,df.change)

  my.mult.coef=cbind(my.mult.coef,exp(my.mult.coef[,1]),1/exp(my.mult.coef[,1]))

  my.mult.coef=my.mult.coef[,c(1,2,5,6,3,4)]
  if(round){
    my.table=quick.tasks::quick.p.val(my.mult.coef,6)
    my.table[,c(1,2,3,4,5)]=sapply(my.table[,c(1,2,3,4,5)],function(x){round(as.numeric(x),3)})
  }
  #my.table[,2]=formatC(my.table[,2], format = "e", digits = 2)
  colnames(my.table)=c("Estimate","Std. Error","Odds Ratio","1/OR","z Value",'Pr(>|z|)')

  my.anode.tableR=cbind(my.mult.anode,round({{my.mult.anode[,1]/my.mult.dev}},4))
  my.anode.tableR=rbind(my.anode.tableR,c(abs(my.mult.dev.change),df.change,my.p,my.mult.r21[4]))
  rownames(my.anode.tableR)[{length(rownames(my.anode.tableR))}]="Total Change"
  colnames(my.anode.tableR)[4]="McFadden"
  if(round){
    my.anode.tableR=quick.tasks::quick.p.val(my.anode.tableR,3)
    my.anode.tableR[,1]=round(as.numeric(my.anode.tableR[,1]),2)
  }
  my.stats=c(my.mult.dev,df.orig,my.mult.dev.change,df.change,my.p)
  names(my.stats)=c("Deviance","df","Deviance change","Df Change","p-val")

  my.return=NULL
  my.return$ANODE=my.anode.tableR
  my.return$zTable=my.table
  my.return$dev=my.stats
  my.return$pR2=my.mult.r21
  my.return$analyses=my.mult$analyses
  if(include.null){
  my.return$null=my.mult1$analyses
  class(my.return)=c("list","pooled","multinom","null")
  }else{
  class(my.return)=c("list","pooled","multinom")
  }
  return(my.return)
}



#' Print method for multinom tables
#'
#' @param my.print My print
#'
#' @return printed table
#' @export
#'
print.multinom=function(my.print){
  print(my.print$zTable)
}




#' Create multinomial omnibus tables
#'
#' @param ... List of models to be combined. One must have null model.
#'
#' @return Omnibus table
#' @export
#'
quick.table.multi.model=function(...){
  arguments=list(...)
  flag=0

  #Get null model
  for(i in 1:length(arguments)){
    if(sum(class(arguments[[i]])=="null")==1){
      null.models=arguments[[i]]$null
      i=length(arguments)
      flag=1
    }
  }
  if(flag==0){
    stop("No model has the null model.")
  }
  null.dev=sapply(null.models,deviance)
  null.df=null.models[[1]]$edf


  my.table=NULL
  ### Make list of deviance, df to order
  for(i in 1:length(arguments)){
    dev=mean(sapply(arguments[[i]]$analyses,deviance))
    df=arguments[[i]]$analyses[[1]]$edf
    if(i==1){
      my.table=c(i,dev,df,NA,NA,NA)
    }else{
      my.table=rbind(my.table,c(i,dev,df,NA,NA,NA))
    }
  }
  my.table=as.data.frame(my.table)
  my.table=my.table[order(my.table[[3]]),]

  ### Make actual table
  for(i in 1:length(arguments)){
    if(i==1){
      dev.change=my.table[i,2]-mean(null.dev)
      df.change=my.table[i,3]-null.df
    }else{
      dev.change=my.table[i,2]-my.table[i-1,2]
      df.change=my.table[i,3]-my.table[i-1,3]
    }
    my.table[i,4:6]=c(dev.change,df.change,pchisq(dev.change,df.change))
  }
  my.table=rbind(c(0,mean(null.dev),null.df,NA,NA,NA),my.table)

  colnames(my.table)=c("Model","Deviance","df","Dev.Change","df.Change","p-val")
  rownames(my.table)=my.table[[1]]
  my.table[[1]]=NULL

  my.table=quick.tasks::quick.p.val(my.table,5)
  my.table[,3]=round(my.table[,3],2)


  return(my.table)
}
#' #' NOT NEEDED
#' #'
#' #' Has some good code that can be reused later. Not needed atm.
#' #'
#' #' @param imputed.df Mids object from MICE
#' #' @param my.new.df Data frame that has been through label.explor.r
#' #'
#' #' @return Data frame with imputed values averaged and replaced.
#' #' @export
#' #'
#' #' @examples
#' quick.MICE=function(imputed.df,my.new.df){
#'   #### Get rowsums for non-binary variables
#'   #### Unfortunately, does not store variable names.
#'   for(i in 1:length(names(imputed.df$imp))){
#'     if(!is.null(imputed.df$imp[[i]])){
#'       #### Check for factor
#'       if(length(grep("_F$",names(imputed.df$imp)[i]))>0 | length(grep("_O$",names(imputed.df$imp)[i]))>0){
#'         this.grep=grep(paste("^",substr(names(imputed.df$imp)[i],1,{nchar(names(imputed.df$imp)[i])-2}),"$",sep=""),names(my.new.df))
#'         the.names=attr(my.new.df[[this.grep]],"labels")
#'         the.levels=trimws(names(the.names))
#'         my.sums=tryCatch(round(rowMeans(
#'           sapply(imputed.df$imp[[i]][1,],function(y){
#'             sapply(imputed.df$imp[[i]][,y],function(x){
#'               as.numeric(the.names[grep(trimws(as.character(x)),the.levels)][1])})}),na.rm=T)),
#'           error = function(e){warning(e);NULL})
#'         #print(my.sums)
#'       }else{
#'         #### If regular number
#'         this.grep=grep(paste("^",names(imputed.df$imp)[i],"$",sep=""),names(my.new.df))
#'         my.sums=rowMeans(imputed.df$imp[[i]],na.rm=T)
#'       }
#'       if(!is.null(my.sums)){
#'         names(my.sums)=rownames(imputed.df$imp[[i]])
#'         for(j in 1:length(my.sums)){
#'           copy.grep=grep(as.numeric(names(my.sums)[j]),as.numeric(rownames(my.new.df)))
#'           my.new.df[copy.grep,this.grep]=as.numeric(my.sums[j])
#'         }
#'       }else{
#'         print(paste(names(imputed.df$imp)[[i]], "did not work. Please do by hand."))
#'       }
#'     }
#'   }
#'   return(my.new.df)
#' }

#' #' Title
#' #'
#' #' @param my.formula
#' #' @param imputed.df
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' quick.multi=function(my.formula,imputed.df){
#'   library(mice)
#'   library(nnet)
#'   library(AER)
#'   library(pscl
#'   )
#'
#'   #### Make intercept formula
#'   imputed.df.complete=complete(imputed.df)
#'   my.formula=as.character(my.formula)
#'   my.base.form=paste(my.formula[2],my.formula[1],"1")
#'   my.formula=paste(my.formula[2],my.formula[1],my.formula[3])
#'
#'
#'   my.mult1=with(imputed.df,nnet::multinom(eval(parse(text=as.character(my.base.form))),weights=weight,trace=F))
#'   null.dev=mean(sapply(my.mult1$analyses,deviance))
#'   null.df=nnet::multinom(rschl_F~1,trace=F,data=imputed.df.complete)$edf
#'
#'
#'   my.mult=with(imputed.df,nnet::multinom(eval(parse(text=as.character(my.formula))),weights=weight),trace=F)
#'   my.mult4.anode=quick.tasks::quick.multi.anode(my.mult,task = "ANODE")
#'   my.mult4.coef=quick.tasks::quick.multi.anode(my.mult,task = "coef")
#'   #quick.tasks::quick.multi.anode(my.mult2,task = "R2")
#'   my.mult4.r2=lapply(my.mult$analyses,pscl::pR2)
#'   my.mult4.r2t=t(as.data.frame(my.mult4.r2))
#'   my.mult4.r21=apply(my.mult4.r2t,2,mean)
#'
#'   my.mult4.dev=mean(sapply(my.mult$analyses,deviance))
#'   df.change4=sum(my.mult4.anode[,2])
#'   my.mult4.dev.change=my.mult4.dev-null.dev
#'   df.change4=df.change4-null.df
#'   my.p=pchisq(my.mult4.dev.change,df.change4)
#'
#'   my.mult4.coef[5]=exp(my.mult4.coef[,1])
#'   my.mult4.coef[6]=1/exp(my.mult4.coef[,1])
#'   my.mult4.coef=my.mult4.coef[,c(1,2,5,6,3,4)]
#'   my.table=quick.tasks::quick.p.val(my.mult4.coef,6)
#'   my.table[,c(1,3,4,5)]=sapply(my.table[,c(1,3,4,5)],function(x){round(as.numeric(x),3)})
#'   my.table[,2]=formatC(my.table[,2], format = "e", digits = 2)
#'   colnames(my.table)=c("Estimate","Std. Error","Odds Ratio","1/OR","z Value",'Pr(>|z|)')
#'
#'   return(list(my.mult4.anode,my.mult4.coef,my.mult4.r21,list(my.mult4.dev.change,df.change4,my.p)))
#' }
