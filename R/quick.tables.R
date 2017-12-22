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
                     test.stat = "Wilks",
                     my.factor = NULL,
                     do.glance=T,
                     adjustment = "bonferroni",
                     show.contrasts=F,
                     show.y.contrasts=F,
                     show.latent=F,
                     show.intercept=F) {
  library(pixiedust)
  library(broom)
  library(car)
  library(tidyr)
  library(phia)
  library(quick.tasks)

  #### Find type
  my.reg.type=quick.type(my.model)

  if(type=="ord"){
    ab.len=30
    library(ordinal)
  }else{
    ab.len=15
  }
  #### Get data frame from parent environment
  my.found.df = eval(parse(text=capture.output(my.model$call$data)),envir = parent.frame())
  #print(dim(my.found.df))
  if (is.null(my.found.df)) {
    stop(paste("No data frame found"))
  }

  #### Make factor list
  #### MANOVA ####
  if (type == "manova2" | type == "stats::manova2") {
    x3 = capture.output(car::Anova(my.model, type = SS.type, test = test.stat))
    my.manova.test = data.frame(matrix(ncol = 7, nrow = 1))
    my.var.temp = 4

    while (my.var.temp < {
      length(x3) - 1
    }) {
      test = strsplit(x3[my.var.temp], "\\s+")

      if (length(test[[1]]) == 9) {
        test2 = test[[1]][-9]
        test2 = test2[-7]

      } else if (length(test[[1]]) == 8) {
        test2 = test[[1]][-8]

      } else{
        test2 = test[[1]]
      }

      my.manova.test[{
        my.var.temp - 3
      }, ] = test2
      my.var.temp = my.var.temp + 1

    }

    my.manova.test[[2]] = as.numeric(my.manova.test[[2]])
    my.manova.test[[3]] = as.numeric(my.manova.test[[3]])
    my.manova.test[[4]] = as.numeric(my.manova.test[[4]])
    my.manova.test[[5]] = as.numeric(my.manova.test[[5]])
    my.manova.test[[6]] = as.numeric(my.manova.test[[6]])
    my.manova.test[[7]] = as.numeric(my.manova.test[[7]])

    options(pixie_interactive = pix.int)
    my.manova.pixie = pixiedust::dust(my.manova.test) %>%
      sprinkle_print_method(pix.method) %>%
      sprinkle(cols = "X7", fn = quote(pvalString(
        value, digits = 3, format = "default"
      ))) %>%
      sprinkle(cols = "X3", round = 3) %>%
      sprinkle_colnames(
        "",
        "df",
        paste(test.stat, " <br /> Statistic"),
        "approx <br /> F-value",
        "num df",
        "den df",
        "Pr(>F)"
      ) %>%
      sprinkle(
        cols = 1:7,
        rows = {
          length(x3) - 5
        },
        border = c("bottom", "left", "right")
      ) %>%
      sprinkle(cols = 1:7, pad = 10) %>%
      sprinkle(
        cols = 1:7,
        rows = 1:{
          length(x3) - 5
        },
        border = c("left", "right")
      ) %>%
      sprinkle(
        cols = 1:7,
        rows = 1,
        border = c("top", "bottom", "left", "right"),
        part = "head"
      ) %>%
      sprinkle(
        cols = 1,
        rows = 1,
        border = "left",
        part = "head"
      ) %>%
      sprinkle(
        cols = 7,
        rows = 1,
        border = "right",
        part = "head"
      ) %>%
      sprinkle_width(
        cols = 1,
        rows = 1:2,
        width = 90,
        width_units = "pt"
      ) %>%
      sprinkle_width(
        cols = 2,
        rows = 1:2,
        width = 30,
        width_units = "pt"
      ) %>%
      sprinkle_width(cols = 3,
                     width = 60,
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
                     width = 70,
                     width_units = "pt") %>%
      sprinkle(rows = 1,
               halign = "center",
               part = "head")


    if (pix.int) {
      return(my.manova.pixie)
    } else{
      my.manova.pixie = print(my.manova.pixie, quote = F)[1]
      return(my.manova.pixie)
    }
  } else if(type=="ord1" | type=="glm1"){
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

  }else{
    #### Use car::Anova to get SS Type 3
    if (type == "lm") {
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
      #### Calculate total SS
      my.total = sum(my.III.summary$`Sum Sq`[2:length(my.III.summary$`Sum Sq`)])
      my.df.total = sum(my.III.summary$Df[2:length(my.III.summary$Df)])
      total.intercepts = 1
      my.rownames = c(abbreviate(rownames(my.summary$coefficients), minlength = abbrev.length),
                      "Residuals",
                      "Total")
      treat.SS=sum(my.III.summary$`Sum Sq`[2:{length(my.III.summary$`Sum Sq`)-1}])
      treat.df=sum(my.III.summary$Df[2:{length(my.III.summary$Df)-1}])
    } else if(type == "manova"){

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

      # options(contrasts=c("contr.sum","contr.poly"))
      # my.new.model=manova(data=Electric3,cbind(HT58,WT58)~FIRSTCHD_F+DAYOFWK)

      for(i in 1:length(my.model$call)){
        my.formula.grep=grep("cbind",my.model$call[[i]])
        if(length(my.formula.grep)>0){
          my.formula=my.model$call[[i]]
        }
      }
      my.vars=strsplit(as.character(my.formula),"~")


      my.new.df=my.model$model
      my.y.levels=dim(my.model$model[[1]])[2]


      #### NULL MODEL ####
      my.null.model=manova(my.new.df[[1]]~1)
      my.seq.model=my.null.model
      my.new.df.inc=my.new.df
      my.new.model=NULL
      my.new.model[[1]]=my.null.model



      #### Type I & II SS ####
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
          #my.j=k
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

      ## Null change
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
      #### Make latent
      #### Make latent variables (~equivalent to ANOVAs) of Full matrix
      my.SSP.type.2.change.fa.eigen=eigen(my.SS.type.2.change[[1]][[1]])
      my.SSP.type.2.change.fa.values=1/sqrt(my.SSP.type.2.change.fa.eigen$values)
      for(n in 2:length(my.SSP.type.2.change.fa.values)){
        my.SSP.type.2.change.fa.values=rbind(my.SSP.type.2.change.fa.values,my.SSP.type.2.change.fa.values)
      }
      my.SSP.type.2.change.fa.values.t=t(my.SSP.total.fa.values)

      ###This
      my.SSP.type.2.change.fa[[1]]=my.SSP.type.2.change.fa.eigen$vectors%*%my.SSP.type.2.change.fa.values.t

      ###This
      my.latent.SSP.type.2.change[[1]]=my.SS.type.2.change[[1]][[1]]%*%my.SSP.type.2.change.fa[[1]]

      # my.SSP.treat.fa.eigen=NULL
      # my.SSP.treat.fa.values=NULL
      # my.SSP.treat.fa.values.t=NULL
      # my.SSP.treat.fa=NULL
      # my.latent.SSP.treat=NULL
      # for(k in 2:{length(my.SS.type.2.change)-1}){
      #   my.SSP.treat.fa.eigen[[k]]=eigen(my.SS.type.2.change[[1]])
      #   my.SSP.treat.fa.values[[k]]=tryCatch(as.matrix({1/sqrt(my.SSP.treat.fa.eigen[[k]]$values)}),warning=function(w){
      #     my.SSP.treat.fa.values[[k]]=matrix(ncol=1,nrow=length(my.SSP.treat.fa.eigen[[k]]$values))
      #     for(j in 1:length(my.SSP.treat.fa.eigen[[k]]$values)){
      #       my.SSP.treat.fa.values[[k]][j]={1/sqrt(max(my.SSP.treat.fa.eigen[[k]]$values[j],0))}
      #       if(my.SSP.treat.fa.values[[k]][j]==Inf){
      #         my.SSP.treat.fa.values[[k]][j]=0
      #       }
      #
      #     }
      #     return(my.SSP.treat.fa.values[[k]])
      #   })
      #   for(j in 2:length(my.SSP.treat.fa.values[[k]])){
      #     my.SSP.treat.fa.values[[k]]=cbind(my.SSP.treat.fa.values[[k]],my.SSP.treat.fa.values[[k]])
      #   }
      #   #my.SSP.treat.fa.values.t[[i]]=t(my.SSP.treat.fa.values[[i]])
      #   my.SSP.type.2.treat.fa[[1]][[k]]=my.SSP.treat.fa.eigen[[k]]$vectors%*%my.SSP.treat.fa.values[[k]]
      #
      #   my.latent.SSP.type.2.change[[1]][[k]]=my.SS.type.2.change[[k]][[1]]%*%my.SSP.type.2.treat.fa[[1]][[k]]
      # }


      for(l in 2:{dim(my.new.df)[2]-1}){
        my.gamma=kappa*l+1
        my.SS.type.2.change[[l]]=summary(my.new.model[[my.gamma]])$SS[the.var.length]
        my.SSP.type.2.change.fa.eigen=eigen(my.SS.type.2.change[[l]][[1]])
        my.SSP.type.2.change.fa.values=1/sqrt(my.SSP.type.2.change.fa.eigen$values)
        for(n in 2:length(my.SSP.type.2.change.fa.values)){
          my.SSP.type.2.change.fa.values=rbind(my.SSP.type.2.change.fa.values,my.SSP.type.2.change.fa.values)
        }
        my.SSP.type.2.change.fa.values.t=t(my.SSP.total.fa.values)

        ###This
        my.SSP.type.2.change.fa[[l]]=my.SSP.type.2.change.fa.eigen$vectors%*%my.SSP.type.2.change.fa.values.t

        ###This
        my.latent.SSP.type.2.change[[l]]=my.SS.type.2.change[[l]][[1]]%*%my.SSP.type.2.change.fa[[l]]

      }


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

      my.y.levels=dim(my.SSP.total)[1]


      #### LATENT MODEL ####
      #### Make latent variables (~equivalent to ANOVAs) of Full matrix
      my.SSP.total.fa.eigen=eigen(my.SSP.total)
      my.SSP.total.fa.values=1/sqrt(my.SSP.total.fa.eigen$values)
      for(i in 2:length(my.SSP.total.fa.values)){
        my.SSP.total.fa.values=rbind(my.SSP.total.fa.values,my.SSP.total.fa.values)
      }
      my.SSP.total.fa.values.t=t(my.SSP.total.fa.values)

      ###This
      my.SSP.total.fa=my.SSP.total.fa.eigen$vectors%*%my.SSP.total.fa.values.t

      ###This
      my.latent.SSP.total=my.SSP.total%*%my.SSP.total.fa

      my.SSP.treat.fa.eigen=NULL
      my.SSP.treat.fa.values=NULL
      my.SSP.treat.fa.values.t=NULL
      my.SSP.treat.fa=NULL
      my.latent.SSP.treat=NULL
      for(i in 1:length(my.SSP.treat)){
        my.SSP.treat.fa.eigen[[i]]=eigen(my.SSP.treat[[i]])
        my.SSP.treat.fa.values[[i]]=tryCatch(as.matrix({1/sqrt(my.SSP.treat.fa.eigen[[i]]$values)}),warning=function(w){
          my.SSP.treat.fa.values[[i]]=matrix(ncol=1,nrow=length(my.SSP.treat.fa.eigen[[i]]$values))
          for(j in 1:length(my.SSP.treat.fa.eigen[[i]]$values)){
            my.SSP.treat.fa.values[[i]][j]={1/sqrt(max(my.SSP.treat.fa.eigen[[i]]$values[j],0))}
            if(my.SSP.treat.fa.values[[i]][j]==Inf){
              my.SSP.treat.fa.values[[i]][j]=0
            }

          }
          return(my.SSP.treat.fa.values[[i]])
        })
        for(j in 2:length(my.SSP.treat.fa.values[[i]])){
          my.SSP.treat.fa.values[[i]]=cbind(my.SSP.treat.fa.values[[i]],my.SSP.treat.fa.values[[i]])
        }
        #my.SSP.treat.fa.values.t[[i]]=t(my.SSP.treat.fa.values[[i]])
        my.SSP.treat.fa[[i]]=my.SSP.treat.fa.eigen[[i]]$vectors%*%my.SSP.treat.fa.values[[i]]

        my.latent.SSP.treat[[i]]=my.SSP.treat[[i]]%*%my.SSP.treat.fa[[i]]
      }

      #### Get SSP treatment and error for Type II error SS

      my.SSP.treat=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSP
      my.SSP.err=car::Anova(my.new.model[[my.y.levels+1]],type=SS.type,test=test.stat)$SSPE

      ### Get total SSP
      my.SSP.treat.total=my.SSP.treat[[1]]
      for(i in 2:{length(my.SSP.treat)}){
        my.SSP.treat.total=my.SSP.treat.total+my.SSP.treat[[i]]
      }
      my.SSP.total=my.SSP.treat.total+my.SSP.err

      my.y.levels=dim(my.SSP.total)[1]

      #### Have incremental model
      #### Make incremental NULL SS
      # my.NULL.cumulative.SS=NULL
      # for(i in 2:dim(my.new.df)[2]){
      #   my.NULL.cumulative.SS[[i-1]]=anova(my.new.model[[i]],my.new.model[[1]])
      # }




      #my.dep.var=my.vars[[2]]
      # my.other.vars=names(new.model$model)
      # my.var.grep=grep(paste("^",my.dep.var,"$",sep=""),my.other.vars)
      # my.other.vars=my.other.vars[-my.var.grep]
      # my.formula=paste("~",my.other.vars[1])
      # for(i in 2:length(my.other.vars)){
      #   my.formula=paste(my.formula,"*",my.other.vars[i],sep="")}




      #### Make error matrices
      my.SSP.err.fa.eigen=eigen(my.SSP.err)
      my.SSP.err.fa.values=1/sqrt(my.SSP.err.fa.eigen$values)
      for(i in 2:length(my.SSP.err.fa.values)){
        my.SSP.err.fa.values=rbind(my.SSP.err.fa.values,my.SSP.err.fa.values)
      }
      my.SSP.err.fa.values.t=t(my.SSP.err.fa.values)
      my.SSP.err.fa=my.SSP.err.fa.eigen$vectors%*%my.SSP.err.fa.values.t

      my.latent.SSP.err=my.SSP.err%*%my.SSP.err.fa



      #### Get all df
      my.SSP.total.df=my.y.levels*dim(my.model$model)[1]-1
      my.SSP.err.df=my.model$df.residual
      my.SSP.treat.df.total=my.SSP.total.df-my.SSP.err.df
      #my.latent.SSP.treat.df.total=my.latent.SSP.total.df-my.latent.SSP.err.df

      my.SSP.treat.df=1
      for(i in 2:length(my.SSP.treat)){
        manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
        if(length(manova.grep)>0){
          my.SSP.treat.df=c(my.SSP.treat.df,{length(my.model$xlevels[[manova.grep]])-1})
        }else{
          my.SSP.treat.df=c(my.SSP.treat.df,1)
        }
      }

      #### Get treatment change totals
      my.SSP.treat.change.total=0
      my.SSP.treat.change=as.data.frame(matrix(ncol=my.y.levels,nrow=my.y.levels))
      my.latent.SSP.treat.change.total=0
      my.latent.SSP.treat.change=as.data.frame(matrix(ncol=my.y.levels,nrow=my.y.levels))
      for(i in 2:length(my.SSP.treat)){
        #my.SSP.treat.change=my.SSP.treat.change+my.SSP.treat[[i]]
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

      my.SSP.treat.change.total=quick.tr(my.SSP.treat.change)
      my.latent.SSP.treat.change.total=quick.tr(my.latent.SSP.treat.change)


      #### Get residual stuff
      the.resid.SS=sum(diag(my.SSP.err))
      the.resid.df=my.model$df.residual

      #### Get total change stuff
      my.total.change=my.SSP.err+my.SSP.treat.change
      the.total.change.SS=the.resid.SS+sum(diag(my.SS.type.1.change.total))
      the.total.change.df=the.resid.df+{my.SSP.treat.change.df/my.y.levels}


      #my.treat.err=solve(my.SSP.err)%*%my.SSP.treat[[length(my.SSP.treat)]]
      #quick.m.test(my.treat.err,"Wilks")

      #### Make basic table ####
      my.table.names=c("var","test.stat","f.val","SS","df","mult.df","resid df","p.val")

      my.manova.table=as.data.frame(matrix(ncol=v.p.len,nrow=1))
      names(my.manova.table)=my.table.names
      #### add intercept
      #my.manova.table[1,]=c("Intercept",NA,NA,sum(eigen(my.SSP.treat[[1]])$values),my.y.levels,my.y.levels*my.SSP.err.df,NA)

      my.line.var=1
      for(i in 1:length(my.SSP.treat)){
        my.treat.err=solve(my.SSP.err)%*%my.SSP.treat[[i]]
        my.test.stat=quick.m.test(my.treat.err,test.stat)
        # my.s=sqrt({my.y.levels*{my.y.levels*my.SSP.treat.df[i]}^2-4}/{my.y.levels^2+{my.y.levels*my.SSP.treat.df[i]}^2-5})
        #
        # my.m=my.SSP.err.df+my.y.levels*my.SSP.treat.df[i]-.5*{my.y.levels+my.y.levels*my.SSP.treat.df[i]+1}
        #
        # my.R.val.part={1-my.test.stat^{1/my.s}}/{my.test.stat^{1/my.s}}
        #
        # my.R.val.part2={my.m*my.s-{my.y.levels*my.y.levels*my.SSP.treat.df[i]}/3}/{my.y.levels*my.y.levels*my.SSP.treat.df[i]}
        #
        # my.R.val=my.R.val.part*my.R.val.part2

        my.SS=quick.tr(my.SSP.treat[[i]])
        my.df=my.y.levels*my.SSP.treat.df[i]
        my.resid.df=min({my.SSP.err.df*my.SSP.treat.df[i]-my.SSP.treat.df[i]},my.y.levels*my.SSP.err.df)
        my.f.val={my.SS/my.df}/{the.resid.SS/my.resid.df}
        my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

        my.manova.table[my.line.var,]=c(names(my.SSP.treat)[i],my.test.stat,my.f.val,my.SS,ifelse(i==1,NA,my.df/my.y.levels),ifelse(i==1,my.y.levels,my.df),my.resid.df,my.p.val)
        my.line.var=my.line.var+1




        #### Put in decomposed factors (y constrasts)
        #### Intercept
        if({show.intercept & i==1}){
          for(y in 1:my.y.levels){
            my.name=paste(rownames(my.SSP.total)[y],sep="")
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
            my.name=paste(y,"|",names(my.SSP.treat)[i],sep="")
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

        #### Put in contrasts
        if(show.contrasts){
          #### Check length
          if({my.SSP.treat.df[i]>1}){
            other.manova.grep=grep(paste("^",names(my.SSP.treat)[i],"$",sep=""),names(my.model$xlevels))
            my.phia=phia::testInteractions(my.model,pairwise = names(my.model$xlevels)[[other.manova.grep]],adjustment = adjustment)

            for(k in 1:{my.SSP.treat.df[i]}){
              my.name=rownames(my.phia)[k]
              my.test.stat=my.phia$`test stat`[k]
              my.f.val=my.phia$`approx F`[k]
              my.SS=sum(my.phia[k,1:{dim(my.phia)[2]-6}])
              my.df=my.phia$`num Df`[k]
              my.p.val=my.phia$`Pr(>F)`[k]
              my.resid.df=my.phia$`den Df`[k]
              my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,my.df,my.resid.df,my.p.val)
              my.line.var=my.line.var+1
            }
          }
        }


        #### Put in latents (ANOVA)
        my.counter=NULL
        for(r in 2:my.y.levels){
          my.counter=c(my.counter,r)
        }
        my.counter=c(my.counter,1)
        if(show.latent &i!=1){

          my.i.temp=my.counter[i-1]
          for(y in 1:{my.y.levels}){
            # my.y=NULL
            # if(y>my.y.levels){
            #   my.y=1
            # }else{
            #   my.y=y
            # }
            my.y=y
            my.name=paste(y,"|",names(my.SSP.treat)[i],sep="")
            my.test.stat=NA
            my.resid.df=my.SSP.err.df


            my.SS=my.latent.SSP.type.2.change[[my.i.temp]][my.y,1]
            my.df=my.SSP.treat.df[my.i.temp]
            my.f.val={my.SS/my.df}/{my.latent.SSP.err[my.y,my.y]/my.resid.df}
            my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

            my.manova.table[my.line.var,]=c(my.name,my.test.stat,my.f.val,my.SS,my.df,NA,my.resid.df,my.p.val)
            my.line.var=my.line.var+1
          }
        }
        #### Put in treatment
        #### From Type I statistics
        if(i==1){
          my.test.stat=quick.m.test(my.SS.type.1.change.total,test.stat)

          my.SS=sum(diag(my.SS.type.1.change.total))
          my.df=my.SSP.treat.change.df
          #### Should really be a min statement, but for later...
          my.resid.df=my.y.levels*my.SSP.err.df

          my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
          my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)

          my.manova.table[my.line.var,]=c("Treatment",my.test.stat,my.f.val,my.SS,NA,my.df,my.resid.df,my.p.val)
          my.line.var=my.line.var+1

          #### Y Contrasts
          if(show.y.contrasts){
            for(b in 1:my.y.levels){
              my.test.stat=NA

              my.SS=my.SS.type.1.change.total[b,b]
              my.df=my.SSP.treat.change.df
              #### Should really be a min statement, but for later...
              my.resid.df=my.y.levels*my.SSP.err.df

              my.f.val={my.SS/my.df}/{sum(diag(my.SSP.err))/my.resid.df}
              my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
              my.manova.table[my.line.var,]=c(paste(b,"|Treatment",sep=""),my.test.stat,my.f.val,my.SS,{my.df/my.y.levels},NA,{my.resid.df/my.y.levels},my.p.val)
              my.line.var=my.line.var+1
            }
          }

          #### Show latent treatments (ANOVAs)
          #### Need to change latents to type II
          my.counter=NULL
          for(r in 2:my.y.levels){
            my.counter=c(my.counter,r)
          }
          my.counter=c(my.counter,1)
          if(show.latent){
            for(b in 1:{my.y.levels}){
              my.SS=sum(weighted.residuals(my.null.model)[,b]^2)-sum(weighted.residuals(my.model)[,b]^2)
              my.df=my.SSP.treat.change.df/my.y.levels
              #### Should really be a min statement, but for later...
              my.resid.df=my.SSP.err.df
              my.f.val={my.SS/my.df}/{my.latent.SSP.err[y,y]/my.resid.df}
              my.p.val=pf(my.f.val,my.df,my.resid.df,lower.tail = F)
              my.manova.table[my.line.var,]=c(paste(b,"|Treatment",sep=""),NA,my.f.val,NA,my.SS,my.df,{my.resid.df},my.p.val)
              my.line.var=my.line.var+1
            }
          }
        }


      }

      #### Put in residuals, total change, total ss
      my.manova.table[my.line.var,]=c("Total Residuals",NA,NA,the.resid.SS,the.resid.df,NA,my.y.levels*the.resid.df,NA)
      my.line.var=my.line.var+1
      if(show.y.contrasts){
        for(i in 1:my.y.levels){
          my.manova.table[my.line.var,]=c(paste(i,"|Residuals",sep=""),NA,NA,my.SSP.err[i,i],the.resid.df,NA,NA,NA)
          my.line.var=my.line.var+1
        }
      }
      if(show.latent){
        for(i in 1:my.y.levels){
          my.manova.table[my.line.var,]=c(paste(i,"|Residuals",sep=""),NA,NA,my.latent.SSP.err[i,i],NA,NA,NA,NA)
          my.line.var=my.line.var+1
        }
      }
      my.manova.table[my.line.var,]=c("Total Change",NA,NA,the.total.change.SS,the.total.change.df,NA,my.y.levels*the.total.change.df,NA)
      my.line.var=my.line.var+1
      my.manova.table[my.line.var,]=c("Total SS",NA,NA,quick.tr(my.SSP.total),the.total.change.df,NA,my.y.levels*the.total.change.df,NA)


      for(i in 2:dim(my.manova.table)[2]){
        my.manova.table[[i]]=as.numeric(my.manova.table[[i]])
      }
      tmp.change=0
      if(show.latent & show.y.contrasts){
        tmp.change=4
      }else if(show.latent & !show.y.contrasts){
        tmp.change=2
      }else if(!show.latent & show.y.contrasts){
        tmp.change=2
      }
      #my.line.var

      #### Make Dusted table ####
      options(pixie_interactive = pix.int,
              pixie_na_string = "")
      my.dust=pixiedust::dust(my.manova.table)%>%
        sprinkle_na_string()%>%
        sprinkle_print_method(pix.method)%>%
        sprinkle_border(cols=1,border="left")%>%
        sprinkle_border(cols={8+v.p.rep},border="right")%>%
        sprinkle_border(rows=my.line.var,boder="bottom")%>%
        sprinkle_border(rows=my.line.var-{2+tmp.change},border="top")%>%
        sprinkle_border(cols=1,border="left",part="head")%>%
        sprinkle_border(cols={8+v.p.rep},border="right",part="head")%>%
        sprinkle_border(rows=1,border=c("top","bottom"),part="head")%>%
        sprinkle_border(rows={1+ifelse(show.intercept,my.y.levels,0)},border="bottom")%>%
        sprinkle_border(rows={1+ifelse(show.y.contrasts | show.latent,my.y.levels+1,1)+
            ifelse(show.intercept,my.y.levels,0)},
                        border="bottom")%>%
        sprinkle_round(cols=2:v.p.len,round=3)%>%
        sprinkle_colnames("Variable",paste(test.stat, "<br /> Test Statistic",sep=""),
                          "F-Value",paste("Sums of <br /> Squares",sep=""),"dF",
                          "Mult. <br /> dF","Resid <br /> dF","P-value")%>%
        sprinkle_align(rows=1,halign="center",part="head")%>%
        sprinkle_pad(rows=1:{my.line.var+v.p.rep},pad=5)%>%
        sprinkle(cols = "p.val", fn = quote(pvalString(
          value, digits = 3, format = "default"
        )))%>%
        sprinkle_border(rows={1+ifelse(show.intercept,my.y.levels,0)},border="bottom")

      ##### Make glance stats
      my_glance_stats=as.data.frame(matrix(ncol=v.p.len,nrow=1))
      my_glance_stats[1,]=c(paste("Method: QR decomposition",if(show.contrasts){paste(" <br />Adjustment: ", adjustment,sep="")},if(show.latent){paste(" <br /> Latent Contrasts")}),rep(NA,{7+v.p.rep}))


      my.dust=pixiedust::redust(my.dust,my_glance_stats,part="foot")%>%
        sprinkle(merge=T,halign="center",part="foot")

      if (pix.int) {
        return(my.dust)
      } else{
        my.dust.print = print(my.dust, quote = F)[1]
        return(my.dust.print)
      }


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
    }else if (type == "glm") {
      new.df=my.model$model
      new.model=update(my.model,data=new.df)
      null.model=update(new.model,~1)

      total.intercepts = 1

      vars.df=new.model$df.null-new.model$df.residual
      # vars.df=new.model$df.null-new.model$df.residual
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

    #### Make table ####
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
