###################### QUICK.PLOTS #########################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################### 9/2017 ###########################
############################################################


#' Saves plot for later use
#'
#' Turns whatever plot or graph into a plot and saves
#' the object.
#'
#' @return Saved Plot
#' @keywords Explore
#' @export
#' @examples
#' my.plot=quick.grab()
#'
quick.grab=function(){
  Scree=grid::grid.grab()
  Scree
  Scree=recordPlot()
  return(Scree)
}




#' Heat Map of Correlations
#'
#' Helper function to quickly get list of column names
#' based on partial piece.
#'
#' @param myDF Dataframe
#' @param use use option passed to cor()
#' @return GGplot of heat map
#' @keywords Explore
#' @export
#' @examples
#' quick.heat(myDF)
#'
quick.heat=function(myDF,use="pairwise.complet.obs"){
  library(reshape2)
  library(grid)
  library(ggplot2)#for adjusting plot margins
  PreviousCorr=as.data.frame(cor(myDF,use="pairwise.complete.obs"))
  PreviousCorr$name=names(myDF)
  PreviousCorr.m=reshape::melt(PreviousCorr,id="name",variable.name="variable",value.name="value")

  PreviousCorr.m$name=factor(PreviousCorr.m$name,levels=names(myDF))

  #Reverse Order to match
  PreviousCorr.m$variable=factor(PreviousCorr.m$variable,levels=rev(names(myDF)))

  #place the tests on the x- and y-axes,
  #fill the elements with the strength of the correlation
  PreviousP1=ggplot(PreviousCorr.m, aes(name, variable, fill=abs(value))) +
    geom_tile() + #rectangles for each correlation
    #add actual correlation value in the rectangle
    geom_text(aes(label = round(value, 2)), size=2.5) +
    theme_bw(base_size=10) + #black and white theme with set font size
    #rotate x-axis labels so they don't overlap,
    #get rid of unnecessary axis titles
    #adjust plot margins
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.margin = unit(c(3, 1, 0, 0), "mm")) +
    #set correlation fill gradient
    scale_fill_gradient(low="white", high="red") +
    guides(fill=F) #omit unnecessary gradient legend
  return(PreviousP1)
}
