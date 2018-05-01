
##################### LABEL.EXPLOR.R #######################
#################### BY CHRIS KRANER #######################
############## NORTHERN ILLINOIS UNIVERSITY ################
######################### 9/2017 ###########################
############################################################




#' Write table of labels
#'
#' For label.explor.r
#'
#' @param myDF Dataframe
#' @param myxnum Column number of variable of interest
#' @return Table
#' @keywords Explore
#' @examples
#' myfinaltable()

myfinaltable=function(myxnum,myDF){
  #if(attr(myDF[[myxnum]],"labels")[1]){
  mytable2=table(myDF[[myxnum]],useNA="always")
  mylabs=rownames(mytable2)
  j=1
  myfinaltable=matrix(c(as.numeric(mylabs[j]),mytable2[[j]]),ncol=2)
  j=j+1
  while(j<{length(mylabs)+1}){

    myfinaltable=rbind(myfinaltable,c(as.numeric(mylabs[j]),mytable2[[j]]))
    #myfinaltable=rbind(myfinaltable,c(as.numeric(mylabs[j]),mytable2[[j]]))
    j=j+1
  }

  return(myfinaltable)
}

#' Write 'em both!
#'
#' For label.explor.r. Outputs as mytable1, mytable2 from tempnum
#'
#' @param tempnum Number for temp
#' @param myvar List of variable names
#' @param myLabels List of Label Names
#' @param myDF Actual DF
#' @param myDFname Name of DF
#' @param cutoff Number of variables no longer a factor
#' @param should.comma needs a comma or is last?
#' @return File
#' @keywords Explore
#' @examples
#' mytempfactor.ui.server()
mytempfactor.ui.server=function(tempnum,myvar,myLabels,myDF,myDFname,cutoff,should.comma){


  my.temp.ui=tempfile()


  tablename=capture.output(cat("mytable",tempnum,sep=""))
  textInputName=capture.output(cat("NA",tempnum,sep=""))
  add_to_list_num=capture.output(cat("add_to_list_",tempnum,sep=""))
  add_to_factor_list_num=capture.output(cat("add_to_factor_list_",tempnum,sep=""))
  make_ordinal_num=capture.output(cat("make_ordinal_",tempnum,sep=""))
  labelname=capture.output(cat("mylabel",tempnum,sep=""))
  wellname=capture.output(cat("mywell",tempnum,sep=""))
  factorname=capture.output(cat("myfactor",tempnum,sep=""))
  mynum4table=grep(paste("^",myvar[tempnum],"$",sep=""),colnames(myDF))


  myDFname4labels=capture.output(cat(myDFname,"$",myvar[tempnum],sep=""))

  myinputs=NULL
  mylabs=NULL
  my.catch=NULL
  my.label.grep=0
  mytable=table(myDF[[mynum4table]])
  mylabs=rownames(mytable)
  # library(bigtabulate)
  # mytable=bigtable(myDF,myvar[tempnum])
  # mytable=table(myDF[[mynum4table]])
  # mylabs=rownames(mytable)
  if(length(attr(myDF[[mynum4table]],"labels"))>0){
    mylabs2=attr(myDF[[mynum4table]],"labels")
    attributes(mylabs2)=NULL
    my.label.grep=tryCatch(grep(mylabs[1],mylabs2),warning=function(w){NULL}, error=function(e){NULL})
    #my.label.grep=grep(mylabs[1],mylabs2)
  }




  # mynum4table=12
  # myDF=HBSAC
  # cutoff=5
  # tempnum=1
  #
  #### Parse Variable-level responses
  if(length(mylabs)>0){

    if(length(mylabs)==1){
      myinputs=capture.output(cat(",p(class=\"list\",
                                  textInput(\"labsInput",tempnum,"_",1,"\",\"",mylabs[1],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[my.label.grep])[1],
                                  "\"))",sep=""))
    }else if(length(mylabs)<{cutoff+1}){
      wtemp=1
      xtemp=1
      flag.this=0
      while(wtemp<{length(mylabs)}){
        if(wtemp==1){
          if(is.numeric(attr(myDF[[mynum4table]],"labels")[[1]])){
            if(length(as.numeric(attr(myDF[[mynum4table]],"labels")))>0){
              xtemp=my.label.grep
              myinputs=capture.output(cat(",p(class=\"list\",
                                          textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[as.numeric(xtemp)])[1],
                                          "\")),",sep=""))
              flag.this=1
            }else{
              myinputs=capture.output(cat(",p(class=\"list\",
                                          textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[as.numeric(wtemp)])[1],
                                          "\")),",sep=""))
            }}else{
              myinputs=capture.output(cat(",p(class=\"list\",
                                          textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[xtemp])[1],
                                          "\")),",sep=""))
            }
        }else{
          #this.Label=capture.output()
          myinputs=capture.output(cat(myinputs,"textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[as.numeric(xtemp)])[1],"\"),",sep=""))
        }
        wtemp=wtemp+1
        xtemp=xtemp+1
      }
      if(flag.this==1){
        myinputs=capture.output(cat(myinputs,"textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[xtemp])[1],"\")",sep=""))
      }else{
        myinputs=capture.output(cat(myinputs,"textInput(\"labsInput",tempnum,"_",wtemp,"\",\"",mylabs[wtemp],"\",\"",capture.output(attr(myDF[[mynum4table]],"labels")[wtemp])[1],"\")",sep=""))
      }
    }else{
      myinputs=NULL
    }
  }else{
    myinputs=NULL
  }
  #### Parse DF-level responses
  if(tempnum==1){
    my.temp.ui.file=capture.output(cat("fluidRow(shinyjs::hidden(span(id=\"",wellname,"\",
                                       box(title=\"",myvar[tempnum],"\", status=\"",if(length(mylabs[1])>0){if(mylabs[1]<0){"danger"}else if(length(mytable)<{cutoff+1}){"success"}else{"primary"}}else{"danger"},"\",solidHeader=T,collapsible=T,collapsed=T,
                                       tabBox(width=12,tabPanel(\"Info\",{div(textInput(\"nameInput",tempnum,"\",\"\",\"",myLabels[tempnum],"\"),
                                       tableOutput(\"",tablename,"\"),
                                       tableOutput(\"",labelname,"\"),
                                       div(capture.output(cat(\"Please enter one value at a time\")),
                                       textInput(\"",textInputName,"\",\"\",width=\'20%\'),
                                       actionButton(\"",add_to_list_num,"\",\"Change to NA\"),align=\"center\"),
                                       align=\"center\"
                                       )}),
                                       tabPanel(\"Factors\",",
                                       "h4(\"Factors\"),p(\"Please label any missing as Missing if removing variables\")",myinputs,",
                                       checkboxInput(\"",make_ordinal_num,"\",\"Is Ordered?\"),
                                       div(actionButton(\"",add_to_factor_list_num,"\",\"Make Factor\"),align=\"right\")
                                       ))))),",sep=""))
  }else if(should.comma){
    my.temp.ui.file=capture.output(cat("shinyjs::hidden(div(id=\"",wellname,"\",
                                       box(title=\"",myvar[tempnum],"\", status=\"",if(length(mylabs[1])>0){if(mylabs[1]<0){"danger"}else if(length(mytable)<{cutoff+1}){"success"}else{"primary"}}else{"danger"},"\",solidHeader=T,collapsible=T,collapsed=T,
                                       tabBox(width=12,tabPanel(\"Info\",{div(textInput(\"nameInput",tempnum,"\",\"\",\"",myLabels[tempnum],"\"),
                                       tableOutput(\"",tablename,"\"),
                                       tableOutput(\"",labelname,"\"),
                                       div(capture.output(cat(\"Please enter one value at a time\")),
                                       textInput(\"",textInputName,"\",\"\",width=\'20%\'),
                                       actionButton(\"",add_to_list_num,"\",\"Change to NA\"),align=\"center\"),
                                       align=\"center\"
                                       )}),
                                       tabPanel(\"Factors\",",
                                       "h4(\"Factors\"),p(\"Please label any missing as Missing if removing variables\")",myinputs,",
                                       checkboxInput(\"",make_ordinal_num,"\",\"Is Ordered?\"),
                                       div(actionButton(\"",add_to_factor_list_num,"\",\"Make Factor\",class=\"factors\"),align=\"right\"))
                                       )))),",sep=""))
  }else{
    my.temp.ui.file=capture.output(cat("shinyjs::hidden(div(id=\"",wellname,"\",
                                       box(title=\"",myvar[tempnum],"\", status=\"",if(length(mylabs[1])>0){if(mylabs[1]<0){"danger"}else if(length(mytable)<{cutoff+1}){"success"}else{"primary"}}else{"danger"},"\",solidHeader=T,collapsible=T,collapsed=T,
                                       tabBox(width=12,tabPanel(\"Info\",{div(textInput(\"nameInput",tempnum,"\",\"\",\"",myLabels[tempnum],"\"),
                                       tableOutput(\"",tablename,"\"),
                                       tableOutput(\"",labelname,"\"),
                                       div(capture.output(cat(\"Please enter one value at a time\")),
                                       textInput(\"",textInputName,"\",\"\",width=\'20%\'),
                                       actionButton(\"",add_to_list_num,"\",\"Change to NA\"),align=\"center\"),
                                       align=\"center\"
                                       )}),
                                       tabPanel(\"Factors\",",
                                       "h4(\"Factors\"),p(\"Please label any missing as Missing if removing variables\")",myinputs,",
                                       checkboxInput(\"",make_ordinal_num,"\",\"Is Ordered?\"),
                                       div(actionButton(\"",add_to_factor_list_num,"\",\"Make Factor\",class=\"factors\"),align=\"right\"))
                                       )))))",sep=""))
  }


  writeLines(my.temp.ui.file,my.temp.ui)


  myinput.labs1=NULL
  my.new.var.names=NULL


  my.temp.server.output=tempfile()
  my.new.names=tempfile()




  wtemp=1

  my.input.labs1=capture.output(cat("shinyjs::onclick(\"add_to_factor_list_",tempnum,"\",{label.explor.r.temp.list=c(",sep=""))
  if(length(mylabs)<{cutoff+1}){
    while(wtemp<{length(mylabs)}){
      my.input.labs1=capture.output(cat(my.input.labs1,"input$labsInput",tempnum,"_",wtemp,",",sep=""))
      wtemp=wtemp+1
    }}else{
      my.input.labs1=capture.output(cat(my.input.labs1,"input$labsInput",tempnum,"_",wtemp,",",sep=""))
      wtemp=wtemp+1
    }
  my.input.labs1=capture.output(cat(my.input.labs1,"input$labsInput",tempnum,"_",wtemp,");if(!input$make_ordinal_",tempnum,"){
                                    label.explor.r.my.factor.list$Name<<-c(label.explor.r.my.factor.list$Name,\"",myvar[tempnum],"\");
                                    label.explor.r.my.factor.list$Value<<-c(label.explor.r.my.factor.list$Value,list(label.explor.r.temp.list[]))}else{
                                    label.explor.r.my.ordered.list$Name<<-c(label.explor.r.my.ordered.list$Name,\"",myvar[tempnum],"\");
                                    label.explor.r.my.ordered.list$Value<<-c(label.explor.r.my.ordered.list$Value,list(label.explor.r.temp.list[]))
                                    }})",sep=""))


  my.temp.server.file=capture.output(cat("observe({
                                         if(input$add_to_list_",tempnum,"){
                                         #print(\"I'm here\")
                                         isolate(
                                         if(length(label.explor.r.myList)>0){
                                         label.explor.r.myList<<-rbind(label.explor.r.myList,c(\"",myvar[tempnum],"\",as.numeric(input$NA",tempnum,")))
                                         #print(myList)
                                         }else{
                                         label.explor.r.myList<<-matrix(c(\"",myvar[tempnum],"\",as.numeric(input$NA",tempnum,")),ncol=2)
                                         })
                                         shinyjs::disable(capture.output(cat(\"input$labsInput1_",tempnum,"\")))

                                         #   print(label.explor.r.myList)
                                         }})","\n","observeEvent(input$nameInput",tempnum,",{
                                         if(length(attr(myDF[[",tempnum,"]],\"this.flag\")>0)){
                                         attr(myDF[[",tempnum,"]],\"this.flag\")<<-",{tempnum},"
                                         nameInput",{tempnum},"<<-input$nameInput",tempnum,"
                                         }else if(length(attr(myDF[[",tempnum,"]],\"this.flag\"))==0){
                                             attr(myDF[[",tempnum,"]],\"this.flag\")<<-99
                                         }})","\n","output$mytable",tempnum,"=renderTable({
                                         myxnum=grep(paste(\"^\",\"",myvar[tempnum],"\",\"$\",sep=\"\"),colnames(",myDFname,"))
                                         myfinaltable(myxnum,",myDFname,")
                                         }, rownames = F,colnames=F
                                         )","\n",
                                         "output$mylabel",tempnum,"=renderTable({
                                         i=1
                                         if(!is.null(attr(",myDFname,"$",myvar[tempnum],",\"labels\", exact=T)[1])){
                                         while(i<{length(attr(",myDFname,"$",myvar[tempnum],",\"labels\", exact=T))+1}){
                                         if(i==1){
                                         x=capture.output(attr(",myDFname,"$",myvar[tempnum],",\"labels\", exact=T)[1])
                                         x=as.matrix(x)
                                         }else{
                                         x=cbind(x,capture.output(attr(",myDFname,"$",myvar[tempnum],",\"labels\", exact=T)[i]))
                                         }
                                         i=i+1
                                         }
                                         x[2,]=as.numeric(x[2,])
                                         x
                                         }else{
                                         return(NULL)
                                         }
                                         },colnames=F)
                                         output$myfactor",tempnum,"=renderTable({
                                         mytable=table(",myDFname,"$",myvar[tempnum],")
                                         mylabs=rownames(mytable)},colnames=F)","\n",my.input.labs1,sep=""))
  writeLines(my.temp.server.file,my.temp.server.output)

  mytemps=NULL
  mytemps[1]=my.temp.ui
  mytemps[2]=my.temp.server.output
  return(mytemps)
  }



#' Write temp table Server interface
#'
#' For label.explor.r. Outputs as mytable1, mytable2 from tempnum Render table and add to list
#'
#' @param tempnum Number for temp
#' @param myvar string variable name
#' @param myDFname string data frame name
#' @return List of files, first being UI, second being server
#' @keywords Explore
#' @examples
#' make.me.shiny()


make.me.shiny=function(myvar,myDF.name,myLabels,myDF2,mycutoff){
  tempnum=1
  while(tempnum<{length(myvar)})
  {
    if(tempnum==1){
      my.temp.files1=mytempfactor.ui.server(tempnum,myvar,myLabels,myDF2,myDF.name,cutoff=mycutoff,should.comma=T)
      my.temp.ui=my.temp.files1[1]
      my.temp.server=my.temp.files1[2]
      # my.temp.ui=mytempfactor.ui(tempnum,T,myvar,myLabels,myDF2)
      # my.temp.server=mytempfactor.server(tempnum,myvar[tempnum],myDF.name,myDF2)
    }else{
      my.temp.files1=mytempfactor.ui.server(tempnum,myvar,myLabels,myDF2,myDF.name,cutoff=mycutoff,should.comma=T)
      my.temp.ui[tempnum]=my.temp.files1[1]
      my.temp.server[tempnum]=my.temp.files1[2]
      # my.temp.ui[tempnum]=mytempfactor.ui(tempnum,T,myvar,myLabels,myDF2)
      # my.temp.server[tempnum]=mytempfactor.server(tempnum,myvar[tempnum],myDF.name,myDF2)
    }
    tempnum=tempnum+1
  }

  my.temp.files1=mytempfactor.ui.server(tempnum,myvar,myLabels,myDF2,myDF.name,cutoff=mycutoff,should.comma=F)
  my.temp.ui[tempnum]=my.temp.files1[1]
  my.temp.server[tempnum]=my.temp.files1[2]

  tempnum2=1
  while(tempnum2<{length(myvar)+1})
  {
    if(tempnum2==1){
      my.main.server.file=readLines(my.temp.server[1])
      my.main.ui.file=readLines(my.temp.ui[1])
    }else{
      my.main.server.file=c(my.main.server.file,readLines(my.temp.server[tempnum2]))
      my.main.ui.file=c(my.main.ui.file,readLines(my.temp.ui[tempnum2]))
    }
    tempnum2=tempnum2+1
  }

  # my.ui.source=tempfile()
  # my.server.source=tempfile()
  # writeLines(my.main.server.file,my.server.source)
  # writeLines(my.main.ui.file,my.ui.source)

  my.source=tempfile()
  my.source[2]=tempfile()
  writeLines(my.main.ui.file,my.source[1])
  writeLines(my.main.server.file,my.source[2])

  return(my.source)
}




#' SPSS Label Interface (label.explor.R)
#'
#' Interface for working with dataframes to prepare for various
#' analyses using stats.explor.r. Creates labels, can remove
#' erroneous variables, makes factors, etc.
#'
#' @param myDF Dataframe
#' @param caseid String name of column of case ids
#' @param mycutoff Optional. Change how many levels a factor can have
#' @return Dataframe with factors, added attributes. Ready for use in stats.explor.r
#' @keywords Explore
#' @export
#' @examples
#' label.explor.r()


label.explor.r=function(myDF,caseid,mycutoff=5,objections=NULL){
  library(shiny)
  library(shinyjs)
  library(car)
  library(shinydashboard)



  #### Function from stats.explor.r.utils, copied here to get rid of dependency ####
  writetotemp_ui = function(inputType,mylabels,myvars,order,inputVariable,desclabel,inlinetest)
  {

    if(order[1]==2){
      order[1]=3
    }
    if(order[2]==2){
      order[2]=3
    }

    mytempfile=tempfile()

    mylabels2=mapply(c,mylabels, "=" ,myvars)

    beginning=switch(inputType,
                     "radio"="radioButtons(",
                     "select"="selectInput(",
                     "check"="checkboxGroupInput(")

    j=1
    writeLines(c(capture.output(cat(beginning,
                                    "\"",inputVariable,"\",",
                                    "\"",desclabel,"\",",
                                    "c(",sep="")),
                 capture.output(while(j < dim(mylabels2)[2]){
                   cat("\"",mylabels2[order[1],j],"\"",mylabels2[2,j],"\"",mylabels2[order[2],j],"\"",",","\n",sep="")
                   j=j+1
                 }),
                 capture.output(cat("\"",mylabels2[order[1],j],"\"",mylabels2[2,j],"\"",mylabels2[order[2],j],"\"",")",sep="")),
                 {if(inlinetest){capture.output(cat(",inline=T"))}},
                 capture.output(cat(")"))
    ),

    mytempfile)

    return(mytempfile)
  }

  #### INITS ####
  myVars2=colnames(myDF)
  myLabels=NULL
  label.explor.r.my.factor.list=NULL
  label.explor.r.my.factor.list.f=NULL
  label.explor.r.my.factor.list.f.labels=NULL
  label.explor.r.my.ordered.list=NULL
  label.explor.r.my.ordered.list.f=NULL
  label.explor.r.my.ordered.list.f.labels=NULL
  #my.temp.list=NULL
  myDF.name="myDF2"
  label.explor.r.my.factor.list=NULL
  label.explor.r.my.factor.list$Name="Name"
  mystartlist=c("Value 1","Value 2")
  label.explor.r.my.factor.list$Value=list(mystartlist[])
  label.explor.r.myList=NULL
  label.explor.r.temp.list=NULL
  flag=F
  flag=as.data.frame(flag)
  label.explor.r.my.selected.names=NULL
  flag.factor=NULL
  labs.flag=NULL
  myDF2=myDF
  myVars=colnames(myDF2)
  caseidnumber=grep(caseid,myVars2)
  for(i in 1:length(myVars2)){
    eval(parse(text=paste("nameInput",i,"=NULL",sep="")))
  }

  #### Fail Checks ####
  if(length(caseidnumber)==0){
    stop("caseid not found")
  }



  #### Deal with Y/N Factors Problem ####
  mytemp3=1
  while(mytemp3<{length(myVars2)+1}){
    if({typeof(myDF[[mytemp3]])=="character"} && {myDF[[mytemp3]][1]=="Y" || myDF[[mytemp3]][1]=="N"}){
      tempvar=car::recode(myDF[[mytemp3]],"'N'=0;'Y'=1",as.factor.result=F)
      templabels=attr(myDF[[mytemp3]],"labels",exact=T)

      #### Make new named labels
      if(templabels[[1]]=="N"){
        new.templabels=c(0,1)
      }else{
        new.templabels=c(1,0)
      }
      attr(new.templabels,"names")=attr(templabels,"names")
      templabels=new.templabels

      if(length(attr(myDF[[mytemp3]],"label",exact=T))>0){
        templabel=attr(myDF[[mytemp3]],"label",exact=T)
      }else{
        templabel=names(myDF)[mytemp3]
      }
      tempname=colnames(myDF[mytemp3])
      myDF2=myDF[-mytemp3]
      myDF2[[{length(myDF)}+1]]=tempvar
      names(myDF2)[length(myDF2)]=tempname
      attr(myDF2[[length(myDF2)]],"labels")=templabels
      attr(myDF2[[length(myDF2)]],"label")=templabel
    }
    if(!is.null(attr(myDF[[mytemp3]],"label",exact=T))){
      myLabels[mytemp3]=attr(myDF[[mytemp3]],"label", exact=T)[1]
    }else{
      myLabels[mytemp3]=myVars2[mytemp3]
    }
    mytemp3=mytemp3+1
  }
  #myVars=colnames(myDF2)

  #### Make sure all numbers are numbers ####
  #### Sometimes they come in as strings....
  my.int.list=1:dim(myDF)[2]
  if(length(objections)>0){
    for(i in 1:length(objections)){
      objection.grep=grep(objections[i],names(myDF))
      my.int.list=my.int.list[-objection.grep]
    }
  }
  for(i in my.int.list){
    myDF[[i]]=as.numeric(myDF[[i]])
  }

  #### Make the big sidepanel checkbox ####
  myVars3=myVars2[-caseidnumber]
  myLabels3=myLabels[-caseidnumber]
  mybigcheckbox=writetotemp_ui("check",myLabels3,myVars3,c(1,2),"check1","Variables",inlinetest = F)

  #### Make the rest of the sources ####
  myTemps=make.me.shiny(myVars,"myDF2",myLabels,myDF2,mycutoff)
  my.ui.source=myTemps[1]
  my.server.source=myTemps[2]

  library(shinyjqui)

  jsCode="shinyjs.ClickVars=function(){var x=document.getElementsByClassName(\"factors\");for(i = 0; i<x.length; i++){x[i].click();}}"
  ##### Actual Shiny App ####
  label.explor.r.my.selected=runApp(list(
    ui=dashboardPage(dashboardHeader(title="label.explor.R"),
                     dashboardSidebar(


                       tags$style(HTML("#checkdiv {overflow-y:auto;height:100vh}")),
                       # div(style="display:inline-box",align="right",
                       #      actionButton("commit","Select")),
                       #      div(style="display:inline-box",actionButton("reset","Reset")),
                       # p(),
                       tags$style("#check1 {font-size:18px;height:18px}"),
                       div(id="checkdiv",source(mybigcheckbox[1],local=T)[1])
                       # div(actionButton("reset","Reset"),
                       # actionButton("commit","Select"),align="right"),
                       # hr(),
                       # div(actionButton("save","Save"),align="center")
                     ),
                     dashboardBody(useShinyjs(),
                                   extendShinyjs(text=jsCode),
                                   id="ViewPanel",
                                   tags$style(HTML("#commit {background-color: white;color: black;border: 2px solid #008CBA;padding: 16px 32px;
                                                   text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;
                                                   -webkit-transition-duration: 0.4s; /* Safari */
                                                   transition-duration: 0.4s;cursor: pointer;}")),
                                   tags$style(HTML("#commit:hover {    background-color: #008CBA;color: white;}")),
                                   tags$style(HTML("#selectall {background-color: white;color: black;border: 2px solid #008CBA;padding: 16px 32px;
                                                   text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;
                                                   -webkit-transition-duration: 0.4s; /* Safari */
                                                   transition-duration: 0.4s;cursor: pointer;}")),
                                   tags$style(HTML("#selectall:hover {    background-color: #008CBA;color: white;}")),
                                   tags$style(HTML("#allfactors {background-color: white;color: black;border: 2px solid #008CBA;padding: 16px 32px;
                                                   text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;
                                                   -webkit-transition-duration: 0.4s; /* Safari */
                                                   transition-duration: 0.4s;cursor: pointer;}")),
                                   tags$style(HTML("#allfactors:hover {    background-color: #008CBA;color: white;}")),
                                   tags$style(HTML("#reset {background-color: white;color: black;border: 2px solid #f44336;padding: 16px 32px;
                                                   text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;
                                                   -webkit-transition-duration: 0.4s; /* Safari */
                                                   transition-duration: 0.4s;cursor: pointer;}")),
                                   tags$style(HTML("#reset:hover {    background-color: #f44336;color: white;}")),
                                   tags$style(HTML("#save {background-color: white;color: black;border: 2px solid #4CAF50;padding: 16px 32px;
                                                   text-align: center;text-decoration: none;display: inline-block;font-size: 16px;margin: 4px 2px;
                                                   -webkit-transition-duration: 0.4s; /* Safari */
                                                   transition-duration: 0.4s;cursor: pointer;}")),
                                   tags$style(HTML("#save:hover {    background-color: #4CAF50;color: white;}")),
                                   div(style="display:inline-box",span(align="right",actionButton("commit","Update"),actionButton("selectall","Select All"),actionButton("allfactors","Make All Factors")),
                                       span(align="right",
                                            actionButton("save","Save"),actionButton("reset","Reset")
                                       )),
                                   p(),
                                   hr(),
                                   source(my.ui.source[1],local=T)[1]
                                   )
                                   ),

    server=function(input,output,session){

      #jqui_resizable('.main-sidebar',options = list(alsoResize=".content-wrapper"))
      #jqui_resizable('.content-wrapper',options=list(alsoResize="#ViewPanel"))
      #jqui_resizable('#ViewPanel')
      currentChecks=reactive({input$check1})

      shinyjs::onclick("reset",
                       {shinyjs::reset("check1")
                       })
      shinyjs::onclick("selectall",
                       {updateCheckboxGroupInput(session,"check1","Variables:",selected=myVars3)
                       })
      shinyjs::onclick("allfactors",{
        js$ClickVars()
      })
      shinyjs::onclick("save",{
        print("START SAVE")
        label.explor.r.my.selected.variables=input$check1
        # for(i in 1:length(label.explor.r.my.selected.variables)){
        #   new.name.grep=grep(paste("^",label.explor.r.my.selected.variables[i],"$",sep=""),names(myDF2))
        #   eval(parse(text=paste("nameInput",new.name.grep,sep="")))=eval(parse(text=paste("input$nameInput",new.name.grep,sep="")))
        # }
        stopApp(label.explor.r.my.selected.variables)



      })





      #### Hide and Show Wells for correct variables ####
      observe({
        if(input$commit>0){
          isolate({
            thistemp=1
            while(thistemp<{length(myVars)+1}){
              thisarg=capture.output(cat("mywell",thistemp,sep=""))
              shinyjs::hide(thisarg)
              thistemp=thistemp+1
            }
            i=1
            while(i<{length(input$check1)+1}){
              mynum=grep(paste("^",input$check1[[i]],"$",sep=""),myVars)
              myarg=capture.output(cat("mywell",mynum,sep=""))
              shinyjs::show(myarg)
              i=i+1
            }
          })

        }

      })

      #### Rest of the server code ####
      source(my.server.source[1],local=T)
    }))

  #### Make lists global ####
  label.explor.r.my.selected <<- label.explor.r.my.selected
  label.explor.r.myList<<-label.explor.r.myList
  label.explor.r.my.factor.list<<-label.explor.r.my.factor.list

  #my.file=RCurl::clone(file("stdin"))
  #my.file2=file.choose()
  #writeLines(my.file,my.file2)
  #### Make Subset ####
  if(length(grep(caseid,label.explor.r.my.selected)>0)){
    myDF2=myDF2[label.explor.r.my.selected]
  }else{
    myDF2=myDF2[c(label.explor.r.my.selected,caseid)]
  }


  #### Apply any NA Operations ####
  myList2=as.data.frame(label.explor.r.myList)

  i=1

  while(i<{dim(myList2)[1]+1}){
    da.other=myList2[i,1]
    da.other=levels(da.other)
    da.num=myList2[i,2]
    da.num=as.numeric(levels(da.num))
    da.grep=grep(paste("^",da.other,"$",sep=""),names(myDF2))


    TFList=lapply(myDF2[[da.grep]], function(r) any(r %in% da.num))
    #TFList

    last.temp=1

    while(last.temp<{length(TFList)+1}){

      if(TFList[[last.temp]]){

        myDF2[last.temp,da.grep]=NA

      }

      last.temp=last.temp+1

    }
    this.stinks={length(flag)+1}
    flag[this.stinks]=T
    colnames(flag)[this.stinks]=paste(da.other)
    i=i+1
  }


  #### Change Variable Labels if applicable ####
  for(i in 1:length(label.explor.r.my.selected)){
    whomp.grep=grep(paste("^",label.explor.r.my.selected[i],"$",sep=""),names(myDF2))
    if(length(attr(myDF2[[i]],"this.flag"))>0){
      if(attr(myDF2[[i]],"this.flag")!=99){
        if(attr(myDF2[[i]],"this.flag")>1){
          attr(myDF2[[i]],"label")=eval(parse(text=paste("nameInput",attr(myDF2[[i]],"this.flag"),sep="")))
          label.explor.r.temp.list=c(label.explor.r.temp.list,paste("nameInput",attr(myDF2[[i]],"this.flag"),sep=""))
        }
      }
      attr(myDF2[[i]],"this.flag")=NULL
    }
  }


  #### Make Factors ####
  last.last.temp=2
  while(last.last.temp<{length(label.explor.r.my.factor.list[[1]])+1}){

    da.last.grep=grep(paste("^",label.explor.r.my.factor.list[[1]][last.last.temp],"$",sep=""),names(myDF2))
    if(length(da.last.grep)>0){
      da.levels=label.explor.r.my.factor.list[[2]][[last.last.temp]]
      darn.lab.temp=0
      for(darn.i in 1:length(da.levels)){
        darn.lab.temp=darn.lab.temp+ifelse(da.levels[darn.i]=="NULL" | is.null(da.levels[darn.i]),0,1)
      }
      #print(da.levels)
      darn.length.temp=length(table(myDF2[[da.last.grep]]))
      #print(darn.lab.temp)
      #print(darn.length.temp)
      if(length(darn.lab.temp) >0){

        if(darn.lab.temp==darn.length.temp){
          if(length(da.levels)>0 & length(da.levels)<{mycutoff+1}){
            flag.grep=grep(paste("^",label.explor.r.my.factor.list[[1]][last.last.temp],"$",sep=""),names(flag))

            if(length(flag.grep)>0){
              da.last.grep2=grep(paste("^",label.explor.r.my.factor.list[[1]][last.last.temp],"$",sep=""),names(myDF))
              mytable=table(myDF[[da.last.grep2]])
              mylabs=as.numeric(rownames(mytable))
              mytable2=table(myDF2[[da.last.grep]])
              mylabs2=as.numeric(rownames(mytable2))
              my.diff=setdiff(mylabs,mylabs2)

              for.caution=1
              while(for.caution<{length(my.diff)+1}){
                remove.grep=grep(my.diff,mylabs)
                da.levels=da.levels[-remove.grep]
                for.caution=for.caution+1
              }
            }

            #Check for missing
            da.levels2=tolower(da.levels)
            missing.grep=grep("missing",da.levels2)
            if(length(missing.grep)>0){
              da.levels=da.levels[-missing.grep]
            }


            end.num=length(colnames(myDF2))+1
            myDF2[end.num]=factor(myDF2[[da.last.grep]],labels=da.levels,ordered = F)

            colnames(myDF2)[end.num]=paste(label.explor.r.my.factor.list[[1]][last.last.temp],"_F",sep="")
            attr(myDF2[[end.num]],"label")=attr(myDF2[[da.last.grep]],"label")


            #### My factor variable list ####
            label.explor.r.my.factor.list.f=c(label.explor.r.my.factor.list.f,paste(label.explor.r.my.factor.list[[1]][last.last.temp],"_F",sep=""))
            if(length(attr(myDF2[[da.last.grep]],"label",exact=T))>0){
              label.explor.r.my.factor.list.f.labels=c(label.explor.r.my.factor.list.f.labels,attr(myDF2[[da.last.grep]],"label",exact=T))
            }else{
              label.explor.r.my.factor.list.f.labels=c(label.explor.r.my.factor.list.f.labels,names(myDF2)[da.last.grep])
            }
          }}else{
            print(paste(label.explor.r.my.factor.list[[1]][last.last.temp],"did not have labels to make a factor."))
          }}}
    last.last.temp=last.last.temp+1
  }

  #### Make Ordinal Factors ####
  #print(label.explor.r.my.ordered.list)
  ordinal.temp=1
  while(ordinal.temp<{length(label.explor.r.my.ordered.list[[1]])+1}){

    ordered.grep=grep(paste("^",label.explor.r.my.ordered.list[[1]][ordinal.temp],"$",sep=""),names(myDF2))
    ordered.levels=label.explor.r.my.ordered.list[[2]][[ordinal.temp]]
    ordered.flag=grep(paste("^",label.explor.r.my.ordered.list[[1]][ordinal.temp],"$",sep=""),names(flag))

    if(length(ordered.flag)>0){
      ordered.grep2=grep(paste("^",label.explor.r.my.ordered.list[[1]][ordinal.temp],"$",sep=""),names(myDF))
      ordered.table=table(myDF[[ordered.grep2]])
      ordered.labs=as.numeric(rownames(ordered.table))
      ordered.table2=table(myDF2[[ordered.grep]])
      ordered.labs2=as.numeric(rownames(ordered.table2))
      ordered.diff=setdiff(ordered.labs,ordered.labs2)

      ordered.caution=1
      while(ordered.caution<{length(ordered.diff)+1}){
        ordered.remove=grep(ordered.diff,ordered.labs)
        ordered.levels=ordered.levels[-ordered.remove]
        ordered.caution=ordered.caution+1
      }
    }

    #Check for missing
    ordered.levels2=tolower(ordered.levels)
    ordered.missing=grep("missing",ordered.levels2)
    if(length(ordered.missing)>0){
      ordered.levels=ordered.levels[-ordered.missing]
    }


    ordered.end.num=length(colnames(myDF2))+1
    myDF2[ordered.end.num]=factor(myDF2[[ordered.grep]],labels=ordered.levels,ordered = T)

    colnames(myDF2)[ordered.end.num]=paste(label.explor.r.my.ordered.list[[1]][ordinal.temp],"_O",sep="")
    attr(myDF2[[ordered.end.num]],"label")=attr(myDF2[[ordered.grep]],"label")


    #### My factor variable list ####
    label.explor.r.my.ordered.list.f=c(label.explor.r.my.ordered.list.f,paste(label.explor.r.my.ordered.list[[1]][ordinal.temp],"_O",sep=""))
    if(length(attr(myDF2[[ordered.grep]],"label",exact=T))>0){
      label.explor.r.my.ordered.list.f.labels=c(label.explor.r.my.ordered.list.f.labels,attr(myDF2[[ordered.grep]],"label",exact=T))
    }else{
      label.explor.r.my.ordered.list.f.labels=c(label.explor.r.my.ordered.list.f.labels,names(myDF2)[ordered.grep])
    }

    ordinal.temp=ordinal.temp+1
  }
  #### Bug Check ID
  if(length(grep(caseid,label.explor.r.my.factor.list.f.labels))>0){
    this.thang2=grep(caseid,label.explor.r.my.factor.list.f.labels)
    label.explor.r.my.factor.list.f.labels=label.explor.r.my.factor.list.f.labels[-this.thang2]
  }


  #### Set label names, set to variable name if none provided ####
  names.temp=1
  while(names.temp<{length(label.explor.r.my.selected)+1}){
    names.grep=grep(paste("^",label.explor.r.my.selected[names.temp],"$",sep=""),names(myDF2))
    if(!is.null(eval(parse(text=paste("nameInput",names.grep,sep=""))))){
      #print(eval(parse(text=paste("nameInput",names.grep,sep=""))))
      attr(myDF2[[names.grep]],"label")=eval(parse(text=paste("nameInput",names.grep,sep="")))
      label.explor.r.my.selected.names=c(label.explor.r.my.selected.names,attr(myDF2[[names.grep]],"label",exact = T))
    }else if(length(attr(myDF2[[names.grep]],"label",exact=T))>0){
      label.explor.r.my.selected.names=c(label.explor.r.my.selected.names,attr(myDF2[[names.grep]],"label",exact = T))
    }else{
      label.explor.r.my.selected.names=c(label.explor.r.my.selected.names,names(myDF2)[names.grep])
    }
    names.temp=names.temp+1
  }


  #### Remove caseid if selected in my.selected ####
  if(length(grep(caseid,label.explor.r.my.selected)>0)){
    caseidnumber2=grep(caseid,label.explor.r.my.selected)
    label.explor.r.my.selected=label.explor.r.my.selected[-caseidnumber2]
    label.explor.r.my.selected.names=label.explor.r.my.selected.names[-caseidnumber2]
  }

  #### Remove lines where all entries are NA ####
  myDF2=myDF2[rowSums(is.na(myDF2))!=ncol(myDF2),]

  #### Set dataframe attributes ####
  attr(myDF2,"selected")=label.explor.r.my.selected
  attr(myDF2,"selected.labels")=label.explor.r.my.selected.names
  attr(myDF2,"factors")=label.explor.r.my.factor.list.f
  attr(myDF2,"ordered")=label.explor.r.my.ordered.list.f
  attr(myDF2,"ordered.labels")=label.explor.r.my.ordered.list.f.labels
  attr(myDF2,"factors.labels")=label.explor.r.my.factor.list.f.labels
  attr(myDF2,"caseid")=caseid

  #### Remove all global variables created ####
  rm(label.explor.r.my.factor.list,envir = .GlobalEnv)
  rm(label.explor.r.my.selected,envir = .GlobalEnv)
  rm(label.explor.r.myList,envir = .GlobalEnv)
  rm(list=label.explor.r.temp.list,envir=.GlobalEnv)

  #savehistory(file=file.choose())
  return(myDF2)

  }
