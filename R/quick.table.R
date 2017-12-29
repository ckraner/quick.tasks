#' Table Printer
#'
#' Helper function to turn tables into html tables
#'
#' @param my.table Final table ready to be turned into html
#' @param type Type (lm, glm, manova, ordinal)
#' @param test.stat Test stat used in MANOVA
#' @param print.type "full" for full html, "part" for just <div>
#' @param the.caption Caption for top of table
#' @param the.footer
#' @return List of lists with 3 pieces: null model, pre models, and full models
#' @keywords Explore
quick.table=function(my.table,
                     type,
                     test.stat="Pillai",
                     print.type="full",
                     the.caption=NULL,
                     the.footer=NA,
                     abbrev.length=ab.len,
                     SS.type=2,
                     new.rownames.int=NULL,
                     new.rownames.treat=NULL,
                     swap.na=NULL,
                     round.num=2,
                     col.names=my.colnames,
                     print.now=T,
                     show.footer=T,
                     make.red=NULL,
                     make.black=NULL){

  #### Inits ####

  if(type=="ord"){
    ab.len=30
    library(ordinal)
  }else{
    ab.len=15
  }

  if(type=="manova" | type=="stats::manova"){
    attr(my.table,"quick.test.stat")=test.stat
  }
  attr(my.table,"quick.print.type")=print.type
  attr(my.table,"quick.abbrev.length")=abbrev.length
  attr(my.table,"quick.round")=round.num
  attr(my.table,"quick.type")=type
  attr(my.table,"quick.footer")=the.footer
  attr(my.table,"quick.SS.type")=SS.type
  attr(my.table,"class")=c(attr(my.table,"class"),"quick.table")


  if(type=="manova" | type=="stats::manova"){
    round.rows=c(2,3,4)
    p.val.row=8
    my.colnames=c("Variable",paste(test.stat,"<br />Statistic"),"F-Value",
                  paste("Type ",ifelse(SS.type==2,"II","III"),"<br />Sums of<br />Squares"),
                  "dF","Mult<br />dF","Resid<br />dF","Pr(>F)")
  }
  attr(my.table,"quick.col.names")=col.names


  #### Find Intercept, Treatment, Total Change Locations ####
  int.loc=grep("Intercept",my.table[[1]])
  treat.loc=grep("Treatment Change",my.table[[1]])
  total.loc=grep("Total Change",my.table[[1]])
  end.loc=grep("^Total$",my.table[[1]])
  my.table2=my.table
  #### Swap value for NA (i.e. remove it) ####
  if(!is.null(swap.na)){
    na.vals=strsplit(swap.na)[[1]]
  }

  #### Replace rownames ####
  if(!is.null(new.rownames.int)){
    #### Do later
  }

  #### Round ####
  for(i in 1:length(round.rows)){
    my.table2[[round.rows[i]]]=round(as.numeric(my.table2[[round.rows[i]]]),digits =round.num)
  }

  #### P-val ####
  my.table2=quick.p.val(my.table2,p.val.row)


  #### Change NA to &nbsp; ####
  my.table2=replace(my.table2,is.na(my.table2),"&nbsp;")

  #### Add style and basic tag structure####
  attr(my.table,"quick.doctype")=paste("<!DOCTYPE html --- Created with quick.tasks by Christopher Kraner>")
  attr(my.table,"quick.full.start")=paste("<html><head><style>table{border: 1px solid black;border-collapse: collapse;}",
                                          "th{padding: 15px;}td {padding: 5px;}#red {border: 2px solid red;}",
                                          "#black {border: 2px solid black;}#change {border-top: 1px solid black;text-align: left;}",
                                          "tr:hover {background-color: #f5f5f5;}#int {border-top: 1px solid black;}#col {border-bottom: 1px solid black;}</style></head>",sep="")
  attr(my.table,"quick.part.start")=paste("<div style=\"overflow-x:auto;\"><table style=\"width:100%\">")
  attr(my.table,"quick.caption")=ifelse(is.null(the.caption),NA,paste("<caption>",the.caption,"</caption>"))
  attr(my.table,"quick.part.end")=paste("</table></div>")
  attr(my.table,"quick.full.end")=paste("</body></html>")

  ##### Start table ####
  if(print.type=="full"){
    my.html.table=paste(attr(my.table,"quick.doctype"),attr(my.table,"quick.full.start"),attr(my.table,"quick.part.start"))
  }else{
    my.html.table=paste(attr(my.table,"quick.part.start"),attr(my.table,"quick.doctype"))
  }

  #### Put in caption
  if(!is.null(the.caption)){
    my.html.table=paste(my.html.table,attr(my.table,"quick.caption"))
  }

  #### Put in Column Headings
  if(type!="glm"){
    col.headings="<tr id=\"col\", align=\"center\">"
    for(i in 1:length(my.colnames)){
      col.headings=paste(col.headings,"<th>",my.colnames[i],"</th>")
    }
    col.headings=paste(col.headings,"</tr>")
  }else{
    stop("Sorry, not to GLM yet.")
  }
  my.html.table=paste(my.html.table,col.headings)

  #### Put in Rows
  for(i in 1:end.loc){

    #### Variable name
    if(i==1 & my.table2[1,1]=="Intercept Change"){
      #### GLM stuff
      #### Make th, add id="change"
      if(i %in% make.red | i %in% make.black){
        if(i %in% make.red){
          my.line=paste("<tr id=\"red\"><th>",my.table2[i,1],"</th>")
        }else{
          my.line=paste("<tr id=\"black\"><th>",my.table2[i,1],"</th>")
        }
      }else{
      my.line=paste("<tr id=\"int\"><th>",my.table2[i,1],"</th>")
      }
    }else if(i==treat.loc | i==total.loc){
      if(i %in% make.red | i %in% make.black){
        if(i %in% make.red){
          my.line=paste("<tr id=\"red\"><td align=\"left\"><b>",my.table2[i,1],"</b></td>")
        }else{
          my.line=paste("<tr id=\"black\"><td align=\"left\"><b>",my.table2[i,1],"</b></td>")
        }
      }else{
      my.line=paste("<tr id=\"change\"><td align=\"left\"><b>",my.table2[i,1],"</b></td>")
      }
    }else if(i==1 & i==int.loc){
      if(i %in% make.red | i %in% make.black){
        if(i %in% make.red){
          my.line=paste("<tr id=\"red\"><td>",my.table2[i,1],"</td>")
        }else{
          my.line=paste("<tr id=\"black\"><td>",my.table2[i,1],"</td>")
        }
      }else{
      my.line=paste("<tr id=\"int\"><td>",my.table2[i,1],"</td>")
      }
    }else if(i>total.loc){
      if(i %in% make.red | i %in% make.black){
        if(i %in% make.red){
          my.line=paste("<tr id=\"red\"><td><b>",my.table2[i,1],"</b></td>")
        }else{
          my.line=paste("<tr id=\"black\"><td><b>",my.table2[i,1],"</b></td>")
        }
      }else{
        my.line=paste("<tr><td><b>",my.table2[i,1],"</b></td>")
      }
    }else{
      if(i %in% make.red | i %in% make.black){
        if(i %in% make.red){
          my.line=paste("<tr id=\"red\"><td>",my.table2[i,1],"</td>")
        }else{
          my.line=paste("<tr id=\"black\"><td>",my.table2[i,1],"</td>")
        }
      }else{
      my.line=paste("<tr><td>",my.table2[i,1],"</td>")
      }
    }

    #### Rest of row
    for(j in 2:p.val.row){
      my.line=paste(my.line,"<td>",my.table2[i,j],"</td>")
    }

    my.line=paste(my.line,"</tr>")

    my.html.table=paste(my.html.table,my.line)
    if(i==1){
      attr(my.table,"quick.rows")=my.line
    }else{
      attr(my.table,"quick.rows")=paste(attr(my.table,"quick.rows"),my.line)
    }
  }

  #### End table
  my.html.table=paste(my.html.table,"</table>")


  #### Put in custom bottom
  if(!is.na(the.footer) & show.footer){
    my.html.table=paste(my.html.table,"<p align=\"center\">",the.footer,"</p>")
  }
  #### Put in end
  my.html.table=paste(my.html.table,"</div>")
  if(print.type=="full"){
    my.html.table=paste(my.html.table,"</body></html>")
  }

  if(print.now){
    tempDir <- tempfile()
    dir.create(tempDir)
    htmlFile <- file.path(tempDir, "index.html")
    writeLines(my.html.table,htmlFile)
    viewer <- getOption("viewer")
    viewer(htmlFile)
  }

  return(my.table)
}

#' Table Check
#'
#' Checks to make sure that table and html table are the same before update
#'
#' @param q.tab quick.table
#' @return Logical of whether all matches
#' @keywords Explore
quick.table.check=function(q.tab){

  #### Turn HTML into something easily useable
  row.check=attr(q.tab,"quick.rows")
  row.split=strsplit(row.check,"</tr>")
  col.split=lapply(row.split,strsplit,"</td>")
  col.split=lapply(col.split[[1]],strsplit,"<td>")

  my.comp.table=NULL
  for(i in 1:length(col.split)){
    my.temp.row=unlist(col.split[[i]])
    my.temp.row=my.temp.row[my.temp.row != " "]
    if(length(grep("<b>",my.temp.row[1]))==0){
    my.temp.row=my.temp.row[-1]
    }else{
      my.temp.row[1]=strsplit(my.temp.row[1],"<b>")[[1]][2]
      my.temp.row[1]=strsplit(my.temp.row[1],"</b>")[[1]][1]
    }
    my.temp.row=trimws(my.temp.row)

    if(i==1){
      my.comp.table=my.temp.row
    }else{
      my.comp.table=rbind(my.comp.table,my.temp.row)
    }
  }

  #### Get round number and round table
  round.num=attr(q.tab,"quick.round")
  if(attr(q.tab,"quick.type")=="manova" | attr(q.tab,"quick.type")=="stats:manova"){
    q.tab[[2]]=round(as.numeric(q.tab[[2]]),digits=round.num)
    q.tab[[3]]=round(as.numeric(q.tab[[3]]),digits=round.num)
    q.tab[[4]]=round(as.numeric(q.tab[[4]]),digits=round.num)
  }

  #### P-val
  p.row=grep("p.val",colnames(q.tab))
  q.tab=quick.p.val(q.tab,p.row)
  #### Check against values in table
  my.map.table=NULL
  for(i in 2:dim(q.tab)[2]){
    if(i==2){
      my.map.table=map2(q.tab[[i]],my.comp.table[,i],quick.eq.check)
    }else{
      my.map.table=cbind(my.map.table,map2(q.tab[[i]],my.comp.table[,i],quick.eq.check))
    }
  }
  no.false=T
  for(i in 1:dim(my.map.table)[2]){
  if(length(my.map.table[which(my.map.table[,i]==F),])>0){
    no.false=F
  }
  }

  return(no.false)
}

#' Table Update
#'
#' Update unchanged table
#'
#' @param q.tab quick.table
#' @return Logical of whether all matches
#' @export
#' @keywords Explore
quick.table.update=function(q.tab,make.red=NULL,make.black=NULL,the.caption=my.caption,show.footer=T,new.rownames.int=NULL,
                            new.rownames.treat=NULL,swap.na=NULL,the.round=my.round,print.type="full",
                            print.now=T,do.return=F){
  my.check=quick.table.check(q.tab)
  if(my.check){
    my.caption=attr(q.tab,"quick.caption")
    my.round=attr(q.tab,"quick.round")
    my.type=attr(q.tab,"quick.type")
    my.ab.len=attr(q.tab,"abbrev.length")
    my.SS.type=attr(q.tab,"quick.SS.type")
    if(my.type=="manova" | my.type=="stats::manova"){
    my.test.stat=attr(q.tab,"quick.test.stat")
    }else{
      my.test.stat=NULL
    }
    new.q.tab=quick.table(q.tab,type,test.stat=my.test.stat,print.type=print.type,
               the.caption=the.caption,the.footer=the.footer,
               abbrev.length=my.ab.len,
               SS.type=my.SS.type,
               new.rownames.int=new.rownames.int,
               new.rownames.treat=new.rownames.treat,
               swap.na=swap.na,
               round.num=the.round,
               col.names=my.colnames,
               print.now=T,
               show.footer=show.footer,
               make.red=make.red,
               make.black=make.black)
    if(do.return){
    return(new.q.tab)
    }
  }else{
    stop("This table has changed.")
  }

}
