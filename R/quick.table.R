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
                     the.footer=NULL,
                     abbrev.length=ab.len,
                     SS.type=2,
                     new.rownames.int=NULL,
                     new.rownames.treat=NULL,
                     swap.na=NULL,
                     round.num=2,
                     col.names=my.colnames,
                     print.now=T){

  #### Inits ####

  if(type=="ord"){
    ab.len=30
    library(ordinal)
  }else{
    ab.len=15
  }

  attr(my.table,"quick.print.type")=print.type
  attr(my.table,"quick.abbrev.length")=abbrev.length
  attr(my.table,"class")=c(attr(my.table,"class"),"quick.table")


  if(type=="manova" | type=="stats::manova"){
    round.rows=c(2,3,4)
    p.val.row=8
    my.colnames=c("Variable",paste(test.stat,"<br />Statistic"),"F-Value",
                  paste("Type ",ifelse(SS.type==2,"II","III"),"<br />Sums of<br />Squares"),
                  "dF","Mult<br />dF","Resid<br />dF","Pr(>F)")
  }

  #### Find Intercept, Treatment, Total Change Locations ####
  int.loc=grep("Intercept",my.table[[1]])
  treat.loc=grep("Treatment Change",my.table[[1]])
  total.loc=grep("Total Change",my.table[[1]])
  end.loc=grep("^Total$",my.table[[1]])

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
    my.table[[round.rows[i]]]=round(as.numeric(my.table[[round.rows[i]]]),digits =round.num)
  }

  #### P-val ####
  for(i in 1:end.loc){
    if(!is.na(my.table[i,p.val.row])){
      my.table[i,p.val.row]=round(as.numeric(my.table[i,p.val.row]),digits=3)
      if(my.table[i,p.val.row]==0){my.table[i,p.val.row]="<.001"}
      if(my.table[i,p.val.row]==1){my.table[i,p.val.row]=">.999"}
    }
  }

  #### Change NA to &nbsp; ####
  my.table=replace(my.table,is.na(my.table),"&nbsp;")

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
    if(i==1 & my.table[1,1]=="Intercept Change"){
      #### GLM stuff
      #### Make th, add id="change"
      my.line=paste("<tr id=\"int\"><th>",my.table[i,1],"</th>")
    }else if(i==treat.loc | i==total.loc){
      my.line=paste("<tr id=\"change\"><td align=\"left\"><b>",my.table[i,1],"</b></td>")
    }else if(i==1 & i==int.loc){
      my.line=paste("<tr id=\"int\"><td>",my.table[i,1],"</td>")
    }else if(i>total.loc){
      my.line=paste("<tr><td align=\"left\"><b>",my.table[i,1],"</b></td>")
    }else{
      my.line=paste("<tr><td>",my.table[i,1],"</td>")
    }

    #### Rest of row
    for(j in 2:p.val.row){
      my.line=paste(my.line,"<td>",my.table[i,j],"</td>")
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
  if(!is.null(the.footer)){
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
