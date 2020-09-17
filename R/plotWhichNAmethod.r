#' Which missing method to choose ?
#' 
#' Plots the impact of increasing missing data on RGCCA
#' @param x A list resulting of whichNAmethod (see \link{whichNAmethod})
#' @param type ="rv": Can be also "a" for correlations between axes, "bm" for the percent of similar biomarkers, "rvComplete" if the RV is calculated only on complete dataset, or "rmse" for Root Mean Squares Error.
#' @param ylim =c(0.8,1) y limits
#' @param block ="all" or a number indicating the position of the chosen block in the initial list
#' @param bars ="sd" or "stderr". Indicates which error bar to build
#' @param main =NULL Name of the file
#' @param ylab label of y-axis
#' @param legend If TRUE, the legend is displayed
#' @param namesForLegend names of the used methods to be plotted on the legend
#' @param ... Further graphical parameters in plot
#' @examples 
#' set.seed(42);X1=matrix(rnorm(350),70,5);X2=matrix(rnorm(280),70,4);X1[1,1]=NA;X2[2,]=NA
#' colnames(X1)=paste("A",1:5);colnames(X2)=paste("B",1:4); 
#' rownames(X1)=rownames(X2)=paste0("S",1:70);A=list(X1,X2);
#' listResults=whichNAmethod(blocks=A,patternNA=c(0.1,0.2),
#' listMethods=c("mean","complete","nipals","knn4"))
#' plot(x=listResults,ylim=c(0,1),type="a")
#' @importFrom grDevices graphics.off 
#' @importFrom graphics plot.new
#' @export


plot.whichNAmethod=function(x,type="rv",ylim=NULL,block=length(x[[1]][[1]][[1]]),bars="sd",main=NULL,ylab="",legend=TRUE,namesForLegend=NULL,...)
{ #type : "rv", "pct" ou "a"
  #bars="sd" or "stderr"
    
    # TODO: par(new=TRUE)
  match.arg(bars,c("sd","stderr"))
  match.arg(type,c("rv","rvComplete","a","rmse","bm"))
  if(is.null(main)){main=type}
  #graphics.off()
 nameData= names(x)
  abscisse=as.numeric(substr(nameData,5,7));names(abscisse)=nameData
  abscisse=1:length(x[[1]])
  pas=1 
  par(las=1)
  J=length(x[[1]][[1]][[1]]) #nblock
 # close.screen(all.screens=TRUE)
  if(block=="all")
  { 
      par(mfrow=c(floor(sqrt(J)+1),floor(sqrt(J)+1)));toPlot=1:J
    #  split.screen(c(floor(sqrt(J)+1),floor(sqrt(J)+1)));toPlot=1:J
  }
  else
  {
      toPlot=block:block
  }
  # print(toPlot)

  namesMethod=names(x[[1]]) 
  if(is.null(namesForLegend))
  {
      namesForLegend=namesMethod
  }

  print(namesMethod)
  #colMethod=rainbow(5)[1:length(namesMethod)]
  colMethod=c("cornflowerblue","chocolate1","chartreuse3","red","blueviolet","darkturquoise","darkgoldenrod1","coral","bisque4","darkorchid1","deepskyblue1")[1:length(namesMethod)]
  nMeth=0:length(namesMethod)
  names(colMethod)=names(nMeth)=namesMethod
  for(j in toPlot)
  {
    #if(block=="all"){screen(j)}
    par(mar=c(5, 4, 4, 2) + 0.1)
    par(mgp=c(3,1,0))
 
    moyenne=rep(NA,length(namesMethod));names(moyenne)=namesMethod
    ecartType=rep(NA,length(namesMethod));names(ecartType)=namesMethod
    
    for(rg in namesMethod)
    {
      result=sapply(x,function(x){return(x[[rg]][[type]][[j]])})
      moyenne[rg]=mean(result)
      if(bars =="no"){ecartType=0}
      if(bars=="sd"){ecartType[rg]=sd(result)}
      if(bars=="stderr"){ecartType[rg]=sd(result)/sqrt(length(result))}
    } 
    if(is.null(ylim))
    { 
          minim=min(moyenne-ecartType)
          maxim=max(moyenne+ecartType)
        if(!is.na(minim))
        {
          Ylim=c(minim,maxim)
        }
        else{Ylim=ylim}
    }
    else
    {
        Ylim=ylim
    }
    plot(NULL,main=paste(main,": Block",j),xlim=c(0,length(namesMethod)-1),ylim=Ylim,xlab="Methods",ylab=ylab,bty="n",xaxt="n",...)
    axis(side = 1,col="grey",line=0)
    axis(side = 2,col="grey",line=0)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
           "#e9ebec",border="#e9ebec")
    # grid(nx = NULL, ny = NULL, col = "white", lty = "dotted",  lwd = par("lwd"), equilogs = TRUE)

       for(rg in namesMethod)
      {
        if(!(rg=="complete"&&type=="rv"))
        {
          points(pas*nMeth[rg],moyenne[rg],pch=16,col=colMethod[rg])
          segments(pas*nMeth[rg],moyenne[rg]-ecartType[rg],pas*nMeth[rg],moyenne[rg]+ecartType[rg],col=colMethod[rg])
          
        }
     }
  }
  if(legend)
  {
      if(block=="all")
      {
          #  screen(J+1)
          plot.new()
          par(cex=0.8)
          legend("center",legend=namesForLegend,fill=colMethod,box.lwd=0,,bty="n")
      }
      if(is.numeric(block))
      {
          par(cex=0.8)
          legend("bottomleft",legend=namesForLegend,fill=colMethod,box.lwd=0,bty="n")
      }
      
      par(mfrow=c(1,1))
      par(mar=c(5,4,3,3))
      par(cex=1)  
  }

}
