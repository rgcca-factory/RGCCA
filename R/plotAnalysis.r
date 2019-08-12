plotAnalysis=function(listFinale,output="rv",fileName=NULL,ylim=c(0.8,1),block="all",barType="sd",namePlot=NULL,width=480,height=480)
{ #output : "rv", "pct" ou "a"
  #barType="sd" or "stdErr"
  if(is.null(namePlot)){namePlot=output}
  graphics.off()
  if(!is.null(fileName)){png(paste(fileName,".png",sep=""),width=width,height=height)}
  nameData= names(listFinale)
  abscisse=as.numeric(substr(nameData,5,7));names(abscisse)=nameData
  abscisse=1:length(listFinale[[1]])
  pas=1 
  par(las=1)
  J=length(listFinale[[1]][[1]][[1]]) #nblock
  if(block=="all"){ split.screen(c(2,2));toPlot=1:J}else{toPlot=block:block}
  # print(toPlot)
  namesMethod=names(listFinale[[1]])
  #colMethod=rainbow(5)[1:length(namesMethod)]
  colMethod=c("cornflowerblue","chocolate1","chartreuse3","red","blueviolet")[1:length(namesMethod)]
  nMeth=0:length(namesMethod)
  names(colMethod)=names(nMeth)=namesMethod
  for(j in toPlot)
  {
    if(block=="all"){screen(j)}
    par(mar=c(5, 4, 4, 2) + 0.1)
    par(mgp=c(3,1,0))
    plot(NULL,main=paste(namePlot,": Block",j),xlim=c(0,length(namesMethod)-1),ylim=ylim,xlab="Methods",ylab="Correlation",bty="n")
    axis(side = 1,col="grey",line=0)
    axis(side = 2,col="grey",line=0)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = 
           "#e9ebec",border="#e9ebec")
    # grid(nx = NULL, ny = NULL, col = "white", lty = "dotted",  lwd = par("lwd"), equilogs = TRUE)

      for(rg in namesMethod)
      {
        result=sapply(listFinale,function(x){return(x[[rg]][[output]][j])})
        moyenne=mean(result)
        if(!barType %in% c("sd","stderr")){ecartType=0}
        if(barType=="sd"){ecartType=sd(result)}
        if(barType=="stderr"){ecartType=sd(result)/sqrt(length(result))}
        
        points(pas*nMeth[rg],moyenne,pch=16,col=colMethod[rg])
        segments(pas*nMeth[rg],moyenne-ecartType,pas*nMeth[rg],moyenne+ecartType,col=colMethod[rg])
      }
  }
  if(block=="all")
  {
    screen(4)
    legend("center",legend=namesMethod,fill=colMethod,box.lwd=0)
  }
  if(is.numeric(block))
  {
    legend("bottomleft",legend=namesMethod,fill=colMethod,box.lwd=0)
  }
  if(!is.null(fileName)){dev.off()}
}
