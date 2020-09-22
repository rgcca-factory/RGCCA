#' plot.patternNA
#' 
#' Plot an object of class patternNA
#' @param x A get_patternNA object (see \code{\link[RGCCA]{get_patternNA}})
#' @param type "all" if all blocks should be plotted, ifelse an integer corresponding to the position of the block to be plotted in the initial list
#' @param legend if TRUE the legend is plotted
#' @param color color of the non-missing data. If rainbow, depends of the values of the non-missing data
#' @param outlierVisible if FALSE, the outliers will be -2* standard deviations if negative, 2 standard deviations if positive
#' @param scale if TRUE all the variables are scaled before graphical representation
#' @param ... Further graphical parameters (see \code{\link[graphics]{plot}})
#' @export
#' @examples
#' X1=matrix(rnorm(150),30,5)
#' X2=matrix(rnorm(150),30,5)
#'X3=matrix(rnorm(150),30,5)
#'X1[1:10,1]=NA
#'X1[2,2]=NA
#'X1[11,]=NA
#'X2[1:10,1]=NA
#'colnames(X1)=paste0("A",1:5)
#'colnames(X2)=paste0("B",1:5)
#'colnames(X3)=paste0("C",1:5)
#'A=list(bloc1=X1,bloc2=X2,bloc3=X3)
#'p=get_patternNA(A)
#'plot(p)
#'@seealso \link{get_patternNA},\link{whichNAmethod}

plot.patternNA=function(x,type="all",color="springgreen4",legend=FALSE,scale=TRUE,outlierVisible=FALSE,...)
{
    completeSubjectByBlock <- NULL
    blocks=x$blocks
    if(type=="all")
    {
        nvar=sapply(blocks,NCOL)
        par(mfrow=c(1,ifelse(legend,length(blocks)+1,length(blocks))))
        par(mar=c(2,1,4,1))
        for(i in 1:length(blocks))
        {
            mat=apply(blocks[[i]],2,rev)
            mat2=apply(mat,2,function(x) scale(x,scale=scale))
            
            minimum=-2*max(apply(mat2,2,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
            maximum=2*max(apply(mat2,2,function(x) sd(x,na.rm=TRUE)),na.rm=TRUE)
            if(!outlierVisible)
            {
                mat2[mat2< minimum]=minimum
                mat2[mat2>maximum]=maximum
            }
          #  par(fg="black")
            plot(NULL,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty="n",main=paste0(names(blocks)[i],"\n",nvar[i]," var.\n",sum(x$completeSubjectByBlock[,i]), "/",NROW(blocks[[i]])," complete ind."),...)
            rect(0,0,1,1,col="black")
          
            if(color=="rainbow")
            { 
                image(z=t(mat2),breaks=seq(minimum,maximum,length.out=51),xaxt="n",yaxt="n",col=rainbow(50,start=0,end=0.5),add=TRUE,xlim=c(0,1),ylim=c(0,1))
            }
            else
            {
                image(t(mat2),xaxt="n",yaxt="n",col=color,add=TRUE,xlim=c(0,1),ylim=c(0,1))
            }
           
        }
    }
    if(type %in% names(blocks))
    {
        nvar=NCOL(blocks[[type]])
        par(mfrow=c(1,1))
        mat=apply(blocks[[type]],2,rev)
      #  par(bg="black")
        image(t(mat),main=paste0(names(blocks)[type],"\n(",nvar," var.,",sum(completeSubjectByBlock[,type]), "/",NROW(blocks[[type]])," complete ind.)"),xaxt="n",yaxt="n",col=c("light blue","black"))
        par(bg="white")
    }
    if(legend & color=="rainbow")
    {
        plot(NULL,xlim=c(0,1),ylim=c(0,50),bty="n",xaxt="n",yaxt="n")
        lapply(1:50,function(i){points(0.5,i,col=rainbow(50,start=0,end=0.5)[i],pch=15)})
        text(0.6, 1, "min");   text(0.6, 50, "max")
        
    }
}