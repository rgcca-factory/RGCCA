#' plots a list_rgcca object 
#' 
#' @param x Result of MIRGCCA function (see\code{\link[RGCCA]{MIRGCCA}} ) or a list containing 
#' \itemize{\item{rgcca0}{ RGCCA results for the reference dataset}
#' \item{data}{ list of imputed data obtained}
#'  \item{rgccaList}{ list of RGCCA obtained}}  
#' @inheritParams plot.rgcca
#' @examples
#' set.seed(42);X1=matrix(rnorm(500),100,5);
#' set.seed(22);X2=matrix(rnorm(400),100,4);
#' set.seed(2);X3=matrix(rnorm(700),100,7);
#'rownames(X1)=rownames(X2)=rownames(X3)=paste("S",1:100,sep="")
#'colnames(X1)=paste("A",1:5)
#'colnames(X2)=paste("B",1:4)
#'colnames(X3)=paste("C",1:7)
#'X1[1,]=NA
#'X2[7,1]=NA
#'X2[5,1]=NA
#'X3[3,1:2]=NA
#'A=list(X1,X2,X3)
#'res=MIRGCCA(A,k=3,ni=5,scale=TRUE,scale_block=TRUE,tau=rep(1,3))
#'plot(res,type="ind")

#' @param errorbar ("CImean","CIscores","sd")
#' @importFrom gridExtra grid.arrange
#' @export
plot.list_rgcca=function(x,type="ind",resp=rep(1, NROW(x$Y[[1]])),block=1,comp=1:2,remove_var=FALSE,text_var=TRUE,text_ind=TRUE,response_name= "Response",no_overlap=FALSE,title=NULL,n_mark=100,collapse=FALSE,cex=1,cex_sub=10,cex_main=14,cex_lab=12,colors=NULL,errorbar="ci",...)
{

    if(length(comp)==1){comp=rep(comp,2)}
    compx=comp[1]
    compy=comp[2]
    if(length(block)==1){block=rep(block,2)}
    i_block=block[1]
    i_block_y=block[2]
    
    lower_band <- NULL -> upper_band
    rgcca_res=x$rgcca0
    list_rgcca=x$rgccaList
    nRgcca=length(list_rgcca)
    match.arg(errorbar,c("cim","ci","sd","stderr"))
    match.arg(type,c("ind","var","cor","weight"))
     n=dim(rgcca_res$Y[[i_block]])[1]
         colors=c(rainbow(10),rainbow(10,s=0.7),rainbow(10,v=0.7),rainbow(10,s=0.5),rainbow(10,v=0.5),rainbow(max(n-50,0),s=0.3))
  if(type=="ind")
  {
      if(is.null(title)){title=paste0(names(rgcca_res$call$blocks)[i_block],": Sample space")}
      
      resp=1:n
      p1<-plot_ind(rgcca_res,resp=resp, i_block=i_block,i_block_y = i_block_y,compx=compx,compy=compy,legend=FALSE,colors=colors[1:n],cex=cex,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab) 
      colt=c()
      for(i in 1:(length(list_rgcca)))
      { 
          df1 <- get_comp(
              rgcca_res = list_rgcca[[i]],
              compx = compx,
              compy = compy,
              i_block_x = i_block,
              i_block_y = i_block_y,
              predicted = NULL
          )
          if(cor(df1[,1], rgcca_res$Y[[i_block]][,compx])<0){df1[,1]=-df1[,1]}
          if(cor(df1[,2], rgcca_res$Y[[i_block_y]][,compy])<0){df1[,2]=-df1[,2]}
          
          colnames(df1)=c("axis1","axis2")
          if(dim(df1)[1]!=length(resp)){stop_rgcca("two rgcca have two different numbers of subjects")}
          if(all.equal(rownames(df1),rownames(rgcca_res$Y[[i_block]]))!=TRUE){stop_rgcca("not same names in rgcca")}
          
          if(i==1){dft=df1;}else{    dft<-rbind(dft,df1)}
          colt=c(colt,colors[1:n])
          
      }                
      
      p2<- p1+geom_point(data=dft,aes(x=dft[,1],y=dft[,2]),colour=colt,size=0.8)
      #plot(p2)
      
  }
     
  if(type=="var")
     {
        if(is.null(title)){title=paste0(names(rgcca_res$call$blocks)[i_block],": Variable correlations")}
        
          colt=c()
         nvar=dim(rgcca_res$a[[i_block]])[1]
         resp=1
         p1 <- plot_var_2D(rgcca_res,resp=resp,i_block=i_block,compx=compx,compy=compy,cex_sub=cex_sub,cex_main=cex_main,cex_lab=cex_lab,remove_var=remove_var,text=text_var,no_overlap=no_overlap,title=title,n_mark = n_mark,collapse=collapse,colors=colors[1:nvar])
         
        
         for(i in 1:(length(list_rgcca)))
         { 
             if(cor(list_rgcca[[i]]$Y[[i_block]][,compx], rgcca_res$Y[[i_block]][,compx])<0){list_rgcca[[i]]$a[[i_block]][,compx]=-list_rgcca[[i]]$a[[i_block]][,compx]}
             if(cor(list_rgcca[[i]]$Y[[i_block_y]][,compy], rgcca_res$Y[[i_block_y]][,compy])<0){list_rgcca[[i]]$a[[i_block_y]][,compy]=-list_rgcca[[i]]$a[[i_block_y]][,compy]}
             
             df1 <- get_ctr2(
                 rgcca_res = list_rgcca[[i]],
                 compx = compx,
                 compy = compy,
                 i_block = i_block,
                 type = "cor",
                 n_mark = n_mark,
                 collapse = collapse,
                 remove_var = remove_var,
                 resp=resp
             )
             colnames(df1)=c("axis1","axis2")

             if(dim(df1)[1]!=dim(rgcca_res$a[[i_block]])[1]){stop_rgcca("two rgcca have two different numbers of subjects")}
             if(all.equal(rownames(df1),rownames(rgcca_res$a[[i_block]]))!=TRUE){stop_rgcca("not same names in rgcca")}
             
             if(i==1){dft=df1;}else{    dft<-rbind(dft,df1)}
             colt=c(colt,colors[1:nvar])
             
         }                
         
         p2<- p1+geom_point(data=dft,aes(x=dft[,1],y=dft[,2]),colour=colt,size=0.8)
         #plot(p2)
         
     
  }
  if(type=="cor")
  {
      if(is.null(title)){title=paste0(names(rgcca_res$call$blocks)[i_block],": Variable correlations")}
      
           
      attributes=colnames(rgcca_res$A[[i_block]])
     # list_rgcca_sup_a=lapply(1:length(list_rgcca),function(i)
     # {
      list_rgcca2=list_rgcca
      list_rgcca2[[length(list_rgcca2)+1]]=rgcca_res
      list_rgcca_sup_a=list()
     # print(attributes)
      for(i in 1:length(list_rgcca))
      { #print(colnames(list_rgcca[[i]]$call$blocks[[i_block]]))
          if(cor(list_rgcca[[i]]$Y[[i_block]][,compx], rgcca_res$Y[[i_block]][,compx])<0)
          {
                list_rgcca[[i]]$a[[i_block]][,compx]=-list_rgcca[[i]]$a[[i_block]][,compx];
                list_rgcca[[i]]$Y[[i_block]][,compx]=-list_rgcca[[i]]$Y[[i_block]][,compx]
           }
           
          res= sapply(1:length(attributes),function(att)
          {
              return(cor(list_rgcca[[i]]$Y[[i_block]][,compx],list_rgcca[[i]]$call$blocks[[i_block]][,attributes[att]]))
          })
          names(res)=attributes 
          list_rgcca_sup_a[[i]]=res
      }
      res0= sapply(1:length(attributes),function(att)
      {
          return(cor(rgcca_res$Y[[i_block]][,compx],rgcca_res$call$blocks[[i_block]][,attributes[att]]))
      })
      df1=Reduce(cbind, list_rgcca_sup_a)
      df_ordered=df1[order(abs(res0),decreasing=TRUE),]
      #  return(res)
      #}    )
      
  }
  if(type=="weight")
  {
      if(is.null(title)){title=paste0(names(rgcca_res$call$blocks)[i_block],": Variable weight")}
      
      for(i in 1:length(list_rgcca))
      {
          if(cor(list_rgcca[[i]]$Y[[i_block]][,compx], rgcca_res$Y[[i_block]][,compx])<0){list_rgcca[[i]]$a[[i_block]][,compx]=-list_rgcca[[i]]$a[[i_block]][,compx];list_rgcca[[i]]$Y[[i_block]][,compx]=-list_rgcca[[i]]$Y[[i_block]][,compx]}
      }
   
      list_rgcca_sup_a=lapply(list_rgcca, function(v){return(v$a[[i_block]][,compx])})
      list_rgcca_sup_a[[length(list_rgcca_sup_a)+1]]=rgcca_res$a[[i_block]][,compx]
  
      df1=Reduce(cbind, list_rgcca_sup_a)
      df_ordered=df1[order(abs(rgcca_res$a[[i_block]][,compx]),decreasing=TRUE),]
    }
  if(type%in% c("cor","weight"))
  {
            p1=plot_var_1D(rgcca_res, comp = compx, n_mark = n_mark,
                     i_block = i_block, type = type, collapse = collapse,
                     title = title, colors = colors,...)
     
      statT=qt(0.975,df=nRgcca-1)
      if(errorbar=="cim")
      {
          lowerBand=apply(df_ordered,1,function(x){return(mean(x)-statT*sd(x)/sqrt(nRgcca))})
          upperBand=apply(df_ordered,1,function(x){return(mean(x)+statT*sd(x)/sqrt(nRgcca))})
      }
      if(errorbar=="ci")
      {
          lowerBand=apply(df_ordered,1,function(x){return(mean(x)-statT*sd(x))})
          upperBand=apply(df_ordered,1,function(x){return(mean(x)+statT*sd(x))})
      }
      if(errorbar=="sd")
      {
          lowerBand=apply(df_ordered,1,function(x){return(mean(x)-sd(x))})
          upperBand=apply(df_ordered,1,function(x){return(mean(x)+sd(x))})
      }
      if(errorbar=="stderr")
      {
          lowerBand=apply(df_ordered,1,function(x){return(mean(x)-sd(x)/sqrt(nRgcca))})
          upperBand=apply(df_ordered,1,function(x){return(mean(x)+sd(x)/sqrt(nRgcca))})
      }
      df_ordered=as.data.frame(df_ordered)
      df_ordered[,"lower_band"]=lowerBand
      df_ordered[,"upper_band"]=upperBand
      p2<- p1+geom_errorbar(data=df_ordered,aes(x=rev(1:dim(df_ordered)[1]),ymin = lower_band, ymax = upper_band,width=0.5))
      
    }

      
      
  
     
       return(p2)
}    
  