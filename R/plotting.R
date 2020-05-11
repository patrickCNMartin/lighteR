################################################################################
############################ NPQ analysis ######################################
################################################################################

## Plotting Functions

plotScoresOverTimePerPlate <- function(data,file=NULL,combine=FALSE,col=NULL){
   ## testing classes

   if(is(data)[1]=="data.frame"){
      data<-list("Unspecified"=data)
   }

   # Filter single points for plotting
   idx<-sapply(data, singleElememntPos)
   data<-data[!idx]

   ## Initalise plot window

   # set limits
   xlim <- c(1,nrow(data[[1]]))
   ylim <- c(min(sapply(data, function(x){min(as.numeric(x[,match.column.Type(x)]))})),
             max(sapply(data, function(x){max(as.numeric(x[,match.column.Type(x)]))})))

   # set cols
   # could be improved  for higher flexibility
   # thisis pretty static input
   if(is.null(col)){
      cols <- viridis(length(data))
   }
   # file name set up

   if(!combine){
       if(!is.null(file)){
         pdf("test.pdf", width=10, height=4)
       }

       par(xpd=NA)
       par(mar=c(7,7,4,12))
       plot(0,type="n", xlim=xlim, ylim=ylim,xlab="",ylab="",axes=F)
       title(xlab="Time Point", line=5.5)
       title(ylab="Associated Score", line=5.5)
       axis(1,at=round(seq(xlim[1],xlim[2],length.out=10)), labels=paste("Time Point",round(seq(xlim[1],xlim[2],length.out=10))),las=2, cex.axis=0.75)
       axis(2,at=seq(ylim[1],ylim[2],length.out=5), labels=round(seq(ylim[1],ylim[2],length.out=5),2),las=2)
       legend(x=xlim[2]+5,y=ylim[2],legend=names(data),fill=cols,bty="n",cex=1)

       for(i in seq_along(data)){
           lines(seq(xlim[1],xlim[2]),as.numeric(data[[i]][,match.column.Type(data[[i]])]),col=cols[i],lwd=1.6)
       }
       if(!is.null(file)){
         dev.off()
       }
   }



}


## ploting plates per zone once extracted and cleaned

plotScoresOverTimePerZone <- function(data,model=NULL,file=NULL,split=c("none","measure","ID"),
                             quant=c(0.025,0.975),time=c(40,80),combine=TRUE,col=NULL){
     ## alright we are going to assume that they are already cleaned and split

     if(is.null(col)){
         cols <- c("#666666","#000061","black","#cc0000")
     }

     if(combine){
          ## this part might need to be changed if I do split plots

         if(class(data)=="byMeasure"){
           .plotMeasure(data,quant,cols,file, split,model,time)

         }else if(class(data)=="byIDCombined"){
            .plotbyID(data,quant,cols,file, split,model,time)
         } else if(class(data)=="byID"){
            .plotbyID(data,quant,cols,file, split,model,time)
         }


     }


}



##### plot by measure
.plotMeasure <- function(data,quant,cols,file,split,model,time=c(40,80)){
  for(i in seq_along(data)){
    local <-data[[i]][,2:ncol(data[[i]])]
    xlim<-c(1,ncol(local))
    ylim<-c(min(apply(local,1,min)),max(apply(local,1,max)))

    UpperBound<-apply(local,2,quantile,quant[2])
    LowerBound<-apply(local,2,quantile,quant[1])
    med <- apply(local,2,median)

    plot(0,type="n",ylab=names(data)[i],xlab="Time points",main=paste("Summary Score Distribution",names(data)[i]), xlim=xlim, ylim=ylim)

    for(j in seq_len(nrow(local))){

       lines(seq_along(as.vector(as.matrix(local[j,]))),as.vector(as.matrix(local[j,])), col=cols[1])
     }

     lines(seq_along(UpperBound),UpperBound, col=cols[2])
     lines(seq_along(LowerBound),LowerBound, col=cols[2])
     lines(seq_along(med),med, col=cols[3],lwd=2)
   }
}


## plot by ID

## might do a multi pdf thing
## feel like that would be cleaner







## plot by ID batch


.plotbyID <- function(data, quant, cols,file, split,model,time=c(40,80)){

if(grepl("none",split[1], ignore.case=TRUE) &!is.null(file)){
   pdf(paste0(file,"All_ZoneData.pdf"),width=15, height=15)
   par(mfrow=c(3,3))
}
  for(i in seq_along(data)){
    if(grepl("measure",split[1], ignore.case=TRUE) &!is.null(file)){
        pdf(paste0(file,names(data)[i],"_ZoneData.pdf"),width=15, height=15)
        if(length(data[[i]])>8)par(mfrow=c(3,3))
    }
    for(k in seq_along(data[[i]])){
      if(grepl("byID",split[1], ignore.case=TRUE)& !is.null(file)){
         pdf(paste0(file,names(data[[i]])[k],"_ZoneData_.pdf"),width=15, height=15)
         if(length(data[[i]])>8)par(mfrow=c(3,3))
      }
       local<-data[[i]][[k]][,6:ncol(data[[i]][[k]])]
       xlim<-c(1,ncol(local))
       ylim<-c(min(apply(local,1,min)),max(apply(local,1,max)))

       UpperBound<-apply(local,2,quantile,quant[2])
       LowerBound<-apply(local,2,quantile,quant[1])
       med <- apply(local,2,median)
       mainTit<-paste("ID combination : ", names(data[[i]])[k])
       mainTit<-gsub("\\."," ", mainTit)
       plot(0,type="n",ylab=names(data)[i],xlab="Time points",main="", xlim=xlim, ylim=ylim)
       title(main=mainTit, cex.main=1.1)

       for(j in seq_len(nrow(local))){

         points(seq_along(as.vector(as.matrix(local[j,]))),as.vector(as.matrix(local[j,])), col=cols[1])
       }

       lines(seq_along(UpperBound),UpperBound, col=cols[2])
       lines(seq_along(LowerBound),LowerBound, col=cols[2])
       lines(seq_along(med),med, col=cols[3],lwd=2)
       if(!is.null(model) & class(model)=="individualFit"){
         highseq<-seq(1,time[1])
         lowseq<-seq(time[1]+1,time[2])

         tmpFit <- .modelFitPlot(model[[i]][[k]],time,fit=class(model))

         for(j in seq_along(tmpFit[[1]])){
           ## plotting fitted lines
           lines(highseq,tmpFit$High[[j]], lwd=2, col=cols[4])
           lines(lowseq,tmpFit$Low[[j]], lwd=2, col=cols[4])
         }
       } else if(!is.null(model) & class(model)=="combinedFit"){
         highseq<-seq(1,time[1])
         lowseq<-seq(time[1]+1,time[2])
         tmpFit <- .modelFitPlot(model[[i]][[k]],time,fit=class(model))
         ## plotting fitted lines
         lines(highseq,tmpFit$High, lwd=2, col=cols[4])
         lines(lowseq,tmpFit$Low, lwd=2, col=cols[4])


       } else if(!is.null(model) & class(model)=="medianFit"){
         highseq<-seq(1,time[1])
         lowseq<-seq(time[1]+1,time[2])

         tmpFit <- .modelFitPlot(model[[i]][[k]],time,fit=class(model))
         sdHigh <- model[[i]][[k]][["sdHigh"]][[1]]
         sdLow <- model[[i]][[k]][["sdLow"]][[1]]

         for(j in seq_along(tmpFit[[1]])){
           ## plotting fitted lines
           lines(highseq,tmpFit$High[[j]], lwd=2, col=cols[4])

           lines(lowseq,tmpFit$Low[[j]], lwd=2, col=cols[4])

           arrows(highseq[round(seq(1,length(highseq),length.out=10))],
                 (tmpFit$High[[j]]-sdHigh)[round(seq(1,length(highseq),length.out=10))],
                 highseq[round(seq(1,length(highseq),length.out=10))],
                 (tmpFit$High[[j]]+sdHigh)[round(seq(1,length(highseq),length.out=10))],
                 length=0.05,angle=90, code=3,col=cols[4])

           arrows(lowseq[round(seq(1,length(lowseq),length.out=10))],
                 (tmpFit$Low[[j]]-sdLow)[round(seq(1,length(lowseq),length.out=10))],
                 lowseq[round(seq(1,length(lowseq),length.out=10))],
                 (tmpFit$Low[[j]]+sdLow)[round(seq(1,length(lowseq),length.out=10))],
                 length=0.05,angle=90, code=3,col=cols[4])
         }
       }

       if(grepl("byID",split[1], ignore.case=TRUE)& !is.null(file)){
          dev.off()
       }
     }
     if(grepl("measure",split[1], ignore.case=TRUE)&!is.null(file)){
         dev.off()
     }
  }
  if(grepl("none",split[1], ignore.case=TRUE)& !is.null(file)){
      dev.off()
  }

}



##### model fitting
.errorBarsSetup<- function(data){
    if(is.null(dim(data)))
    med <- apply(data,2,median)
    sdd <- apply(data,2,sd)
    if(any(is.na(sdd))) sdd <-rep(0,ncol(data))
    return(list("median"=med,"sd"=sdd))
}
