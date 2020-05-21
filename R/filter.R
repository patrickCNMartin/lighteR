################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Filtering
################################################################################
################################################################################


.dataConsitency<-function(x,dataType=c("NPQ","XE","EF","OE"),time=c(40,80)){
	## clean

    x<-x[which(!names(x)%in%c("diskID","Zone"))]
    x<-as.matrix(x)
    rownames(x) <-NULL
    colnames(x) <- NULL
        # extracting basic param
        # NOTE this can be changed to add more light consistency measures


    boundaries <- apply(x,1,function(dat,time){
                        tmp <- vector("list",length(time))
                        names(tmp) <- as.character(time)
                        # adding first point
                        time <- c(1,time)

                        for(i in seq(1,l=length(time)-1)){
                             tmp[[i]] <- c(dat[time[i]],
                                           dat[time[i+1]])
                        }
                        return(tmp)
                    },time = time)

    ## making the checks

   passVec <- rep(TRUE,nrow(x))
   for(rows in seq_len(nrow(x))){
       tmpPass <-c()
       for(bound in seq_along(boundaries[[rows]])){

           if(dataType[1] != "XE"){
               tmpPass <- c(tmpPass,boundaries[[rows]][[bound]][1] < boundaries[[rows]][[bound]][2])
           } else{
               tmpPass <- c(tmpPass,boundaries[[rows]][[bound]][1] > boundaries[[rows]][[bound]][2])
           }

       }
       ## This is where I need to fix this shit first
       ## At the moment it is just using high light for consistency
       ## Need to find a way to consider
       ## This is something that needs to be fixed
       # passVec[bound] <- !sum(tmpPass) ==0

       passVec[rows] <- !sum(tmpPass) ==0
    }

    return(passVec)
}


.modelThresholdSelection<-function(model,residualThreshold=5,time,modelType){

      high<-model$High
	    low<-model$Low

      resiH<-c()
	    resiL<-c()

      if(modelType=="combinedFit"){
            localResiH<-res(high)
            localResiL<-res(low)
            timeHigh<-length(localResiH)/time[1]
            timeLow<-length(localResiL)/(time[2]-time[1])

            for(i in seq_len(timeHigh)){
                 resiH<-c(resiH,sum((res(high)[seq(i,length(localResiH),by=timeHigh)])^2))

                 resiL<-c(resiL,sum((res(low)[seq(i,length(localResiL),by=timeLow)])^2))

            }

      } else{
        for(i in seq_along(high)){

	        resiH<-c(resiH,sum(res(high[[i]])^2))

		      resiL<-c(resiL,sum(res(low[[i]])^2))

        }
      }


	resH <- resiH <residualThreshold
	resL <- resiL <residualThreshold
  #if(sum(resH & resL)==0)browser()
	return(resH | resL)
}


.modelR2Selection<-function(model,data,R2=0.5,time=c(40,80),modelType){
	## Model set up for pseudo R2
  modelClass<-modelType
	modelType<-model$model
  high<-model$High
	low<-model$Low
 	## Data set up for pseudo R2

	datHigh<-data[,seq(6,by=1, length.out=time[1])]
	datLow<-data[,seq(from=time[1]+6,by=1, length.out=time[2]-time[1])]
        ## set up other shit for this
	seqHigh<-seq_len(time[1])
	seqLow<-seq(1,time[2]-time[1])
  resiH<-c()
	resiL<-c()

  if(modelClass=="combinedFit"){
        localResiH<-fitted(high)
        localResiL<-fitted(low)
        timeHigh<-length(localResiH)/time[1]
        timeLow<-length(localResiL)/(time[2]-time[1])
        #browser()
        for(i in seq_len(timeHigh)){
             resiH<-c(resiH,.intMSE(datHigh[i,],localResiH[seq(i,length(localResiH),by=timeHigh)]))


             resiL<-c(resiL,.intMSE(datLow[i,],localResiL[seq(i,length(localResiL),by=timeLow)]))
        }

  } else{
        for(i in seq_along(high)){
              resiH<-c(resiH,.intMSE(datHigh[i,],.extractFittedModel(high[[i]],seqHigh,modelType[1])))
              resiL<- c(resiL,.intMSE(datLow[i,],.extractFittedModel(low[[i]],seqLow,modelType[2])))

        }

  }
	resiHi <- resiH <R2
	resiLo <- resiL <R2
#if(sum(resiHi)==0 | sum(resiLo)==0)browser()
	return(resiHi & resiLo)
}

## returns index of runs to drop or to retain
selectPlants <- function(plants,measure=c("NPQ","XE","EF","OE"),method=c("MSE",0.005)){

    ## first lets get time
    ## need to add a check here
    if(length(plants@time@timePoints)>0){
        time <- plants@time@timePoints
    } else {
        warning("No Time Points have been set - using default 40 - 80")
        time <- c(40,80)
    }


    ## next check what type of
    ## check we will be doing
    models <- sum(unlist(slotApply(plants@models,length))) == 0
    if(models){
        template <- vector("list", length(slotNames(plants@measures)))
        names(template)<- slotNames(plants@measures)
        for(i in seq_along(measure)){
            template[[measure[i]]] <- .dataConsitency(slot(plants@measures,measure[i]),measure[i],time)
        }
        ## Filtering over measures
        template <- !apply(do.call("cbind",template),1,sum) < length(measure)

        retain <- slotSubset(plants@measures,1,template)
        dropped <- slotSubset(plants@measures,1,!template)

        plants@retain <- slotAssign(plants@retain,retain)
        plants@dropped <- slotAssign(plants@dropped,dropped)
        if(sum(sapply(dropped, nrow))>0 &
           sum(unlist(slotApply(plants@origin, function(origin){lapply(origin,length)})))>1){
            message("Re-generating sample Origin after filtering")
            plants <- getPlants(plants,splitby = plants@originType)
        }

    } else {
message("blabla")
    }

    return(plants)
}










