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

        x<-x[which(!names(x)%in%c("ID","plot","pedigree","line","stem"))]
        x<-as.numeric(x)
        # highlight only
        x<-x[seq(1,time[1])]

	      minX<-x[1]
        ## this is maybe a bit harsh
	## this depends on what has been done experimentally
        #minxnext <-x[2]
        maxX<-x[length(x)]
	if(dataType[1]!="XE"){
        	if(minX>maxX ){
             		return(FALSE)
		} else{
	     		return(TRUE)
		}
	} else{
		return(TRUE)
	}
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
selectDisks <- function(data,model=NULL,dataType=c("NPQ","XE","EF","OE"),method=c("MSE",0.005), time=c(40,80)){
    ## this is crap! can be imporved and changed
    ## There is no need for this change
    ## but it's late and this faster to change

    if(method[1]=="MSE"){
       R2<-c(T,method[2])
       residualThreshold <-F

    } else if(method[1]=="RSS"){
        R2<-c(F,method[2])
        residualThreshold <-c(T,method[2])
    } else{
      stop("Unknown Filter method! Either MSE or RSS")
    }


    if(is.null(model)){
        classType<-class(data)
        buffer<-data
        for(i in seq_along(data)){
            for(j in seq_along(data[[i]])){
                locIdx<-apply(data[[i]][[j]],1,.dataConsitency,dataType=names(data)[i],time=time)
                buffer[[i]][[j]]<-locIdx
            }
        }
    } else{
      buffer<-data
      modelType<-class(model)
      for(i in seq_along(data)){
          for(j in seq_along(data[[i]])){
            if(residualThreshold[1]){
                  locResi <- .modelThresholdSelection(model[[i]][[j]],residualThreshold[2],time,modelType=modelType)

            } else if(R2[1]){
                  locResi <-.modelR2Selection(model[[i]][[j]],data[[i]][[j]],R2[2],time,modelType=modelType)
            }
              buffer[[i]][[j]]<-locResi
          }
      }
    }
    return(buffer)
}




#### actually filtering disks
## we are just writing a wrapped for select disks

filterDisks <- function(data,idx,keep=c("retain","drop")){
        classType<-class(data)
        if(keep[1]=="retain"){
            dropLoc<-vector("list", length(idx[[1]]))
            for(i in seq_along(idx[[1]])){
                 primer <- idx[[1]][[i]]
                 for(j in seq(2,length(idx))){

                     dropLoc[[i]]<-c(primer & idx[[j]][[i]])

                 }
            }
#browser()
            ### filtering out disks
            for(i in seq_along(dropLoc)){
               for(j in seq_along(data)){
                    data[[j]][[i]]<-data[[j]][[i]][dropLoc[[i]],]
               }
            }
            ## removing disks
            for(i in seq_along(data)){
                loc<-sapply(data[[i]],function(x){
                            return(ifelse(nrow(x)>0,T,F))
                            })
                data[[i]]<-data[[i]][loc]
                class(data[[i]])<-classType
            }
            class(data)<-classType
        }else if(keep[1]=="drop"){
          dropLoc<-vector("list", length(idx[[1]]))

          for(i in seq_along(idx[[1]])){
               primer <- !idx[[1]][[i]]
               for(j in seq(2,length(idx))){
                   dropLoc[[i]]<-c(primer | !idx[[j]][[i]])

               }

          }
          ### filtering out disks
          for(i in seq_along(dropLoc)){
             for(j in seq_along(data)){
                  data[[j]][[i]]<-data[[j]][[i]][dropLoc[[i]],]
             }
          }
          ## removing disks
          for(i in seq_along(data)){
              loc<-sapply(data[[i]],function(x){
                          return(ifelse(nrow(x)>0,T,F))
                          })
              data[[i]]<-data[[i]][loc]
              class(data[[i]])<-classType
          }
          class(data)<-classType

        }else{
            stop("Ooops! I am not sure what argument you are parsing to keep
                  Either retain or drop are accepted.")
        }
        return(data)

}
