################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Extraction
################################################################################
################################################################################




.extractAndMergeParam<-function(data, model=NULL,dataType=c("NPQ","XE","EF","OE"),time=c(40,80)){
    if(is.null(model)){
         classType<-class(data)
         for(i in seq_along(data)){
             for(j in seq_along(data[[i]])){

                   if(names(data)[i]!="XE"){

                     if(length(grep("filterStep", colnames(data[[i]][[j]])))>0){

                       tmpH<-data[[i]][[j]][,seq(7,time[1]+6)]
                       tmpL<-data[[i]][[j]][,seq(time[1]+7,time[2]+6)]
                       buffer <- matrix(0,ncol=12,nrow(data[[i]][[j]]))

                     }else{
                       tmpH<-data[[i]][[j]][,seq(6,time[1]+5)]
                       tmpL<-data[[i]][[j]][,seq(time[1]+6,time[2]+5)]
                       buffer <- matrix(0,ncol=12,nrow(data[[i]][[j]]))
                     }

                   } else{

                      if(length(grep("filterStep", colnames(data[[i]][[j]])))>0){

                       tmpH<-data[[i]][[j]][,seq(8,time[1]+7)]
                       tmpL<-data[[i]][[j]][,seq(time[1]+7,time[2]+7)]
                       #dark <-data[[i]][[j]][,7]
                       buffer <- matrix(0,ncol=12,nrow(data[[i]][[j]]))
                     }else{

                       tmpH<-data[[i]][[j]][,seq(7,time[1]+6)]
                       tmpL<-data[[i]][[j]][,seq(time[1]+7,time[2]+6)]
                       #dark <-data[[i]][[j]][,7]
                       buffer <- matrix(0,ncol=12,nrow(data[[i]][[j]]))
                     }

                   }
                   colnames(buffer)<-c("startHighLight","minHighLight","maxHighLight",
                                       "InductionTime","LinearRate","minLowLight","startLowLight",
                                       "OverCompTime","RelaxationTime","startPlateau",
                                       "endPlateau","stableLowLightTime")
                   ## could probably improves on this part
                   buffer[,1]<-tmpH[,1]
                   buffer[,2]<-apply(tmpH,1,min)
                   buffer[,3]<-apply(tmpH,1,max)
                   if(names(data)[i]!="XE"){

                      buffer[,4]<-apply(tmpH,1,function(x){return(which(x==max(x))[1])})
                   } else{

                     buffer[,4]<-apply(tmpH,1,function(x){return(which(x==min(x))[1])})
                   }
                   buffer[,5]<-apply(tmpH,1,.rate)
                   buffer[,6]<-apply(tmpL,1,min)
                   buffer[,7]<-tmpL[,1]
                   buffer[,8]<-apply(tmpL,1,.findDip)
                   buffer[,9]<-apply(tmpL,1,.findDropTime)
                   loc<-t(apply(tmpL,1,.findPlateau))

                   buffer[,10]<-loc[,1]
                   buffer[,11]<-loc[,2]
                   buffer[,12]<-loc[,3]
                    if(length(grep("filterStep", colnames(data[[i]][[j]])))>0){
                     data[[i]][[j]]<-cbind(data[[i]][[j]][,1:6],buffer,data[[i]][[j]][,7:ncol(data[[i]][[j]])])

                   }else{
                     data[[i]][[j]]<-cbind(data[[i]][[j]][,1:5],buffer,data[[i]][[j]][,6:ncol(data[[i]][[j]])])

                   }
             }
         }
         class(data)<-classType
         return(data)
    } else{
           classType<-class(data)
           modelType<-class(model)
           data<-.extractAndMergeParam(data,model=NULL,dataType,time)
           tmpMod<-model
           for(i in seq_along(model)){
               for(j in seq_along(model[[i]])){
                   buffer <- matrix(0,ncol=6+time[2],nrow(data[[i]][[j]]))

                   for(k in seq_len(nrow(data[[i]][[j]]))){
                       local <-c(model[[i]][[j]][[3]][1],
                                 coef(model[[i]][[j]][[1]][[1]]),
                                 model[[i]][[j]][[3]][2],
                                 coef(model[[i]][[j]][[2]][[1]]))


                      loc<-.modelFitPlot(model[[i]][[j]],time,modelType)
                      loc <-c(loc[[1]][[1]],loc[[2]][[1]])
                      #if(length(loc)>80)browser()
                       tag<-c("modelHighLight",names(coef(model[[i]][[j]][[1]][[1]])),
                              "modelLowLight",names(coef(model[[i]][[j]][[2]][[1]])),
                              paste0("Fitted_at_time",seq_along(loc)))

                      if(k==1){
                        buffer <-matrix(c(local,loc), ncol=length(c(local,loc)))
                        colnames(buffer)<-tag
                      } else{
                        buffer <-rbind(buffer,c(local,loc))
                      }


                   }


                   #browser()
                   tmpMod[[i]][[j]] <- cbind(data[[i]][[j]][,1:5],buffer)
               }
           }
           class(data)<-classType
           class(model)<-modelType
           return(list(data=data,model=tmpMod))
    }

}

