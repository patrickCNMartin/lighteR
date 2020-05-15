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

### data extraction
## forzone data
## you havent really done anything for all data

.extractByID <- function(data, splitby=c("plot","pedigree","line","stem"),tagID=c("plot","pedigree","line","stem")){

  if(is(data)[1]=="data.frame"){
    IDs<- list("Zone"=as.character(data[,"Zone"]))
  } else if(class(data)=="byMeasure"){
    IDs<- lapply(data, function(x){return(as.character(x[,"Zone"]))})
  }
  tag <- names(data)
  ## format data check
  ## this is me being lazy
  # i dont't want to have two large chunks of code for the two Options
  ## so i will coerce to list for computing splits
  ## then extract from list element if needed
  if(is(data)[1]=="data.frame"){
    data <-list(data)
    classCheck<-TRUE
  } else {
    classCheck<-FALSE
  }



  # split by ID type

   IDSplit<- lapply(IDs,strsplit," ")

   #IDSplit<-lapply(IDSplit,spaceRemoval)

   IDSplit<-lapply(IDSplit,function(x, tags){return(lapply(x, IDtag,tags))},
                   tags=tagID)


   ## setting up for splitting by adding factor levels for each option
   IDSplit <- lapply(IDSplit,function(dat){return(do.call("rbind",dat))})

   for(i in seq_along(IDSplit)){
       IDSplit[[i]]<-as.data.frame(IDSplit[[i]])
       data[[i]]<-data[[i]][,-c(2)]
       data[[i]]<-cbind("ID"=data[[i]][,1],IDSplit[[i]], data[[i]][,2:ncol(data[[i]])])

   }

   ## now lets split this bad boy
   for(i in seq_along(data)){

       data[[i]]<-split(data[[i]],lapply(splitby,function(split, data){return(data[,split])},data[[i]]), drop=TRUE)
   }

   ## finalising data output
  if(classCheck){
     ## removing that from that list
     # and assigne class for later code
     data<-data[[1]]
     class(data)<-"byIDCombined"
  }else{
     names(data)<-tag
     class(data)<-"byID"
  }


  return(data)


}

### batch extraction for ID

.batchExtractByID <-function(data,splitby=c("plot","pedigree","line","stem"),
                            tagID=c("plot","pedigree","line","stem"),combine=TRUE,norm=TRUE){

        ## step one check for plate errors in data and remove

        ## check for plate errors
        ## select one that are not plate errors
        NoError <- sapply(data, function(x){
                   if(is(x)[1]=="byMeasure"){
                       return(TRUE)
                   }else{
                       return(FALSE)
                   }
                   })

        localData <- data[NoError]

        ## check to see if data sets need to be combined
        #if(combine){


           ## sub setted data / hopefully data will be retrained
           if(length(localData)>1){
                ## creating combined list for each measure
                tmp <- vector("list", length(localData[[1]]))
                names(tmp)<- names(localData[[1]])
                tmp<-lapply(tmp, function(tmp,data){
                            tmp<-vector("list", length(data))
                            #names(tmp)<-names(data)
                     },localData)

                ### first extract all measure independantly

                for(i in seq_along(localData)){

                    for(j in seq_along(localData[[i]])){
                        if(norm){
                            Zone<-as.character(localData[[i]][[j]][,"Zone"])
                            ID<-as.character(localData[[i]][[j]][,"ID"])
			                      loc<-localData[[i]][[j]][,3:ncol(localData[[i]][[j]])]
                            loc<-t(apply(loc,1, function(x){return(x/max(x))}))
                            loc<-data.frame(ID,Zone,loc)


                        }else{
                            loc<-localData[[i]][[j]]
                        }
                        tmp[[j]][[i]]<- loc
                    }

                }
                ## build combined data frames

                if(combine){
                  tmp<-lapply(tmp,function(df){return(do.call("rbind",df))})
                }

                ## extract by ID

                tmp<-lapply(tmp,extractByID, splitby,tagID)

                if(!combine){
                    tmp<-lapply(tmp,function(x,tag){
                              names(x)<-tag
                              return(x)
                         },tag=names(localData))
                }

           } else{
             stop("Cannot combine less than 2 data sets")
           }

    if(combine){
      class(tmp)<-"byIDCombined"
    }else{
      class(tmp)<-"byID"
    }


     return(tmp)
}
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

### data extraction
## forzone data
## you havent really done anything for all data

.extractMeasure <- function(data,ID,type=c("NPQ","XE","EF","OE"),threshold=5){

    datasub <- nonZeroIndex(data,threshold)


     ## lets custom split this
     datasubSplit<- vector("list", length(type))
     names(datasubSplit)<-type

     # ID remap
     Zone <- as.character(datasub[,"Zone"])
     ID <-rep(ID, length(Zone))
        for(i in seq_along(datasubSplit)){
             datasubSplit[[i]]<-datasub[,grep(type[i],colnames(datasub))]
             datasubSplit[[i]]<-cbind(ID,Zone, datasubSplit[[i]])
        }
     class(datasubSplit)<-"byMeasure"
     return(datasubSplit)

}

## extracting evey measure in batches for zone data

batchExtractMeasure <- function(data,type=c("NPQ","XE","EF","OE"),threshold=5){

    zoneOnly<-data$zone
    ID<-names(zoneOnly)
    local<-mapply(function(zoneData,ID,type=type){
                          if(is(zoneData)[1]=="character"){
                             return("Plate Error - Check for salvaging")
                          } else{
                             return(extractMeasure(zoneData,ID, type=type,threshold=threshold))
                          }
    },zoneOnly,ID,MoreArgs=list(type=type))
    class(local)<-"byMeasure"
    return(local)
}

extractByID <- function(data, splitby=c("plot","pedigree","line","stem"),tagID=c("plot","pedigree","line","stem")){

  if(is(data)[1]=="data.frame"){
    IDs<- list("Zone"=as.character(data[,"Zone"]))
  } else if(class(data)=="byMeasure"){
    IDs<- lapply(data, function(x){return(as.character(x[,"Zone"]))})
  }
  tag <- names(data)
  ## format data check
  ## this is me being lazy
  # i dont't want to have two large chunks of code for the two Options
  ## so i will coerce to list for computing splits
  ## then extract from list element if needed
  if(is(data)[1]=="data.frame"){
    data <-list(data)
    classCheck<-TRUE
  } else {
    classCheck<-FALSE
  }



  # split by ID type

   IDSplit<- lapply(IDs,strsplit," ")

   #IDSplit<-lapply(IDSplit,spaceRemoval)

   IDSplit<-lapply(IDSplit,function(x, tags){return(lapply(x, IDtag,tags))},
                   tags=tagID)


   ## setting up for splitting by adding factor levels for each option
   IDSplit <- lapply(IDSplit,function(dat){return(do.call("rbind",dat))})

   for(i in seq_along(IDSplit)){
       IDSplit[[i]]<-as.data.frame(IDSplit[[i]])
       data[[i]]<-data[[i]][,-c(2)]
       data[[i]]<-cbind("ID"=data[[i]][,1],IDSplit[[i]], data[[i]][,2:ncol(data[[i]])])

   }

   ## now lets split this bad boy
   for(i in seq_along(data)){

       data[[i]]<-split(data[[i]],lapply(splitby,function(split, data){return(data[,split])},data[[i]]), drop=TRUE)
   }

   ## finalising data output
  if(classCheck){
     ## removing that from that list
     # and assigne class for later code
     data<-data[[1]]
     class(data)<-"byIDCombined"
  }else{
     names(data)<-tag
     class(data)<-"byID"
  }


  return(data)


}

### batch extraction for ID

batchExtractByID <-function(data,splitby=c("plot","pedigree","line","stem"),
                            tagID=c("plot","pedigree","line","stem"),combine=TRUE,norm=TRUE){

        ## step one check for plate errors in data and remove

        ## check for plate errors
        ## select one that are not plate errors
        NoError <- sapply(data, function(x){
                   if(is(x)[1]=="byMeasure"){
                       return(TRUE)
                   }else{
                       return(FALSE)
                   }
                   })

        localData <- data[NoError]

        ## check to see if data sets need to be combined
        #if(combine){


           ## sub setted data / hopefully data will be retrained
           if(length(localData)>1){
                ## creating combined list for each measure
                tmp <- vector("list", length(localData[[1]]))
                names(tmp)<- names(localData[[1]])
                tmp<-lapply(tmp, function(tmp,data){
                            tmp<-vector("list", length(data))
                            #names(tmp)<-names(data)
                     },localData)

                ### first extract all measure independantly

                for(i in seq_along(localData)){

                    for(j in seq_along(localData[[i]])){
                        if(norm){
                            Zone<-as.character(localData[[i]][[j]][,"Zone"])
                            ID<-as.character(localData[[i]][[j]][,"ID"])
			                      loc<-localData[[i]][[j]][,3:ncol(localData[[i]][[j]])]
                            loc<-t(apply(loc,1, function(x){return(x/max(x))}))
                            loc<-data.frame(ID,Zone,loc)


                        }else{
                            loc<-localData[[i]][[j]]
                        }
                        tmp[[j]][[i]]<- loc
                    }

                }
                ## build combined data frames

                if(combine){
                  tmp<-lapply(tmp,function(df){return(do.call("rbind",df))})
                }

                ## extract by ID

                tmp<-lapply(tmp,extractByID, splitby,tagID)

                if(!combine){
                    tmp<-lapply(tmp,function(x,tag){
                              names(x)<-tag
                              return(x)
                         },tag=names(localData))
                }

           } else{
             stop("Cannot combine less than 2 data sets")
           }

    if(combine){
      class(tmp)<-"byIDCombined"
    }else{
      class(tmp)<-"byID"
    }


     return(tmp)
}
