################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Export and formatting
################################################################################
################################################################################






###### filtering functions




## merging filtered disks


### finalising export function

exportData<- function(data,file=NULL,sep=","){
    for(i in seq_along(data)){
        for(j in seq_along(data[[i]])){

            if(length(grep("fit",class(data),ignore.case=T))>0){
                filename<-paste0(file,names(data)[i],names(data[[i]])[j],"fitted.csv")
            } else{
                filename<-paste0(file,names(data)[i],names(data[[i]])[j],".csv")
            }
            write.csv(data[[i]][[j]], file=filename)
        }
    }

}





generateSummaryFiles <- function(data,file=NULL,sep=","){
         classType<-class(data)

         data<-lapply(data,function(x){
                     return(do.call("rbind",x))
                     })

         for(i in seq_along(data)){
           if(length(grep("fit",classType,ignore.case=T))>0){
               filename<-paste0(file,"_",names(data)[i],"_Summary_fitted.csv")

           } else{
               filename<-paste0(file,"_",names(data)[i],"_Summary.csv")
           }
           write.csv(data[[i]], file=filename)

         }
}



.convertTime<-function(data,imageData){
      if(length(grep("Image",names(imageData)))==0){
        stop("No Image data loaded - batchLoading type should be set to all or image")
      }else{
         time<-lapply(imageData[["Image"]],function(x){
                      x<-x[[grep("OE",names(x))]]
                      time<-as.numeric(x$Time)
                      time<-time-min(time)
                      light<-as.numeric(x$PPFD)

                      return(list("time"=time,"light"=light))
         })
         names(time)<-gsub("ImageData","ZoneData",names(time))
         if(all(names(data)%in%c("data","model"))){
              locDat<-data[["data"]]
              locMod<-data[["model"]]
              for(i in seq_along(locDat)){
                  for(j in seq_along(locDat[[i]])){
                      findTime<-match(as.character(locDat[[i]][[j]][,1]),names(time))
                      localTime<-time[findTime]

                      for(k in seq_along(localTime)){
                          highLight <- localTime[[k]][["time"]][which(
                            localTime[[k]][["light"]]==unique(localTime[[k]][["light"]])[1])]
                          lowLight <- localTime[[k]][["time"]][which(
                              localTime[[k]][["light"]]==unique(localTime[[k]][["light"]])[2])]
                          locDat[[i]][[j]][k,"InductionTime"]<-highLight[locDat[[i]][[j]][k,"InductionTime"]]
                          locDat[[i]][[j]][k,"OverCompTime"]<-lowLight[locDat[[i]][[j]][k,"OverCompTime"]]-lowLight[1]
                          ## this is only when you have a very quick relation time based on other metrics
                          ## this might not be ideal
                          if(length(lowLight[locDat[[i]][[j]][k,"RelaxationTime"]]-lowLight[1])==0){
                            locDat[[i]][[j]][k,"RelaxationTime"]<-0
                          }else{
                            locDat[[i]][[j]][k,"RelaxationTime"]<-lowLight[locDat[[i]][[j]][k,"RelaxationTime"]]-lowLight[1]
                          }

                          locDat[[i]][[j]][k,"startPlateau"]<-lowLight[locDat[[i]][[j]][k,"startPlateau"]]
                          locDat[[i]][[j]][k,"endPlateau"]<-lowLight[locDat[[i]][[j]][k,"endPlateau"]]
                          locDat[[i]][[j]][k,"stableLowLightTime"]<-lowLight[locDat[[i]][[j]][k,"stableLowLightTime"]]-lowLight[1]

                      }

                  }
              }

              return(list("data"=locDat,"model"=locMod))

         } else{
           locDat<-data[["data"]]
           for(i in seq_along(locDat)){
               for(j in seq_along(locDat[[i]])){
                   findTime<-match(as.character(locDat[[i]][[j]][,1]),names(time))
                   localTime<-time[findTime]
                   for(k in seq_along(localTime)){
                       highLight <- localTime[[k]][["time"]][which(
                         localTime[[k]][["light"]]==unique(localTime[[k]][["light"]])[1])]
                       lowLight <- localTime[[k]][["time"]][which(
                           localTime[[k]][["light"]]==unique(localTime[[k]][["light"]])[2])]
                       locDat[[i]][[j]][k,"InductionTime"]<-highLight[locDat[[i]][[j]][k,"InductionTime"]]
                       locDat[[i]][[j]][k,"OverCompTime"]<-lowLight[locDat[[i]][[j]][k,"OverCompTime"]]-lowLight[1]
                       ## this is only when you have a very quick relation time based on other metrics
                       ## this might not be ideal
                       if(length(lowLight[locDat[[i]][[j]][k,"RelaxationTime"]]-lowLight[1])==0){
                         locDat[[i]][[j]][k,"RelaxationTime"]<-0
                       }else{
                         locDat[[i]][[j]][k,"RelaxationTime"]<-lowLight[locDat[[i]][[j]][k,"RelaxationTime"]]-lowLight[1]
                       }
                       locDat[[i]][[j]][k,"startPlateau"]<-lowLight[locDat[[i]][[j]][k,"startPlateau"]]
                       locDat[[i]][[j]][k,"endPlateau"]<-lowLight[locDat[[i]][[j]][k,"endPlateau"]]
                       locDat[[i]][[j]][k,"stableLowLightTime"]<-lowLight[locDat[[i]][[j]][k,"stableLowLightTime"]]-lowLight[1]


                   }


               }
            }
        }
        return(locDat)
    }
}
