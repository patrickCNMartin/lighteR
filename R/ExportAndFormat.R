################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Export and formatting
################################################################################
################################################################################



.match.column.Type <- function(df){
    # Removing meta data
    localdf<-df[,12:ncol(df)]

    # Filthy coersion
    subdf<-as.vector(apply(apply(localdf,2,as.numeric),2, function(x){
                                        if(all(is.na(x)))return(FALSE)
                                        if(!all(is.na(x)))return(TRUE)
                                      }))
    idx<-c(rep(FALSE,11),subdf)
    return(idx)
}



.nullConversion <- function(input){
   if(input=="NULL"){
     return(NULL)
   }else{
     return(input)
   }
}

### this was to remove unwated space but it seems like grep can't find the patern
.spaceRemoval<-function(id){

    d<-lapply(id,function(x){
           loc<-grep("",x,ignore.case=T)
           return(x[-loc])

    })
    return(id)
}

.IDtag<-function(x,tag){
     if(length(x)==length(tag)){
       names(x)<-tag
     } else if(length(x)<length(tag)){
       message("ID does not follow template \n")
       message(c(paste(x,collapse=" "),"\n"))
       x<-c(x,rep("missing",length(tag)-length(x)))
       names(x)<-tag
     }
     return(x)
}


.formatConversion <- function(dat, type=c("NPQ","XE","EF","OE"),individual=TRUE){

     ## the is.true function is a didgy and hacky
     ## it only return true if what ever you pass is true otherwise just
     ## FALSE . even if it is numeric, character or what ever
     if(is.false(individual)){

        df<- melt(dat,id.vars=colnames(dat)[1:5])
     # clean this shit up
        df$variable <- type[1]
        df<-data.frame(df, time=rep(seq_len(ncol(dat[,6:ncol(dat)])),each=nrow(dat)))

     } else if(is.true(individual)){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
          tmp<- melt(dat[i,],id.vars=colnames(dat)[1:5])
        # clean this shit up
          tmp$variable <- type[1]
          tmp<-data.frame(tmp, time=rep(seq_len(ncol(dat[i,6:ncol(dat)])),each=nrow(dat[i,])))
          df[[i]]<-tmp
        }


      } else if(individual=="median"){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
          tmp<- melt(dat[i,],id.vars=colnames(dat)[1:5])
        # clean this shit up
          tmp$variable <- type[1]
          tmp<-data.frame(tmp, time=rep(seq_len(ncol(dat[i,6:ncol(dat)])),each=nrow(dat[i,])))
          df[[i]]<-tmp
        }
         if(length(df)==0)browser()
         df<-.setMedian(df)

      }
     return(df)
}






###### filtering functions

#### clean this shit up over the week end make it pretty and robust


## merging filtered disks

mergeFilteredDisks <- function(df1, df2){
    tag1<-names(df1[[1]])
    tag2<-names(df2[[1]])

    merged <-vector("list", length(df1))
    names(merged)<-names(df1)
    merged <-lapply(merged,function(x,tag1,tag2){
                    common<-match(tag1,tag2)
                    diff<-c(match(tag1,tag2),match(tag2,tag1))
                    x<-vector("list", sum(!is.na(common))+sum(is.na(diff)))
                    names(x)<-c(tag1[!is.na(common)],tag1[is.na(common)],tag2[is.na(match(tag2,tag1))])
                    return(x)
                  },tag1=tag1,tag2=tag2)

    ## merging
    for(i in seq_along(merged)){
        locNames<-names(merged[[i]])
        for(j in seq_along(merged[[i]])){
            tmp1<-grep(locNames[j],tag1)
            tmp2<-grep(locNames[j],tag2)

            if(length(tmp1)>0 & length(tmp2)>0){
                #print(paste("i=",i,"j=",j,"binding"))
                merged[[i]][[j]]<-rbind(df1[[i]][[tmp1]],df2[[i]][[tmp2]])
                merged[[i]][[j]]<-cbind(merged[[i]][[j]][,1:5],
                                        "filterStep"=rep("merged",nrow(merged[[i]][[j]])),
                                        merged[[i]][[j]][,6:ncol(merged[[i]][[j]])])
            }else if(length(tmp1)>0 & length(tmp2)==0){
              #print(paste("i=",i,"j=",j,"premod"))
                merged[[i]][[j]]<-df1[[i]][[tmp1]]
                merged[[i]][[j]]<-cbind(merged[[i]][[j]][,1:5],
                                        "filterStep"=rep("predMod",nrow(merged[[i]][[j]])),
                                        merged[[i]][[j]][,6:ncol(merged[[i]][[j]])])

            }else if(length(tmp1)==0 & length(tmp2)>0){
              #print(paste("i=",i,"j=",j,"postmod"))
               merged[[i]][[j]]<-df2[[i]][[tmp2]]
               merged[[i]][[j]]<-cbind(merged[[i]][[j]][,1:5],
                                       "filterStep"=rep("postMod",nrow(merged[[i]][[j]])),
                                       merged[[i]][[j]][,6:ncol(merged[[i]][[j]])])

            }
        }
    }
    return(merged)
}


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
