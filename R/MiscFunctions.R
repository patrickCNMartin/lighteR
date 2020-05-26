################################################################################
############################ NPQ analysis ######################################
################################################################################



# Misc functions


.singleElememntPos <- function(data){
    rows<- nrow(data)
    if(rows<=2){
        return(TRUE)
    } else{
        return(FALSE)
    }
}

.nonZeroIndex <- function(data,threshold=5){
    ## first let's remove artifacts
    area <- data[,2]
    area <- which(area >= threshold)
    return(data[area,])

}

## the is.true function is a didgy and hacky
## it only return true if what ever you pass is true otherwise just FALSE
## even if it is numeric, character or what ever
.is.true<-function(x){
    if(x==TRUE){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

.is.false<-function(x){
    if(x==FALSE){
        return(TRUE)
    } else {
        return(FALSE)
    }
}



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
        miss <- c(paste(x,collapse=" "))
        warning(paste("MAP ID ",miss, "does not follow template \n"),call. = FALSE)

        x<-c(x,rep("missing",length(tag)-length(x)))
        names(x)<-tag
    }
    return(x)
}


.polyRefit<-function(model,seq){

    div<- length(model)/length(seq)
    model <- model[seq(1,length(model), by=div)]
    return(model)
}

### Plate error checks
.plateError <- function(roots){
    ## First need to check the nature of the roots
    ZoneError <- length(roots@Zone) != 0
    ImageError <- length(roots@Image) != 0

    if(ZoneError & ImageError){
        Zone <- sapply(roots@Zone,is.character)
        Image <- sapply(roots@Image,is.character)
        if(sum(Zone)!=0| sum(Image)!=0){
            warning("Plate Error while Loading - check seed meta data for failed plates")
            plateError <- list("ZoneError" = names(roots@Zone[Zone]),
                               "ImageError" = names(roots@Image[Image]))
            seed.meta.data <- plateError

            roots@Zone <- roots@Zone[!Zone & !Image]
            roots@Image <- roots@Image[!Zone & !Image]
        } else {
            seed.meta.data <- list("Plate Error" = "All Plates Loaded Succesfully")
        }
    } else if(ZoneError & !ImageError){
        Zone <- sapply(roots@Zone,is.character)
        if(sum(Zone)!=0){
            warning("Plate Error while Loading - check seed meta data for failed plates")
            plateError <- list("ZoneError" = names(roots@Zone[Zone]),
                               "ImageError" = "No Image data Loaded")
            seed.meta.data <- plateError
            roots@Zone <- roots@Zone[!Zone]
        } else {
            seed.meta.data <- list("Plate Error" = "All Plates Loaded Succesfully")
        }
    } else if(!ZoneError & ImageError){
        Image <- sapply(roots@Image,is.character)
        if(sum(ImageError)!=0){
            warning("Plate Error while Loading - check seed meta data for failed plates")
            plateError <- list("ZoneError" = "No Zone data Loaded",
                               "ImageError" = names(roots@Image[ImageError]))
            seed.meta.data <- plateError
            roots@Image <- roots@Image[!Image]
        } else {
            seed.meta.data <- list("Plate Error" = "All Plates Loaded Succesfully")
        }
    } else {
        stop("Ooops somthing went wrong - Empty roots - No input data")
    }

    return(list("seed.meta.data" = seed.meta.data,"roots"=roots))
}

.slotApply <- function(x,FUN,...){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),...)
    }
    result
}

.slotUnlist <- function(x){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        if(length(slot(x,i))==1){
            result[[i]] <- slot(x,i)[[1]]
        } else {
            next()
        }
    }
    result
}

.slotExtractParam <- function(x,FUN,time){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),time,i)
    }
    result
}

.slotAssign <- function(x,y){
    cl <- class(x)

    for(i in slotNames(cl)){
        if(is.null(y[[i]]))next()
        slot(x,i)<-y[[i]]
    }
    return(x)
}

.slotSubset <- function(x,dim =1,template){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        if(dim ==1){
            result[[i]] <- slot(x,i)[template,]
        } else if(dim ==2){
            result[[i]] <- slot(x,i)[,template]
        } else {
            stop("Error in number of dimensions")
        }

    }
    result
}

.slotSubsetWithTag <- function(x,template,tag){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
            tmp <-slot(x,i)

            tagInt <- tmp[,colnames(tmp) %in% c("diskID","Zone")]
            tagInt <- apply(tagInt,1,paste, collapse ="")
            tagInt <- gsub(" ","",tagInt)
            names(tagInt) <- NULL
            names(tag)<- NULL
            reorder <- match(tag,tagInt)
            tmp <- tmp[reorder,]
            result[[i]] <-tmp[template,]



    }
    result
}
.slotListSub<- function(x,template){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
            if(length(slot(x,i))!=0){
                tmp <-slot(x,i)
                result[[i]] <- tmp[template]
            } else {
                result[[i]] <- NULL
            }

    }
    result
}


.slotCombine <- function(x,y){
    cl <- class(x)

    for(i in slotNames(cl)){
        tmp <- slot(x,i)
        tmp <- rbind(tmp,y[[i]])
        slot(x,i)<-tmp
    }
    return(x)
}


.slotAddTraits <- function(x,y){
    cl <- class(x)

    for(i in slotNames(cl)){
        tmp <- slot(x,i)

        ## Skipping empty measures
        ## this makes sense as you wont always model every measure
        if(is.null(y[[i]])) next()

        if(sum(c("V1","V2","V3","V4","V5") %in% colnames(y[[i]]))==5){
            tmp <- cbind(tmp,y[[i]][,seq(6,ncol(y[[i]]))])

        } else if(sum(c("V1","V2","V3","V4","V5") %in% colnames(y[[i]]))==2){
            tmp <- cbind(tmp,y[[i]][,seq(3,ncol(y[[i]]))])
        }

        slot(x,i)<-tmp
    }
    return(x)
}
## setting light time points

#' Set Time points used in fluorImage analysis
#'
#' @param seed A seed object
#' @param timePoints a numeric vector containing end points of different light regimes
#' @return returns a seed object with updated time points
setTimePoints <- function(seed,timePoints){
    #time <- new("time")
    seed@meta.param@timePoints <- timePoints

    return(seed)
}

## set Starting gumbel

#' Set Time points used in fluorImage analysis
#'
#' @param seed A seed object
#' @param start a list containing starting values for gumbel optimisation. Please name your list elements a,b, and c
#' @return returns a seed object with updated start points
setStartPoints <- function(seed,start){
    #time <- new("time")
    seed@meta.param@nslrStart <- start

    return(seed)
}

.splitDfs <- function(ranges,cores){

    if(cores>1 & cores < nrow(ranges)){
        SplitRanges <- floor(seq(1,nrow(ranges),length.out=cores+1))
        rangeSet<-vector("list",(cores))

        start <- SplitRanges[1:(length(SplitRanges)-1)]
        end <- c(SplitRanges[2:(length(SplitRanges)-1)]-1,SplitRanges[length(SplitRanges)])

        for(i in seq_len(cores)){
            rangeSet[[i]]<-ranges[start[i]:end[i],]
        }

    } else if(cores>1 & cores==nrow(ranges)) {
        idx<-seq_len(nrow(ranges))
        rangeSet<-vector("list",cores)
        for(i in seq_len(cores)){
            rangeSet[[i]]<-ranges[idx[i],]
        }
    } else {
        rangeSet <- list(ranges)
    }
    return(rangeSet)
}



.formatConversion <- function(dat, type=c("NPQ","XE","EF","OE"),fit.to="plant"){


    if(fit.to[1] == "plant"){

        if(any(colnames(dat) %in% c("plot","pedigree","line","stem"))){
            df<- melt(dat,id.vars= c("diskID","plot","pedigree","line","stem"))
            df$variable <- type[1]
            time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")]))),
                        each=nrow(dat))
            df<-data.frame(df, time=time)
        } else {

            df<- melt(dat,id.vars= c("diskID","Zone"))
            df$variable <- type[1]
            time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","Zone")]))),each=nrow(dat))
            df<-data.frame(df, time=time)
        }




    } else if(fit.to[1] == "allPlants"){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
            if(any(colnames(dat) %in% c("plot","pedigree","line","stem"))){
                tmp<- melt(dat,id.vars= c("diskID","plot","pedigree","line","stem"))
                tmp$variable <- type[1]
                time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")]))),
                            each=nrow(dat))
                tmp<-data.frame(tmp, time=time)
            } else {

                tmp<- melt(dat,id.vars= c("diskID","Zone"))
                tmp$variable <- type[1]
                time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","Zone")]))),each=nrow(dat))
                tmp<-data.frame(tmp, time=time)

            }

            df[[i]]<-tmp
        }


    } else if(fit.to[1]=="medianPlant"){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
            if(any(colnames(dat) %in% c("plot","pedigree","line","stem"))){
                tmp<- melt(dat,id.vars= c("diskID","plot","pedigree","line","stem"))
                tmp$variable <- type[1]
                time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")]))),
                            each=nrow(dat))

                tmp<-data.frame(tmp, time=time)
            } else {

                tmp<- melt(dat,id.vars= c("diskID","Zone"))
                tmp$variable <- type[1]
                time <- rep(seq_len(ncol(data.frame(dat[,!colnames(dat) %in% c("diskID","Zone")]))),each=nrow(dat))
                tmp<-data.frame(tmp, time=time)
            }

            df[[i]]<-tmp
        }


        df<-.setMedian(df)

    } else {
        stop("fit.to should be one of the following: plant, allPlants, medianPlant")
    }
    return(df)
}


.quickSelect<- function(dip,zone){
    tag <- apply(dip[,c("diskID","Zone")],1,paste, collapse="")

    ## Dirty stuff
    tag <- gsub(" ","", tag)
    tag <- gsub("missing","", tag)
    names(tag) <- NULL
    ## Dirty zone
    zone <-gsub(" ","", zone)
    zone <- gsub("missing","", zone)
    names(zone) <- NULL
    ## matching

    mat <- match(zone,tag)
    if(any(is.na(mat)))browser()

    dip <- dip[mat[!is.na(mat)],"OverCompTime"]
    return(dip)
}


.orderModels <- function(models){

    mType <- unique(as.vector(sapply(models, names)))

    template <- vector("list",length(models))
    names(template)<- names(models)

    for(i in seq_along(template)){
        template[[i]] <- vector("list",unique(sapply(models[[i]],length)))
        names(template[[i]]) <- lapply(models[[i]],names)[[1]]

        for(j in seq_along(template[[i]])){
            template[[i]][[j]] <- vector("list", length(models[[i]]))
            names(template[[i]][[j]]) <- names(models[[i]])

            for(k in seq_along(models[[i]])){
                template[[i]][[j]][[k]] <- models[[i]][[k]][[j]]
            }
        }

    }
    return(template )

}


.generateOrigin <- function(origin,mType,measure,template){
    for(i in seq_along(measure)){
        tmp <- slot(origin, measure[i])


        if(any(mType[[measure[i]]]=="No")){
            next()
        }

        for(ori in seq_along(tmp)){

            tmpSub <- !apply(do.call("cbind",lapply(template,"[[",ori)),1,sum) <sum(!grepl("No",mType))
            if(all(tmpSub==FALSE)){
                tmp[[ori]] <- NA
            } else {
                tmp[[ori]]<- tmp[[ori]][tmpSub,]
            }

        }
        tmp <- tmp[!is.na(tmp)]
        slot(origin,measure[i])<-tmp

    }


    return(origin)
}


.modelSub <- function(model,mType,measure,template){
    retMod <- new("retained.models")
    for(i in seq_along(measure)){
        tmp <- slot(model, measure[i])


        if(any(mType[[measure[i]]]=="No")){
            next()
        }

        for(mod in seq_along(tmp)){

            tmpSub <- !apply(do.call("cbind",lapply(template,"[[",mod)),1,sum) <sum(!grepl("No",mType))
            if(all(tmpSub==FALSE)){
                tmp[[mod]] <- NA
            } else {
                tmp[[mod]]<- tmp[[mod]][tmpSub]
            }

        }

        tmp <- tmp[!is.na(tmp)]

        slot(retMod,measure[i])<-tmp

    }


    return(retMod)
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
