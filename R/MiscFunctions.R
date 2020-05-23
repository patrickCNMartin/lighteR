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

slotApply <- function(x,FUN,...){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),...)
    }
    result
}

slotExtractParam <- function(x,FUN,time){
    cl <- class(x)
    result <- list()
    for(i in slotNames(cl)){
        result[[i]] <- FUN(slot(x,i),time,i)
    }
    result
}

slotAssign <- function(x,y){
    cl <- class(x)

    for(i in slotNames(cl)){

        slot(x,i)<-y[[i]]
    }
    return(x)
}

slotSubset <- function(x,dim =1,template){
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

slotCombine <- function(x,y){
    cl <- class(x)

    for(i in slotNames(cl)){
        tmp <- slot(x,i)
        tmp <- rbind(tmp,y[[i]])
        slot(x,i)<-tmp
    }
    return(x)
}

## setting light time points

setTimePoints <- function(seed,timePoints){
    #time <- new("time")
    seed@meta.param@timePoints <- timePoints

    return(seed)
}

## set Starting gumbel
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
            df<-data.frame(df, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")])),
                                        each=nrow(dat)))
        } else {

            df<- melt(dat,id.vars= c("diskID","Zone"))
            df$variable <- type[1]
            df<-data.frame(df, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","Zone")])),
                                        each=nrow(dat)))
        }




    } else if(fit.to[1] == "allPlants"){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
            if(any(colnames(dat) %in% c("plot","pedigree","line","stem"))){
                tmp<- melt(dat,id.vars= c("diskID","plot","pedigree","line","stem"))
                tmp$variable <- type[1]
                tmp<-data.frame(tmp, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")])),
                                            each=nrow(dat)))
            } else {

                tmp<- melt(dat,id.vars= c("diskID","Zone"))
                tmp$variable <- type[1]
                tmp<-data.frame(df, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","Zone")])),
                                            each=nrow(dat)))
            }

            df[[i]]<-tmp
        }


    } else if(fit.to[1]=="medianPlant"){

        df<- vector("list", nrow(dat))
        for(i in seq_along(df)){
            if(any(colnames(dat) %in% c("plot","pedigree","line","stem"))){
                tmp<- melt(dat,id.vars= c("diskID","plot","pedigree","line","stem"))
                tmp$variable <- type[1]
                tmp<-data.frame(tmp, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","plot","pedigree","line","stem")])),
                                              each=nrow(dat)))
            } else {

                tmp<- melt(dat,id.vars= c("diskID","Zone"))
                tmp$variable <- type[1]
                tmp<-data.frame(df, time=rep(seq_len(ncol(dat[,!colnames(dat) %in% c("diskID","Zone")])),
                                             each=nrow(dat)))
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
    ## Dirty zone
    zone <-gsub(" ","", zone)
    zone <- gsub("missing","", zone)
    ## matching

    mat <- match(zone,tag)

    dip <- dip[mat,"OverCompTime"]
    return(dip)
}
