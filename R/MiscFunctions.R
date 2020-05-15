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
