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

    if(ZoneError){

    }
}
ZoneError <- sapply(roots@Zone,is.character)
ImageError <- sapply(roots@Image,is.character)
browser()
if(sum(ZoneError)!=0| sum(ImageError)!=0){
    warning("Plate Error while Loading - check seed meta data for failed plates")
    plateError <- list("ZoneError" = roots@Zone[ZoneError],
                       "ImageError" = roots@Image[ImageError])
    seed@meta.data <- plateError
} else {
    seed@meta.data <- list("Plate Error" = "All Plates Loaded Succesfully")
}

