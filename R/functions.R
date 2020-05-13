################################################################################
############################ NPQ analysis ######################################
################################################################################


## Sourcing other functions
## Split to make it cleaner




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



### formatShift for model fitting
