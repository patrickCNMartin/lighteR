################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Parameter computing
################################################################################
################################################################################


## using median values over replicates to fit model
## might clean things up
.setMedian <- function(df){
    loc<-do.call("rbind",lapply(df,"[[","value"))

    med <-apply(loc,2,median)

    if(nrow(loc)<2){
       sd <- rep(0,ncol(loc))
    } else{
       sd<-apply(loc,2,sd)
    }

    ## replace and add new data
    for(i in seq_along(df)){
       df[[i]]$value <-med
       df[[i]]$sd<-sd
    }
    return(df)
}

.findDip<-function(dat){
     start <-1
     low<-1

     while(dat[low]>dat[low+1]){
        low<-low+1
     }

     count<-low

     while(dat[count]<dat[count+1] & (count+1)<length(dat)){
        count<-count+1
     }


     end <- low+count

     return(end-start)
}

.findDropTime <- function(dat){

     start <-1
     low<-1
     while(dat[low]>dat[low+1]){
        low<-low+1
     }
     return(low-start)
}

.findPlateau<-function(dat){
    start<-.findDip(dat)
    end<-length(dat)
    return(c(start,end,end-start))
}

.rate<-function(dat){

    maxi <- as.numeric(dat[which(dat==max(dat))[1]])
    maxiLoc <- as.numeric(which(dat==max(dat))[1])
    mini <- as.numeric(dat[which(dat==min(dat))[1]])
    miniLoc<-as.numeric(which(dat==min(dat))[1])

    t<-as.numeric(abs(maxiLoc-miniLoc))
    if(maxiLoc>miniLoc){
        rate<-(maxi-mini)/t
    }else{

        rate<-(mini-maxi)/t
    }
    return(rate)
}
