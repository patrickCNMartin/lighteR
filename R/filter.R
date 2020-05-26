################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Filtering
################################################################################
################################################################################


.dataConsitency<-function(x,dataType=c("NPQ","XE","EF","OE"),time=c(40,80)){
	## clean

    x<-x[which(!names(x)%in%c("diskID","Zone"))]
    x<-as.matrix(x)
    rownames(x) <-NULL
    colnames(x) <- NULL
        # extracting basic param
        # NOTE this can be changed to add more light consistency measures


    boundaries <- apply(x,1,function(dat,time){
                        tmp <- vector("list",length(time))
                        names(tmp) <- as.character(time)
                        # adding first point
                        time <- c(1,time)

                        for(i in seq(1,l=length(time)-1)){
                             tmp[[i]] <- c(dat[time[i]],
                                           dat[time[i+1]])
                        }
                        return(tmp)
                    },time = time)

    ## making the checks

   passVec <- rep(TRUE,nrow(x))
   for(rows in seq_len(nrow(x))){
       tmpPass <-c()
       for(bound in seq_along(boundaries[[rows]])){

           if(dataType[1] != "XE"){
               tmpPass <- c(tmpPass,boundaries[[rows]][[bound]][1] < boundaries[[rows]][[bound]][2])
           } else{
               tmpPass <- c(tmpPass,boundaries[[rows]][[bound]][1] > boundaries[[rows]][[bound]][2])
           }

       }
       ## This is where I need to fix this shit first
       ## At the moment it is just using high light for consistency
       ## Need to find a way to consider
       ## This is something that needs to be fixed
       # passVec[bound] <- !sum(tmpPass) ==0

       passVec[rows] <- !sum(tmpPass) ==0
    }

    return(passVec)
}


.rmsq<-function(x,y){
    ## might need to add some checks in here

    x<-as.vector(as.matrix(x))
    y<-as.vector(as.matrix(y))
    return(cor(x,y,method="spearman")^2)
}
.intMSE<-function(gp,chip){
    mse <- sum((1/length(gp))*(gp-chip)^2)
    return(mse)
}


################################################################################
################################################################################
################################################################################
#MSE



.RSSSelection <- function(models,seed,trait,threshold,time,cores=1,origin=FALSE){
    modelType <- all(sapply(models, length) == 3)
    if(modelType){


        if(!origin){
            zone <- seed[,colnames(seed) %in% c("diskID","Zone")]
            tag <- zone
        } else {
            zone <- seed[,colnames(seed) %in% c("diskID","plot","pedigree","line","stem")]
        }
        zone <- apply(zone,1,paste,collapse="")

        dip <- .quickSelect(trait,zone)
        timeLoc <- data.frame("startHighLight" = rep(1,nrow(seed)),
                              "endHighLight" = rep(time[1],nrow(seed)),
                              "startLowLight" = rep(time[1]+1, nrow(seed)),
                              "OverCompTime" = (time[1]+1) + dip,
                              "OverCompTime" = (time[1]+2) + dip,
                              "endLowLight" = rep(time[2],nrow(seed)))

    } else {
        timeLoc <- data.frame("startHighLight" = rep(1,nrow(seed)),
                              "endHighLight" = rep(time[1],nrow(seed)),
                              "startLowLight" = rep(time[1]+1, nrow(seed)),
                              "endLowLight" = rep(time[2],nrow(seed)))
    }
    ## seed split by row


    if(!origin){
        seed <- split(seed, seq(nrow(seed)))
        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))
        models <- mcmapply(.RSSSelect,models,seed,timeLoc,MoreArgs = list(threshold), mc.cores=cores)
    } else {


        seed <- split(seed, seq(nrow(seed)))
        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))
        models <- mapply(.RSSSelect,models,seed,timeLoc,MoreArgs = list(threshold,origin=TRUE))
    }


    return(models)
}

.RSSSelect <- function(model,data,time,threshold=5, origin=FALSE){
    modelLocal <- c()
    if(!origin){
        dataLocal <- as.vector(as.matrix(data[,!colnames(data) %in% c("diskID","Zone")]))

    }else {
        dataLocal <- as.vector(as.matrix(data[,!colnames(data) %in% c("diskID","plot","pedigree","line","stem")]))
    }

    timeLoc <- vector("list", length(model))
    count <- 1
    for(t in seq(1,by=2,length.out =length(model))){
        timeLoc[[count]] <- c(time[t], time[t+1])
        count <- count +1
    }

    if(sum(!is.na(model))<length(model)){
        warning("At least one model failed! Model Goodness of fit will be approximated on remaining models")
    }
    model <- model[!is.na(model)]

    timeLoc <- lapply(timeLoc[!is.na(model)], unlist)
    for(mod in seq_along(model)){


        modelLocal <- c(modelLocal,as.vector(res(model[[mod]]))^2)
    }


    RSSThresh <- sum(modelLocal)

    return(RSSThresh < threshold)
}







.MSESelection <- function(models,seed,trait,threshold,time,cores=1,origin=FALSE){

    ## Time extraction
    if(origin){
        #models <- .orderModels(models)
        modelType <- all(sapply(models, length) == 3)
    } else {
        modelType <- all(sapply(models, length) == 3)
    }

    if(modelType){

        if(!origin){
            zone <- seed[,colnames(seed) %in% c("diskID","Zone")]
            tag <- zone
        } else {
            zone <- seed[,colnames(seed) %in% c("diskID","plot","pedigree","line","stem")]
            tag <- zone
        }

        zone <- apply(zone,1,paste,collapse="")

        dip <- .quickSelect(trait,zone)
        ### Taking care of weird over comp times
        overcomp <-(time[1]+1) + dip
        overcomp2 <-(time[1]+1) + dip +1

        overcomp[overcomp >time[2]] <- time[2]
        overcomp2[overcomp2 >time[2]] <- time[2]

        timeLoc <- data.frame("startHighLight" = rep(1,nrow(seed)),
                              "endHighLight" = rep(time[1],nrow(seed)),
                              "startLowLight" = rep(time[1]+1, nrow(seed)),
                              "OverCompTime" = overcomp,
                              "OverCompTime" = overcomp2,
                              "endLowLight" = rep(time[2],nrow(seed)))

    } else {
        timeLoc <- data.frame("startHighLight" = rep(1,nrow(seed)),
                              "endHighLight" = rep(time[1],nrow(seed)),
                              "startLowLight" = rep(time[1]+1, nrow(seed)),
                              "endLowLight" = rep(time[2],nrow(seed)))
    }
    ## seed split by row
    if(!origin){
        seed <- split(seed, seq(nrow(seed)))
        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))

        models <- mcmapply(.MSESelect,models,seed,timeLoc,MoreArgs = list(threshold), mc.cores=cores)
    } else {


        seed <- split(seed, seq(nrow(seed)))
        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))


        models <- mapply(.MSESelect,models,seed,timeLoc,MoreArgs = list(threshold,origin=TRUE))
    }

    return(models)

}


.MSESelect<-function(model,data,time,threshold=0.05,origin=FALSE){

    ## This is assuming that it comes in one big ass data frame
    ## The question is
    ## Would it be worth it to just do a full merge of both models and then run MSE on
    ## the whole thing instead of having high light and low light? That might be easier
    modelLocal <- c()
    if(!origin){
        dataLocal <- as.vector(as.matrix(data[,!colnames(data) %in% c("diskID","Zone")]))

    }else {
        dataLocal <- as.vector(as.matrix(data[,!colnames(data) %in% c("diskID","plot","pedigree","line","stem")]))
    }


    timeLoc <- vector("list", length(model))
    count <- 1

    for(t in seq(1,by=2,length.out =length(model))){

        timeLoc[[count]] <- c(time[t], time[t+1])
        count <- count +1
    }

    if(sum(!is.na(model))<length(model)){
        warning("At least one model failed! Model Goodness of fit will be approximated on remaining models")
        tmp <- which(is.na(model))

        colRem <- seq(timeLoc[[tmp]][[1]], timeLoc[[tmp]][[2]])
        if(length(colRem)!=1){
            dataLocal <-dataLocal[-colRem]
        }

    }
    model <- model[!is.na(model)]
    timeLoc <- lapply(timeLoc[!is.na(model)], unlist)
    loc <- c()
    for(mod in seq_along(model)){

        ti <- seq(1,(timeLoc[[mod]][[2]]-timeLoc[[mod]][[1]])+1)


        modelLocal <- c(modelLocal,.extractFittedModel(model[[mod]],ti,names(model)[mod]))
    }
     #if(length(dataLocal) != length(modelLocal))browser()

     MSEThresh <- suppressWarnings(.intMSE(dataLocal,modelLocal))



	return(MSEThresh < threshold)
}

#' Extract measures from seed object
#'
#' @param seed a seed object
#' @param measure which measure should be used to select plants/samples ("NPQ","XE","EF","OE")
#' @param method if models have been computed, descirbes which method should be used to assess goodness of fit ("MSE","RSS")
#' @param threshold if models have been computed, threshold at which samples should be removed. Anything above threshold is removed.
#' @param cores number of cores used for selecting plants
#' @return Seed object with filtered data
selectPlants <- function(seed,measure=c("NPQ","XE","EF","OE"),method=c("MSE","RSS"),threshold= c(0.05,0.5), cores= 1){

    ## first lets get time
    ## need to add a check here
    if(length(seed@meta.param@timePoints)>0){
        time <- seed@meta.param@timePoints
    } else {
        message("No Time Points have been set - using default 40 - 80")
        time <- c(40,80)
    }


    ## next check what type of
    ## check we will be doing
    models <- sum(unlist(.slotApply(seed@models,length))) == 0
    if(models){
        template <- vector("list", length(slotNames(seed@measures)))
        names(template)<- slotNames(seed@measures)
        for(i in seq_along(measure)){
            template[[measure[i]]] <- .dataConsitency(slot(seed@measures,measure[i]),measure[i],time)
        }
        ## Filtering over measures
        template <- !apply(do.call("cbind",template),1,sum) < length(measure)

        retain <- .slotSubset(seed@measures,1,template)
        dropped <- .slotSubset(seed@measures,1,!template)

        seed@retain <- .slotAssign(seed@retain,retain)
        seed@dropped <- .slotAssign(seed@dropped,dropped)
        seed@meta.data$dropped <- c("retain"=sum(template),"dropped"=sum(!template))
        if(sum(sapply(dropped, nrow))>0 &
           sum(unlist(.slotApply(seed@origin, function(origin){lapply(origin,length)})))>1){
            message("Re-generating sample Origin after filtering")
            seed <- getOrigin(seed,splitby = seed@meta.param@originType)
        }

    } else {
        models <- seed@models
        modelType <- seed@meta.param@models
        fit.to <- seed@meta.param@fitto

        ### There is a better way of doing this but at the moment i will keep this
        meth <- method[1]
        if(meth == "MSE" & length(threshold)==2){
            message("No MSE threshold has been specificied - using default 0.05")
            threshold <- threshold[1]
        } else if(meth == "RSS" &length(threshold)==2){
            message("No RSS threshold has been specified - using default 0.5")
            threshold <- threshold[2]
        }

        ## extracting models
        template <- vector("list", length(measure))
        names(template) <-measure
        for(i in seq_along(measure)){
            tmp <- modelType[[measure[i]]]

            if(any(grepl("No", tmp, ignore.case=TRUE))){
                warning(paste(measure[i],"models have not been computed - skipping measure"))
                next()
            }
            is.retain.empty <- sum(unlist(.slotApply(seed@retain, length)))==0
            is.origin.empty <- sum(unlist(.slotApply(seed@origin, length)))==0

            if(is.retain.empty){
                plant <- slot(seed@measures,measure[i])
            } else if(!is.retain.empty) {
                plant <- slot(seed@retain, measure[i])
            }
            if(!is.origin.empty){

                plant <- slot(seed@origin, measure[i])
            }
            is.trait.empty <- sum(unlist(.slotApply(seed@traits, length)))==0
            if(is.trait.empty){
                trait <- NULL
            } else{
                trait <- slot(seed@traits, measure[i])
            }
            ## Running threshold selection MSE
            if(meth == "MSE"){

                mods <- slot(models,measure[i])
                is.origin.empty <- sum(unlist(.slotApply(seed@origin, length)))==0
                if(is.origin.empty){
                    template[[i]] <- .MSESelection(mods,plant,trait,threshold,time,cores)


                } else {

                    template[[i]] <- mcmapply(.MSESelection,mods,plant,
                                              MoreArgs = list(trait,threshold,time,fit.to,origin=TRUE),mc.cores =cores)

                }
            } else if(meth == "RSS"){
                mods <- slot(models,measure[i])
                is.origin.empty <- sum(unlist(.slotApply(seed@origin, length)))==0
                if(is.origin.empty){
                    template[[i]] <- .RSSSelection(mods,plant,trait,threshold,time,cores)
                } else {

                    template[[i]] <- mcmapply(.RSSSelection,mods,plant,
                                              MoreArgs = list(trait,threshold,time,fit.to,origin =TRUE),mc.cores=cores)
                }
            } else {
                stop("Unkown model selection method - Select from MSE or RSS")
            }
        }

        ### Subsetting
        ## Checking if origin is empty again
        is.retain.empty <- sum(unlist(.slotApply(seed@retain, length)))==0
        is.origin.empty <- sum(unlist(.slotApply(seed@origin, length)))==0
        if(!is.origin.empty){
            tempSub <- template
           templ <- lapply(template,function(x){
                if(!is.null(x)){
                    return(do.call("c",x))
                } else {
                    return(NULL)
                }})
            tempTag <- lapply(template,function(x){
                    if(!is.null(x)){
                        t <- unlist(lapply(x,names))
                        return(t)
                    } else {
                        return(NULL)
                    }})
            tempTag <- apply(do.call("cbind",tempTag),1,unique)
            template <- !apply(do.call("cbind",templ),1,sum) < sum(!grepl("No",modelType))
        } else {

            template <- !apply(do.call("cbind",template),1,sum) < sum(!grepl("No",modelType))
        }

        ## This is not great but if will have to do for now
        template[which(is.na(template))] <- FALSE
        ## at leat for dropped and retained
        ### Subsetting
        if(is.retain.empty & is.origin.empty){
            retain <- .slotSubset(seed@measures,1,template)
            dropped <- .slotSubset(seed@measures,1,!template)
            seed@retain <- .slotAssign(seed@retain,retain)
            seed@dropped <- .slotCombine(seed@dropped,dropped)
        } else if(!is.retain.empty & is.origin.empty){
            retain <- .slotSubset(seed@retain,1,template)
            dropped <- .slotSubset(seed@retain,1,!template)
            seed@retain <- .slotAssign(seed@retain,retain)
            seed@dropped <- .slotCombine(seed@dropped,dropped)
        } else if(is.retain.empty & !is.origin.empty){
            retain <- .slotSubsetWithTag(seed@measures,template,tempTag)
            dropped <- .slotSubsetWithTag(seed@measures,!template,tempTag)
            seed@retain <- .slotAssign(seed@retain,retain)
            seed@dropped <- .slotCombine(seed@dropped,dropped)
        } else {
            retain <- .slotSubsetWithTag(seed@retain,template,tempTag)
            dropped <- .slotSubsetWithTag(seed@retain,!template,tempTag)
            seed@retain <- .slotAssign(seed@retain,retain)
            seed@dropped <- .slotCombine(seed@dropped,dropped)
        }

        if(!is.trait.empty & is.origin.empty){

            traits <- .slotSubset(seed@traits,1,template)

            seed@traits <- .slotAssign(seed@traits,traits)
        } else if(!is.trait.empty & !is.origin.empty){

            traits <- .slotSubsetWithTag(seed@traits,template,tempTag)

            seed@traits <- .slotAssign(seed@traits,traits)
        }

        if(is.origin.empty){
            models <- .slotListSub(seed@models,template)

            seed@retained.models <- .slotAssign(seed@retained.models,models)
        } else {
            seed@retained.models <- .modelSub(seed@models,modelType,measure,tempSub)
            seed@origin <-.generateOrigin(seed@origin,modelType,measure,tempSub)

        }




        fullDrop <- unique(.slotApply(seed@dropped,nrow))
        seed@meta.data$dropped <- c("retain"=sum(template),"dropped"=fullDrop)


    }

    return(seed)
}









