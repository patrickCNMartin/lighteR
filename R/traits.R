################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Selecting Parameters
################################################################################
################################################################################

#' Extract traits from a seed object - simple trait and model traits
#'
#' @param seed a seed object
#' @param measure light measure that traits should be extracted for
#' @param cores numer of cores used for analysis
#' @return return a seed object with extract traits
getTraits <- function(seed,measure = c("NPQ","XE","EF","OE"),cores=1){
    # Extracting time points
    if(length(seed@meta.param@timePoints)>0){
        time <- seed@meta.param@timePoints
    } else {
        message("No Time Points have been set - using default 40 - 80")
        time <- c(40,80)
    }

    ## check for models
    models <- sum(unlist(.slotApply(seed@models,length))) == 0
    if(models){
        ## Just checking what I will use to extract param
        ## No need to extract param from measure if some disk have already been filtered
        is.retain.empty<- sum(unlist(.slotApply(seed@retain,length))) == 0
        if(is.retain.empty){
            measures <- seed@measures
        } else {
            measures <- seed@retain
        }
        ##
        param <- .slotExtractParam(measures,.extractParam,time)
        seed@traits <- .slotAssign(seed@traits,param)
    } else {
        is.retained.models.empty <-sum(unlist(.slotApply(seed@retained.models,length))) == 0

        if(is.retained.models.empty){

            models <- seed@models
            #models <- slotUnlist(models)
            #models <- slotAssign(seed@models,models)

        } else {
            models <- seed@retained.models
        }

        modelType <- seed@meta.param@models

        is.origin.empty <-sum(unlist(.slotApply(seed@origin,length))) == 0
        is.retain.empty<- sum(unlist(.slotApply(seed@retain,length))) == 0
        is.trait.empty<- sum(unlist(.slotApply(seed@traits,length))) == 0
        template <- vector("list", length(measure))
        names(template)<-measure
        for(i in seq_along(measure)){
            tmp <- modelType[[measure[i]]]

            if(any(grepl("No", tmp, ignore.case=TRUE))){
                 warning(paste(measure[i],"models have not been computed - skipping measure"))
                next()
            }

            if(is.retain.empty){
                plant <- slot(seed@measures,measure[i])
            } else if(!is.retain.empty) {
                plant <- slot(seed@retain, measure[i])
            }
            if(!is.origin.empty){

                plant <- slot(seed@origin, measure[i])
            }

            if(is.trait.empty){
                trait <- NULL
            } else{
                trait <- slot(seed@traits, measure[i])
            }
            mods <- slot(models,measure[i])
            if(is.origin.empty){

                template[[i]] <- .extractModels(mods,plant,trait,time,origin=FALSE,cores)



            } else {

                template[[i]] <- mcmapply(.extractModels,mods,plant,
                                          MoreArgs = list(trait,time,origin=TRUE),mc.cores =cores)

            }
        }

        if(!is.origin.empty){
            ## re-orient df
            #browser()
            template <- lapply(template, function(tmp){
                                if(!is.null(tmp)){
                                    nmax <-range(unique(unlist(sapply(tmp, function(x)sapply(x,length)))))

                                    tmp <- suppressWarnings(lapply(tmp,function(x){

                                        do.call("rbind",matrix(x,ncol=length(x)))
                                        }))

                                    for(i in seq_along(tmp)){
                                        if(ncol(tmp[[i]])== min(nmax)){
                                            tmp[[i]]<- cbind(tmp[[i]],rep(NA,max(nmax)-min(nmax)))
                                        }
                                    }
                                    tmp <- do.call("rbind", tmp)
                                    return(as.data.frame(tmp))
                                }else{
                                    return(data.frame())
                                }
                            })

        } else {
            for(i in seq_along(template)){
                if(!is.null(template[[i]])){

                    template[[i]] <- as.data.frame(do.call("rbind",template[[i]]))
                } else {
                    template[[i]] <- data.frame()
                }
            }
        }

        if(is.trait.empty){
            seed@traits <- .slotAssign(seed@traits,template)
        } else {
            seed@traits <- .slotAddTraits(seed@traits,template)
        }



    }
    return(seed)
}

.extractParam <- function(df,time,measure){

    ##
    paramType <-c("startHighLight","endHighLight","minHighLight","maxHighLight",
                  "InductionTime","LinearRate","startLowLight","endLowLight","minLowLight",
                  "OverCompTime","RelaxationTime","startPlateau",
                  "endPlateau","stableLowLightTime")
    param <- as.data.frame(matrix(0,ncol = length(paramType), nrow = nrow(df)))
    colnames(param) <- paramType


    ## High Light param
    if(any(colnames(df) %in% c("plot","pedigree","line","stem"))){
        dftmp<- df[,!colnames(df) %in% c("diskID","plot","pedigree","line","stem")]
        tags <-df[,colnames(df) %in% c("diskID","plot","pedigree","line","stem")]
        dfh <- as.matrix(dftmp[,seq(1,time[1])])
        dfl <- as.matrix(dftmp[seq(time[1]+1,time[2])])
        ## Cleaning up for debugging purpose
    } else {
        dftmp<- df[,!colnames(df) %in% c("diskID","Zone")]
        tags<- df[,colnames(df) %in% c("diskID","Zone")]
        dfh <- as.matrix(dftmp[,seq(1,time[1])])
        dfl <- as.matrix(dftmp[seq(time[1]+1,time[2])])
    }

    param$startHighLight <- apply(dftmp,1,"[[",1)
    param$endHighLight <- apply(dftmp,1,"[[",time[1])
    param$minHighLight <- apply(dfh,1,min)
    param$maxHighLight <- apply(dfh,1,max)
    if(measure!="XE"){
        param$InductionTime<-apply(dfh,1,function(x){return(which(x==max(x))[1])})
    } else{
        param$InductionTime<-apply(dfh,1,function(x){return(which(x==min(x))[1])})
    }
    param$LinearRate<-apply(dfh,1,.rate)
    param$minLowLight<-apply(dfl,1,min)
    param$startLowLight<-apply(dftmp,1,"[[",time[1]+1)
    param$endLowLight<-apply(dftmp,1,"[[",time[2])
    param$OverCompTime<-apply(dfl,1,.findDip)
    param$RelaxationTime<-apply(dfl,1,.findDropTime)
    loc<-t(apply(dfl,1,.findPlateau))

    param$startPlateau<-loc[,1]
    param$endPlateau<-loc[,2]
    param$stableLowLightTime<-loc[,3]

    ## Adding back names
    param <- cbind(tags,param)
    return(param)

}


.extractModels <- function(models,seed,trait,time,origin =FALSE,cores=1){

    ## Time extraction
    if(origin == TRUE){
        #models <- .orderModels(models)
        modelType <- all(sapply(models, length) == 3)
        zone <- seed[,colnames(seed) %in% c("diskID","plot","pedigree","line","stem")]
        tag <- zone
    } else {

        modelType <- all(sapply(models, length) == 3)
        zone <- seed[,colnames(seed) %in% c("diskID","Zone")]
        tag <- zone
    }


    if(modelType){


        zone <- apply(zone,1,paste,collapse="")

        dip <- .quickSelect(trait,zone)

        ### Taking care of weird over comp times
        overcomp <-(time[1]+1) + median(dip)
        overcomp2 <-(time[1]+1) + median(dip) +1

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
    if(origin==FALSE){

        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))

        tag <- split(tag,seq(nrow(tag)))

        tag <- lapply(tag,function(x){
            as.character(as.vector(as.matrix(x)))})

        models <- mcmapply(.SelectFitted,models,timeLoc,tag,SIMPLIFY = FALSE ,mc.cores=cores)
        models <- lapply(models, unlist)


    } else {

        timeLoc<- split(timeLoc, seq(nrow(timeLoc)))

        tag <- split(tag, seq(nrow(tag)))
        tag <- lapply(tag,function(x){
            as.character(as.vector(as.matrix(x)))})

        models <- mapply(.SelectFitted,models,timeLoc,tag,MoreArgs = list(origin=TRUE),SIMPLIFY= FALSE)

    }

    return(models)

}


.SelectFitted <- function(model,time,tags,origin=FALSE){

    ### Using filtering function as template
    ## they work in similar ways
    modelLocal <- c()
    coefs <- c()
    timeLoc <- vector("list", length(model))
    count <- 1

    for(t in seq(1,by=2,length.out =length(model))){
        timeLoc[[count]] <- c(time[t], time[t+1])
        count <- count +1
    }


    nas <- is.na(model)
    tag <- c()
    for(mod in seq_along(model)){

        ti <- seq(1,(timeLoc[[mod]][[2]]-timeLoc[[mod]][[1]])+1)
        if(any(names(timeLoc[[mod]]) %in% "OverCompTime.1")){
            if(timeLoc[[mod]]$OverCompTime.1 ==timeLoc[[mod]]$endLowLight) next()
        }

        if(nas[mod]){
            modelLocal <- c(modelLocal,rep(NA,length(ti)))
            resi <- NA
            coefs <-c(coefs,resi,rep(NA,3))
        } else {
            modelLocal <- c(modelLocal,.extractFittedModel(model[[mod]],ti,names(model)[mod]))
            resi <- sqrt(sum(model[[mod]]$residuals^2)/(length(model[[mod]]$residuals)-2))
            coefs <- c(coefs,coef(model[[mod]]),resi)
        }
        tag <- c(tag, paste0(rep(names(model)[mod],length(ti)),mod))
    }
    if(length(tag)< length(modelLocal)){
         tag <- c(tag, rep(NA,length(modelLocal)-length(tag)))
    } else if(length(tag)> length(modelLocal)){
        modelLocal <- c(modelLocal, rep(NA,length(tag)-length(modelLocal)))
    }

    names(modelLocal) <- tag

    names(coefs)[names(coefs) %in% ""] <- "RSE"
    mods <- c(as.character(tags),coefs,"fitted_data",modelLocal)

    return(mods)
}

