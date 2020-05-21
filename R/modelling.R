################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Modelling
################################################################################
################################################################################

modelPlants <- function(seed, models = list("NPQ" = c("beta",3)), fit.to=c("plant","allPlants","medianPlant"), cores=1){
    ## Lets get time
    if(length(seed@meta.param@timePoints)>0){
        time <- seed@meta.param@timePoints
    } else {
        message("No Time Points have been set - using default 40 - 80")
        time <- c(40,80)
    }
    ## setting data
    if(sum(unlist(slotApply(seed@origin, length)))==0){
        message("Empty Origin slot - samples will be modelled individually")
        if(sum(unlist(slotApply(seed@retain,length)))==0){
            plants <- seed@measures

        }else {
            plants <- seed@retain
        }
        origin <- FALSE

    } else {
        plants <- seed@origin
        origin <- TRUE
    }

    ## Modelling only what is required
    for(mod in seq_along(models)){

        if(!origin){
            ### This section is based on the assumption that there is not origin slot
            ### This just means that we assume that the "fit.to" argument is set to indiv
            ### No format conversion
            ## setting up modelling
            measure <- slot(plants,names(models)[mod])


            modelsInternal <- .dispatchModel(models[[mod]],seed@meta.param@nlsrStart)
            ### If - model length check
            if(names(models)[mod] != "NPQ" & length(modelsInternal$model)==3){
                stop("Triple model set-up only supported for NPQ")
            } else if(names(models)[mod] == "NPQ" & length(modelsInternal$model)==3){

                #####################################################################
                #####################################################################
                ## Check if Traits have been computed
                if(sum(unlist(slotApply(seed@traits, nrow)))==0){
                    stop("Traits have not yet been computed - use getTraits function first")
                } else {
                #browser()
                    dip <- seed@traits@NPQ$OverCompTime
                    time <- data.frame("startHighLight" = rep(1,nrow(measure)),
                                       "endHighLight" = rep(time[1],nrow(measure)),
                                       "startLowLight" = rep(time[1]+1, nrow(measure)),
                                       "OverCompTime" = time[1] + dip,
                                       "endLowLight" = rep(time[2],nrow(measure))
                                       )

                    timeSplit <- .splitDfs(time,cores)
                #browser()
                    measure <- .splitDfs(measure,cores)
                    models[[mod]] <- mcmapply(.applyTripleModelFit,chunk=measure,timeSplit=timeSplit,
                                              MoreArgs = list(modelsInternal = modelsInternal),mc.cores = cores)
                }
                #####################################################################
                #####################################################################
            }else {
                #####################################################################
                #####################################################################
                ## Two models or less
                measure <- .splitDfs(measure,cores)

                models[[mod]] <- mclapply(measure,.applyModelFit,time =time,
                                        modelsInternal = modelsInternal,mc.cores = cores)
            }
            #####################################################################
            #####################################################################

        } else {
            #####################################################################
            #####################################################################
            ## Section for Origin
            measure <- slot(plants,names(models)[mod])


            modelsInternal <- .dispatchModel(models[[mod]],seed@meta.param@nlsrStart)
            ### If - model length check
            if(names(models)[mod] != "NPQ" & length(modelsInternal$model)==3){
                stop("Triple model set-up only supported for NPQ")
            } else if(names(models)[mod] == "NPQ" & length(modelsInternal$model)==3){

                #####################################################################
                #####################################################################
                ## Check if Traits have been computed
                if(sum(unlist(slotApply(seed@traits, nrow)))==0){
                    stop("Traits have not yet been computed - use getTraits function first")
                } else {

                    models[[mod]] <- mclapply(measure,.applyTripleModelFitOrigin,
                                              seed = seed , time =time,
                                              modelsInternal =modelsInternal,
                                              fit.to =fit.to, mc.cores=cores)

                }
                #####################################################################
                #####################################################################
            }else {
                #####################################################################
                #####################################################################
                ## Two models or less
                #measure <- .splitDfs(measure,cores)
                models[[mod]] <- mclapply(measure,.applyModelFitOrigin,
                                          seed = seed , time =time,
                                          modelsInternal =modelsInternal,
                                          fit.to =fit.to, mc.cores=cores)

            }
            #####################################################################
            #####################################################################
        }


    }
    #####################################################################
    #####################################################################
    ## Done
    return(models)



}



## just cleaning up number of model used
.dispatchModel <- function(model,start){

    if(length(model)==1){
        if(is.numeric(model)){
            tmp<-rep("poly",2)
            start <-rep(model,2)
        } else{
            tmp<-rep(model,2)
            start<-list(start,start)
        }

    } else if(length(model) == 2){

        tmp<-model
        start<-list(start,start)
        ## check which one are poly if any
        classCheck<-grep("([0-9]+).*$",tmp)

        ## change the ones that need to be chnaged
        if(length(classCheck)==1){
            start[[classCheck]]<-as.numeric(tmp[classCheck])
            tmp[classCheck]<-"poly"
        } else if(length(classCheck)==2){
            start[classCheck]<-as.numeric(tmp[classCheck])
            tmp[classCheck]<-"poly"
        }
    }else if(length(model)==3){

        ## set up multi model template
        tmp<- model
        start<-list(start, start, start)
        ## check which one are poly if any
        classCheck<-grep("([0-9]+).*$",tmp)
        ## change the ones that need to be chnaged
        if(length(classCheck)!=0){
            start[[classCheck]]<-as.numeric(tmp[classCheck])
            tmp[classCheck]<-"poly"
        }
    }

    return(list("model"=tmp,"start"=start))
}

## computing metrics of interest

.gumbel <- function(data, time) {

    tmp<-coef(data)

    return((tmp[1]/tmp[2]) * exp( -(time-tmp[3])/tmp[2] - exp(-(time-tmp[3])/tmp[2])))
}

.laplace <- function(data,time) {
    tmp<-coef(data)
    return(tmp[1]/(1+tmp[2]*exp(-tmp[3]*(-time))))
}

.poly<-function(data, time){

    return(.polyRefit(fitted(data),time))
}

.exp<-function(data, time){

    #res <-.polyRefit(fitted(data),time)
    res <- exp(predict(data,list(time=time)))
    return(res)
}

.recip<-function(data, time){

    #res <-.polyRefit(fitted(data),time)
    res <- 1/.polyRefit(fitted(data),time)
    return(res)
}

.Logistic<-function(data, time){
    res<-fitted(data)
    return(res)
}

.negBi<-function(data, time){
    res<-predict(data,list(time=time))
    return(res)
}

.gamma<-function(data, time){
    res<-fitted(data,list(time=time))
    return(res)
}
.beta<-function(data, time){

    res<-.polyRefit(fitted(data),time)
    return(res)
}


.fitGumbel <- function(dat, start){

    return(nlsr::nlxb(value ~ (a/b)*exp(-(time-c)/b -exp(-(time-c)/b)), start=start, data=dat))
}


.fitLaplace <- function(dat, start){
    return(nlsr::nlxb(value ~ a/(1+b*exp(-c*(-time))), start=start, data=dat))

}

.fitExp <- function(dat,start){

    return(lm(log(value)~time,data=dat))
}

.fitPoly <- function(dat, poly){
    return(lm(value~poly(time,poly),data=dat))
}

.fitRecip<-function(dat){
    return(lm(1/value ~time, data= dat))
}

.fitLogistic<-function(dat){
    return(glm(value ~ time, data=dat, family=binomial(link="logit")))
}

.fitNegBi<-function(dat){

    return(glm.nb(value ~ time, data=dat))
}

.fitGamma<-function(dat){
    return(glm(value ~ time, data=dat, family=Gamma()))
}

.fitBeta<-function(dat){
    ## values must be strictly between 0 and 1 non included

    dat$value[dat$value==1] <-0.999999999
    return(betareg(value ~ time, data=dat))
}


.fitModel <- function(dat, fitParam, type){

    res<- switch(type,
                 gumbel=.fitGumbel(dat, fitParam),
                 laplace=.fitLaplace(dat, fitParam),
                 poly=.fitPoly(dat,fitParam),
                 exp=.fitExp(dat),
                 recip=.fitRecip(dat),
                 logistic=.fitLogistic(dat),
                 negBi=.fitNegBi(dat),
                 gamma=.fitGamma(dat),
                 beta=.fitBeta(dat))
    return(res)

}

.extractFittedModel <-function(data,time,type){

    res<-switch(type,
                gumbel=.gumbel(data,time),
                laplace=.laplace(data,time),
                exp=.exp(data,time),
                poly=.poly(data,time),
                recip=.recip(data,time),
                logistic=.Logistic(data,time),
                negBi=.negBi(data,time),
                gamma=.gamma(data,time),
                beta=.beta(data,time))
    return(res)
}



.applyTripleModelFit<-function(chunk,timeSplit,modelsInternal, fit.to ="plant"){

    ## extract some shit
    start <- modelsInternal$start
    model <- modelsInternal$model

    ## adding tags
    if(any(colnames(chunk) %in% c("plot","pedigree","line","stem"))){
        tag <- chunk[,colnames(chunk) %in% c("diskID","plot","pedigree","line","stem")]
    } else {
        tag <- chunk[, colnames(chunk) %in%  c("diskID","Zone")]
    }
    ## further chunking
    fitted<-vector("list",nrow(chunk))
    for(row in seq_len(nrow(chunk))){
        h <- timeSplit[row,"startHighLight"]:timeSplit[row,"endHighLight"]
        dfh <- cbind(tag[row,],chunk[row,h])

        d <- timeSplit[row,"startLowLight"]:timeSplit[row,"OverCompTime"]
        dip <- cbind(tag[row,],chunk[row,d])
        l <- timeSplit[row,"OverCompTime"]:timeSplit[row,"endLowLight"]
        dfl <-cbind(tag[row,],chunk[row,l])

        NewChunk <- list("HighLight"=dfh,"OverCompTime" = dip, "lowLight" = dfl)

        fittedModel <-vector("list", length(model))
        names(fittedModel) <- model

        for(mod in seq_along(model)){
            formattedChunk <- .formatConversion(NewChunk[[mod]],fit.to = fit.to)
            fittedModel[[mod]] <- tryCatch(.fitModel(formattedChunk,fitParam = start[[mod]],type = model[[mod]]),
                                           error = function(cond){
                                               return(NA)
                                           })
        }
        fitted[[row]] <- fittedModel
    }
    return(fitted)
}


.applyModelFit <- function(chunk,time,modelsInternal,fit.to ="plant"){

    ## extract some shit
    start <- modelsInternal$start
    model <- modelsInternal$model
    if(any(colnames(chunk) %in% c("plot","pedigree","line","stem"))){
        tag <- chunk[,colnames(chunk) %in% c("diskID","plot","pedigree","line","stem")]
    } else {
        tag <- chunk[, colnames(chunk) %in%  c("diskID","Zone")]
    }
    ## further chunking
    fitted<-vector("list",nrow(chunk))
    for(row in seq_len(nrow(chunk))){
        h <- seq(1,time[1])
        dfh <- cbind(tag[row,],chunk[row,h])

        l <- seq(time[1]+1,time[2])
        dfl <-cbind(tag[row,],chunk[row,l])

        NewChunk <- list("HighLight"=dfh,"lowLight" = dfl)

        fittedModel <-vector("list", length(model))
        names(fittedModel) <- model
        for(mod in seq_along(model)){
            formattedChunk <- .formatConversion(NewChunk[[mod]],fit.to = fit.to)
            fittedModel[[mod]] <- tryCatch(.fitModel(formattedChunk,fitParam = start[[mod]],type = model[[mod]]),
                                           error = function(cond){
                                               return(NA)
                                           })
        }
        fitted[[row]] <- fittedModel
    }
    return(fitted)
}

.applyTripleModelFitOrigin <- function(measure,seed,time,modelsInternal,fit.to,cores){


    dip <- seed@traits@NPQ
    zone <- measure[,colnames(measure) %in% c("plot","pedigree","line","stem")]
    zone <- apply(zone,1,paste,collapse=" ")
    dip <- .quickSelect(dip,zone)
    timeLoc <- data.frame("startHighLight" = rep(1,nrow(measure)),
                      "endHighLight" = rep(time[1],nrow(measure)),
                      "startLowLight" = rep(time[1]+1, nrow(measure)),
                      "OverCompTime" = time[1] + dip,
                      "endLowLight" = rep(time[2],nrow(measure))
    )


    models <- .applyTripleModelFit(chunk=measure,timeSplit=timeLoc,
                                    modelsInternal = modelsInternal,fit.to=fit.to)

    return(models)
}


.applyModelFitOrigin <- function(measure,seed,time,modelsInternal,fit.to,cores){


    models <- .applyModelFit(measure,time =time,
                              modelsInternal = modelsInternal,fit.to=fit.to)
    return(models)
}
