################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Measure Extraction
################################################################################
################################################################################

lightResponse <- function(seed, measures = c("NPQ","XE","EF","OE"),norm = c("local","global","none")){
    ## making a few checks
    if(any(!norm %in% c("local","global","none"))){
        stop("Unknown norm type - Availbale normalisation : local, global, none")
    } else {
        norm <-norm[1]
    }

    ## First lets check if the measure type are correct
    if(any(!measures %in% c("NPQ","XE","EF","OE"))){
        stop("Unknown measure type - Availbale measures : NPQ, XE, EF, OE")
    }
    ## Extracting Measures
    roots <- seed@roots@Zone
    if(length(roots)==0){
        stop("No Zone data loaded within roots")
    } else if(length(roots) ==1){
        light <- .extractMeasure(roots,names(roots), type = measures)
    } else {
        light <- .batchExtractMeasure(roots, type = measures)
    }
    ## extract measure type from each file
    lightRes <- vector("list", length(measures))
    names(lightRes) <- measures
    for(m in seq_along(lightRes)){
        tmp <- lapply(light,"[[",measures[m])
        tmp <- do.call("rbind", tmp)
        if(norm == "local"){

            tmpLoc <- t(apply(tmp[,!colnames(tmp) %in% c("diskID","Zone")],1,function(x){return(x/max(x))}))
            tmp <- cbind(tmp[,colnames(tmp) %in% c("diskID","Zone")],tmpLoc)
            lightRes[[m]] <- tmp
        } else if(norm =="global"){
            localMax <- max(apply(tmp[,!colnames(tmp) %in% c("diskID","Zone")],1,function(x){return(max(x))}))
            tmpLoc <- t(apply(tmp[,!colnames(tmp) %in% c("diskID","Zone")],1,function(x,l){return(x/l)},localMax))
            light <- cbind(tmp[,colnames(tmp) %in% c("diskID","Zone")],tmpLoc)
            lightRes[[m]] <- tmp
        } else {
            lightRes[[m]] <- tmp
        }

    }



    ### Assinging measures
    types <- c("NPQ","XE","EF","OE")
    localType <- lightRes[names(lightRes) %in% types]
    measureLocal <- match(types,names(localType))
    localEnvir <- environment()

    mapply(function(measureLocal,types,localEnvir,localType){
        if(is.na(measureLocal)){
            assign(types,list(),envir=localEnvir)
        }else{
            assign(types,localType[[measureLocal]],envir=localEnvir)
        }},measureLocal,types,MoreArgs=list(localEnvir,localType))


    lightResp <- new("measures",
                     NPQ = NPQ,
                     XE =XE,
                     EF =EF,
                     OE = OE)

    seed@measures <- lightResp
    return(seed)
}



.extractMeasure <- function(data,ID,type=c("NPQ","XE","EF","OE"),threshold=5,single=TRUE){
    if(single){
        data <- data[[1]]
    }
    datasub <- .nonZeroIndex(data,threshold)


    ## lets custom split this
    datasubSplit<- vector("list", length(type))
    names(datasubSplit)<-type

    # ID remap
    Zone <- as.character(datasub[,"Zone"])
    diskID <-rep(ID, length(Zone))
    for(i in seq_along(datasubSplit)){
        datasubSplit[[i]]<-datasub[,grep(type[i],colnames(datasub))]
        datasubSplit[[i]]<-cbind(diskID,Zone, datasubSplit[[i]])
    }

    measure <- list(datasubSplit)
    names(measure) <- ID
    return(measure)

}

## extracting evey measure in batches for zone data

.batchExtractMeasure <- function(data,type=c("NPQ","XE","EF","OE"),threshold=5){


    ID<-names(data)
    local<-mapply(function(zoneData,ID,type=type){

            return(.extractMeasure(zoneData,ID, type=type,
                                   threshold=threshold,single = FALSE))

    },data,ID,MoreArgs=list(type=type))

    return(local)
}
