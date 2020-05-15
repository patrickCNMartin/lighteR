################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Measure Extraction
################################################################################
################################################################################

lightResponse <- function(seed, measures = c("NPQ","XE","EF","OE")){
    roots <- seed@roots@Zone
    if(length(roots)==0){
        stop("No Zone data loaded within roots")
    } else if(length(roots) == 1){
        measures <- .ExtractMeasure(roots, type = measures)
    } else {
        measures <- .batchExtractMeasure(roots, type = measures)
    }
}



.extractMeasure <- function(data,ID,type=c("NPQ","XE","EF","OE"),threshold=5){

    datasub <- .nonZeroIndex(data,threshold)


    ## lets custom split this
    datasubSplit<- vector("list", length(type))
    names(datasubSplit)<-type

    # ID remap
    Zone <- as.character(datasub[,"Zone"])
    ID <-rep(ID, length(Zone))
    for(i in seq_along(datasubSplit)){
        datasubSplit[[i]]<-datasub[,grep(type[i],colnames(datasub))]
        datasubSplit[[i]]<-cbind(ID,Zone, datasubSplit[[i]])
    }

    return(datasubSplit)

}

## extracting evey measure in batches for zone data

.batchExtractMeasure <- function(data,type=c("NPQ","XE","EF","OE"),threshold=5){


    ID<-names(data)
    local<-mapply(function(zoneData,ID,type=type){
        if(is(zoneData)[1]=="character"){
            return("Plate Error - Check for salvaging")
        } else{
            return(extractMeasure(zoneData,ID, type=type,threshold=threshold))
        }
    },data,ID,MoreArgs=list(type=type))

    return(local)
}
