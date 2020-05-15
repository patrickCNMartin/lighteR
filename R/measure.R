################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Measure Extraction
################################################################################
################################################################################

.extractMeasure <- function(data,ID,type=c("NPQ","XE","EF","OE"),threshold=5){

    datasub <- nonZeroIndex(data,threshold)


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
    class(datasubSplit)<-"byMeasure"
    return(datasubSplit)

}

## extracting evey measure in batches for zone data

.batchExtractMeasure <- function(data,type=c("NPQ","XE","EF","OE"),threshold=5){

    zoneOnly<-data$zone
    ID<-names(zoneOnly)
    local<-mapply(function(zoneData,ID,type=type){
        if(is(zoneData)[1]=="character"){
            return("Plate Error - Check for salvaging")
        } else{
            return(extractMeasure(zoneData,ID, type=type,threshold=threshold))
        }
    },zoneOnly,ID,MoreArgs=list(type=type))
    class(local)<-"byMeasure"
    return(local)
}
