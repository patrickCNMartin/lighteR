################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### splitting data - extracting individuals
################################################################################
################################################################################



.extractByID <- function(light, splitby=c("plot","pedigree","line","stem"),
                         tagID=c("plot","pedigree","line","stem"), norm = c("local","global","none")){





    # extract and clean IDs
    IDs <- as.character(light$Zone)
    IDSplit<- lapply(IDs,strsplit," ")

    IDSplit<-lapply(IDSplit,function(x, tags){return(lapply(x, .IDtag,tags))},
                    tags=tagID)
    IDSplit <- do.call("rbind",lapply(lapply(IDSplit,"[[",1),matrix,ncol=length(tagID)))
    colnames(IDSplit) <- tagID

    ## reinsert ID columns
    light <- data.frame("diskID"=light$diskID,IDSplit,light[,!colnames(light) %in% c("diskID","Zone")])
    rownames(light) <-NULL
    ## now lets split this bad boy


    plants<-split(light,lapply(splitby,function(split, data){return(data[,split])},light), drop=TRUE)




    return(plants)


}




getPlants <- function(seed,splitby=c("plot","pedigree","line","stem"),
                     tagID=c("plot","pedigree","line","stem")){

    origin <- new("origin")
    ### We already have combined them
    ### ADD split and extract function if needed
    ### It will still use the seed object just in case

    if(sum(unlist(slotApply(seed@retain,length)))>0){
        measures <- seed@retain
    } else {
        measures <- seed@measures
    }


    measures <- slotApply(measures,.extractByID,splitby=splitby,tagID=tagID)
    origin <- slotAssign(origin,measures)
    seed@origin <- origin
    seed@originType <- splitby
    return(seed)
}
