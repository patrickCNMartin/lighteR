################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Loading Function
################################################################################
################################################################################


### Loading function
ImageDataLoading <- function(file,colnamesExceptions=c("qP","NPQ")){
    ## special load and split

    file<-read.csv(file, stringsAsFactors=F)

    ## colmun renaming
    temp1 <- colnames(file)
    temp2 <- as.character(file[1,])

    ## not great
    # update to include regex/ and exception shifting
    # see name filter function (work in progress)
    temp1[grep("X", temp1)]<- temp2[grep("X", temp1)]
    colnames(file) <- temp1

    ## For now

    ## splitting for further Analysis
    ## not necessary, it's just cleaner
    ## and I like lists
    if(any(colnames(file)=="Type")){
        file <- split(file,file$Type)
    } else if(any(grepl("NPQ",file[,4]))){
        file <- split(file,file[,4])
    } else {
       stop("Woops! Seems like I can't read the file you are giving me...")
    }

    return(file)
}

ZoneDataLoading <- function(file, mapID=NULL,threshold=5){
    ## load data
    data<-read.csv(file, header=F, skip=3,stringsAsFactors=F)

    ## loading colunm names seperately
    ## there is a weird split going on so this might be easier

    columns <- c("Zone",as.vector(as.matrix(data[1,3:ncol(data)])))
    columns<- gsub("\xb2","",columns)
    columns <- gsub(" ","_",columns)

    ## data cleaning
    data <- as.data.frame(data[2:nrow(data),2:ncol(data)])
    data<- data.frame(data[,1], apply(data[,2:ncol(data)],2, as.numeric))
    colnames(data) <- columns


    if(!is.null(mapID)){
          mapID<-read.csv(mapID, header=F, stringsAsFactors=F)
          data <- mapIDs(data,mapID,threshold)
    }

    return(data)

}

batchLoading<-function(data,mapID=NULL,type=c("zone","image","all"),threshold=5){
    ## this is a bit wonky
    ## loading data from folder
    dataLocal <- paste0(data,dir(data))
    mapIDLocal <- paste0(mapID, dir(mapID))

    ## loading data files
    if(type[1]=="all"){

         ## Finding Zone data
         zoneFiles<-grep("ZoneData",dataLocal, value=T)

         ## reorganising maps if needed
         if(!is.null(mapID)){
             mapID<-lapply(zoneFiles,matchMapToData,map=mapIDLocal,dataLoc=data)

         } else{
            mapID<-vector("list", length(data))
         }

         zone<-mapply(function(zoneFiles,mapID,threshold){
                    mapID<-nullConversion(mapID)

                    zoneData <-tryCatch(ZoneDataLoading(zoneFiles,mapID,threshold),
                    error=function(cond){return("Plate Error - Check for salavaging \n")},
                    warning=function(cond){return("Plate Error - Check for salavaging \n")})
                    return(zoneData)},zoneFiles=zoneFiles,mapID=mapID,threshold=threshold,SIMPLIFY=F)

         names(zone)<-gsub(".csv","",zoneFiles)


        ## loading image files
        ## using setwd just to avoid any issues with differences between operating systems
        ## it's ugly but it should work

        imageFiles<-grep("ImageData",dataLocal, value=T)
        image<-lapply(imageFiles,ImageDataLoading)
        names(image)<-gsub(".csv","",imageFiles)

        localDat<-list("Image"=image,"zone"=zone)


    } else if(type[1]=="zone"){
      ## reorganising maps if needed
      ## Finding Zone data

      zoneFiles<-grep("ZoneData",dataLocal, value=T)

      ## reorganising maps if needed
      if(!is.null(mapID)){
          mapID<-lapply(zoneFiles,matchMapToData,map=mapIDLocal,dataLoc=data)

      } else{
         ## keeping it is list so we can parse NULL elements to the next Functions
         ## otherwise when using vectors ,R removes those elements
         mapID<-vector("list", length(data))
      }

      zone<-mapply(function(zoneFiles,mapID,threshold){
                 mapID<-nullConversion(mapID)
                 print(zoneFiles)
                 zoneData <-tryCatch(ZoneDataLoading(zoneFiles,mapID,threshold),
                 error=function(cond){return("Plate Error - Check for salvaging \n")},
                 warning=function(cond){return("Plate Error - Check for salvaging \n")})
                 return(zoneData)},zoneFiles=zoneFiles,mapID=mapID,threshold=threshold,SIMPLIFY=F)

      names(zone)<-gsub(".csv","",zoneFiles)

      localDat<-list("zone"=zone)

    } else if(type[1]=="image"){

      imageFiles<-grep("ImageData",dataLocal, value=T)
      image<-lapply(imageFiles,ImageDataLoading)
      names(image)<-gsub(".csv","",imageFiles)
      localDat<-list("image"=image)

    }else{
      stop("Unknown input data type")
    }

    return(localDat)

}



#### new map ID function hopefully this one will work

mapIDs <- function(data,mapID,threshold=5){

    mapID <-as.matrix(mapID)
    ## finding map bounds
    boundsX <- apply(mapID,2,function(map){
                             loc<-grepl("([0-9]+).*$",map)
                             if(!all(loc==FALSE)) loc <- TRUE
                             if(all(loc==FALSE)) loc<- FALSE
                             return(loc)
                             })
    boundsY<-apply(mapID,1,function(map){
                              loc<-grepl("([0-9]+).*$",map)
                              if(!all(loc==FALSE)) loc <- TRUE
                              if(all(loc==FALSE)) loc<- FALSE
                              return(loc)
                              })
    ## redcuing map

    mapID<-mapID[boundsY,boundsX]
    locID <- grep("([0-9]+).*$",mapID)
    data <- nonZeroIndex(data,threshold)
    data[,"Zone"]<- as.vector(t(mapID[locID]))
    return(data)

}



matchMapToData<-function(data,map,dataLoc){
     ## match would have been more elgant
     ## but match is only for exact patterns
     ## and there will be issues with OS differences
     ## so good old school looping for now
     tmp<-gsub(dataLoc,"",data)
     tmp<-sapply(strsplit(tmp,"_"),function(x){return(paste0(x[1],"_",x[2]))})
     mapMatch<-grep(tmp,map,value=T)
     if(length(mapMatch)==0){
        return("NULL")
     } else{
        return(mapMatch)
     }
}
