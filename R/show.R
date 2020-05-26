################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Show Methods
################################################################################
################################################################################
setMethod("show",
          signature = "seed",
          definition = function(object){
              widthDisplay<-round(options()$width*0.75)
              cat(rep("*",15),"\n")
              cat("Seed Object \n")
              cat(rep("*",15),"\n")

              RootLengthZone <- length(object@roots@Zone)
              RootLengthImage <- length(object@roots@Image)
              cat("@roots \n")
              cat(paste(RootLengthZone,"Zone data file(s) stored in roots \n"))
              cat(paste(RootLengthImage,"Image data file(s) stored in roots \n"))

              measures <- .slotApply(object@measures, nrow)

              if(sum(unlist(measures))!=0){
                cat("@measures \n")
                cat(paste(measures["NPQ"], "sample(s) for NPQ \n"))
                cat(paste(measures["XE"], "sample(s) for XE \n"))
                cat(paste(measures["OE"], "sample(s) for OE \n"))
                cat(paste(measures["EF"], "sample(s) for EF \n"))
              }
              #origin <- slotApply(object@origin, function(origin){lapply(origin,length)})
              origin <- .slotApply(object@origin, length)

              if(sum(unlist(origin))>=1){
                  cat("@origins \n")
                  cat(paste(origin[["NPQ"]], "sample Origin(s) for NPQ \n"))
                  cat(paste(origin[["XE"]], "sample Origin(s) for XE \n"))
                  cat(paste(origin[["OE"]], "sample Origin(s) for OE \n"))
                  cat(paste(origin[["EF"]], "sample Origin(s) for EF \n"))
              }

              traits <- .slotApply(object@traits, nrow)

              if(sum(unlist(traits))>=1){
                  cat("@traits \n")
                  cat(paste("Traits extracted from",traits[["NPQ"]], "sample(s) in NPQ \n"))
                  cat(paste("Traits extracted from",traits[["XE"]], "sample(s) in XE \n"))
                  cat(paste("Traits extracted from",traits[["OE"]], "sample(s) in OE \n"))
                  cat(paste("Traits extracted from",traits[["EF"]], "sample(s) in EF \n"))
              }

              models <- .slotApply(object@models, length)
              modelType<- object@meta.param@models


              if(sum(unlist(models))>=1){
                  cat("@models \n")
                  cat(paste(paste(modelType[["NPQ"]],collapse = " "),"models ran on NPQ data \n"))
                  cat(paste(paste(modelType[["XE"]],collapse = " "),"models ran on XE data \n"))
                  cat(paste(paste(modelType[["OE"]],collapse = " "),"models ran on OE data \n"))
                  cat(paste(paste(modelType[["EF"]],collapse = " "),"models ran on EF data \n"))

              }
             ### meta data info
              if(length(grep("Zone",names(object@meta.data))) ==1 |
                 length(grepl("Image",names(object@meta.data)))==1){
                  cat("@meta.data \n")
                  cat("Plate Errors stored in meta.data \n")
              }
              if(any(grepl("dropped",names(object@meta.data)))){
                  cat(paste(object@meta.data$dropped["dropped"],"filtered sample(s) \n"))
                  cat(paste(object@meta.data$dropped["retain"],"retained sample(s) \n"))
              }



          })

