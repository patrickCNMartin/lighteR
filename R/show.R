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

              measures <- slotApply(object@measures, nrow)

              if(sum(unlist(measures))!=0){
                cat("@measures \n")
                cat(paste(measures["NPQ"], "sample(s) for NPQ \n"))
                cat(paste(measures["XE"], "sample(s) for XE \n"))
                cat(paste(measures["OE"], "sample(s) for OE \n"))
                cat(paste(measures["EF"], "sample(s) for EF \n"))
              }
              if(length(grep("Zone",names(object@meta.data))) ==1 |
                 length(grepl("Image",names(object@meta.data)))==1){
                  cat("@meta.data \n")
                  cat("Plate Errors stored in meta.data")
              }



          })

