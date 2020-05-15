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

              if(length(grep("Zone",names(object@meta.data))) ==1 |
                 length(grepl("Image",names(object@meta.data)))==1){
                  cat("@meta.data \n")
                  cat("Plate Errors stored in meta.data")
              }

          })

