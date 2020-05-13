################################################################################
############################     lighteR      ##################################
################################################################################

## s4
setClass("roots",slot =c(Image="list",Zone = "list"))


chloroplast <- list("NPQ" = NULL,
                    "XE" = NULL,
                    "OE" = NULL,
                    "EF" = NULL)
class(chloroplast) <- "chloroplast"

leaves <- list("NPQ" = NULL,
                    "XE" = NULL,
                    "OE" = NULL,
                    "EF" = NULL)
class(leaves) <- "leaves"

parameters <- list("NPQ" = NULL,
                    "XE" = NULL,
                    "OE" = NULL,
                    "EF" = NULL)
class(parameters) <- "parameters"

models <- list("NPQ" = NULL,
               "XE" = NULL,
               "OE" = NULL,
               "EF" = NULL)
class(models) <- "models"

retain <- list("NPQ" = NULL,
               "XE" = NULL,
               "OE" = NULL,
               "EF" = NULL)
class(retain) <- "retain"

dropped <- list("NPQ" = NULL,
               "XE" = NULL,
               "OE" = NULL,
               "EF" = NULL)
class(dropped) <- "dropped"






setClass("seed",
         slot = c(roots = "roots",
                  chloroplast = "chloroplast",
                  leaves = "leaves",
                  parameters = "parameters",
                  models = "models",
                  retain = "retain",
                  dropped = "dropped",
                  meta.data = "data.frame",
                  time = "vector")
         )



sowSeed <- function(files,mapID = NULL,type = c("zone", "image"), areaThreshold = 5){
    ## Start by build object
    seed <- new("seed")

    ## Adding roots

    seed@roots <- .rooting(data = files,
                           mapID = mapID,
                           type = type,
                           areaThreshold =areaThreshold)
    return(seed)
}
