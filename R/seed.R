################################################################################
############################     lighteR      ##################################
################################################################################

## s4
setClass("roots",slot =c(Image="list",Zone = "list"))


setClass("measures",slot = c (NPQ = "data.frame",
                    XE = "data.frame",
                    OE = "data.frame",
                    EF = "data.frame"))


setClass("origin",slot = c (NPQ = "list",
                              XE = "list",
                              OE = "list",
                              EF = "list"))
setClass("parameters",slot = c (NPQ = "list",
                              XE = "list",
                              OE = "list",
                              EF = "list"))
setClass("dropped",slot = c (NPQ = "data.frame",
                              XE = "data.frame",
                              OE = "data.frame",
                              EF = "data.frame"))
setClass("retain",slot = c (NPQ = "data.frame",
                              XE = "data.frame",
                              OE = "data.frame",
                              EF = "data.frame"))
setClass("models",slot = c (NPQ = "list",
                              XE = "list",
                              OE = "list",
                              EF = "list"))

setClass("time",slot = c (time ="vector",
                          timePoints = "vector"))

setClass("originType",slot = c(origin = "vector"))









setClass("seed",
         slot = c(roots = "roots",
                  measures = "measures",
                  origin = "origin",
                  parameters = "parameters",
                  models = "models",
                  retain = "retain",
                  dropped = "dropped",
                  meta.data = "list",
                  time = "time",
                  originType = "vector")
         )



sowSeed <- function(files,mapID = NULL,type = c("zone", "image"), areaThreshold = 5){
    ## Start by build object
    seed <- new("seed")

    ## Adding roots

    roots  <- .rooting(data = files,
                       mapID = mapID,
                       type = type,
                       areaThreshold =areaThreshold)


    plateError <- .plateError(roots)

    ## building seed object
    seed@roots <- plateError$roots
    seed@meta.data <- plateError$seed.meta.data
    return(seed)
}
