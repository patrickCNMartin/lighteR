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
setClass("traits",slot = c (NPQ = "data.frame",
                              XE = "data.frame",
                              OE = "data.frame",
                              EF = "data.frame"))
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
setClass("retained.models",slot = c (NPQ = "list",
                            XE = "list",
                            OE = "list",
                            EF = "list"))

setClass("meta.param",slot = c (time ="vector",
                          timePoints = "vector",
                          originType ="vector",
                          models = "list",
                          nlsrStart = "list",
                          fitto = "vector"),
                      prototype = prototype(nlsrStart = list(a=10, b=7.5, c=5),
                                            models = list("NPQ" = "No",
                                                          "XE" = "No",
                                                          "OE" = "No",
                                                          "EF" = "No")))











setClass("seed",
         slot = c(roots = "roots",
                  measures = "measures",
                  origin = "origin",
                  traits = "traits",
                  models = "models",
                  retain = "retain",
                  dropped = "dropped",
                  retained.models = "retained.models",
                  meta.data = "list",
                  meta.param = "meta.param")
         )


#' Create Seed object for lighteR analysis.
#'
#' @param files A directory containing Zone files and Image files
#' @param mapID A directory containing plant IDs
#' @param type Type of data that should be loaded - either zone and/or image
#' @param areaThreshold Area in square milimeters that should be considered as noise if appearing in grid
#' @return Returns a new seed object - contains data provided in \code{files} and \code{mapID}

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
