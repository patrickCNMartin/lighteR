################################################################################
############################ NPQ analysis ######################################
################################################################################


################################################################################
################################################################################
### Export and formatting
################################################################################
################################################################################

#' Export data from seed object
#'
#' @param seed a seed object
#' @param file filename to be used - note that this should only be the first part of the file name as extenssion will be added seperately
#' @param dataType character string describing which data to export. "retain","dropped","traits","measures"
#' @param extension file extension to be used
#' @param sep data seperator
#' @return Create files in specified directory containing extracted data.
exportSeed<- function(seed,file=NULL,dataType = c("retain","dropped","traits","measures"),extension =".csv",sep=","){


      ## Just in case
  ### This makes sense
      if(is.null(file)){

        file <-paste0("SeedExport",gsub(" ","_", Sys.time()),extension)
      }

      ## extracting data sets
      for(dt in seq_along(dataType)){
          tmp <- slot(seed,dataType[dt])

          slots <- slotNames(class(tmp))
          for(sl in slots){
              m <- slot(tmp,sl)
              f <- paste0(file,"_",dataType[dt],"_",sl,extension)
              write.table(m,file=f,sep=sep,col.names = TRUE,row.names = FALSE)
          }
      }
}

#' Plot data from seed object
#'
#' @param seed a seed object
#' @param measure light measure that should be used for plotting
#' @param dropped logical indicating if you wish to plot dropped samples
#' @return plot with sample data -  may include models if they have been computed

plotSeed <- function(seed,measure=c("NPQ","EF","XE","OE"),dropped=FALSE){
    ## checkink what has been computed
    is.origin.empty <- sum(unlist(.slotApply(seed@origin, length))) == 0
    is.retain.empty <- sum(unlist(.slotApply(seed@retain, length))) == 0
    is.traits.empty <- sum(unlist(.slotApply(seed@traits, length))) == 0
    is.dropped.empty <- sum(unlist(.slotApply(seed@dropped, length))) == 0
    is.models.empty <- sum(unlist(.slotApply(seed@models, length))) == 0
    is.retained.models.empty <- sum(unlist(.slotApply(seed@retained.models, length))) == 0

    if(dropped){
      d <- seed@dropped


      for(i in seq_along(measure)){
        tmp <- slot(d,measure[i])
        if(length(tmp) ==0) next()
        tag <- tmp$Zone
        tmp <- tmp[,!colnames(tmp) %in% c("diskID","Zone")]

        for(r in seq_len(nrow(tmp))){
          plot(seq_len(ncol(tmp)), tmp[r,],xlim=c(1,ncol(tmp)),
               ylim=c(0,1),xlab ="Time Points", ylab = measure[i], main = tag[r],pch =19 , col = "#03a5fc")
        }
      }
    }else if(is.origin.empty & is.retain.empty & is.models.empty & is.retained.models.empty ){
         d <- seed@measures

         for(i in seq_along(measure)){
            tmp <- slot(d,measure[i])
            if(length(tmp) ==0) next()
            tag <- tmp$Zone
            tmp <- tmp[,!colnames(tmp) %in% c("diskID","Zone")]

            for(r in seq_len(nrow(tmp))){
                plot(seq_len(ncol(tmp)), tmp[r,],xlim=c(1,ncol(tmp)),
                     ylim=c(0,1),xlab ="Time Points", ylab = measure[i],
                     pch =19 , col = "#03a5fc", main = tag[r])
            }
         }
    } else if(is.origin.empty & !is.retain.empty & is.models.empty &is.retained.models.empty){
      d <- seed@retain
      sl <- slotNames(class(d))
      for(i in seq_along(measure)){
        tmp <- slot(d,measure[i])
        if(length(tmp) ==0) next()
        tag <- tmp$Zone
        tmp <- tmp[,!colnames(tmp) %in% c("diskID","Zone")]
        for(r in seq_len(nrow(tmp))){
          plot(seq_len(ncol(tmp)), tmp[r,],xlim=c(1,ncol(tmp)),
               ylim=c(0,1),xlab ="Time Points", ylab = measure[i], main = tag[r],
               pch =19 , col = "#03a5fc")
        }
      }
    } else if(is.origin.empty & !is.retain.empty & !is.traits.empty & !is.retained.models.empty){
      d <- seed@retain
      m <- seed@traits
      sl <- slotNames(class(d))
      for(i in seq_along(measure)){
        tmpDat <- slot(d,measure[i])
        if(length(tmpDat) ==0) next()
        tagDat <- tmp$Zone
        tmpDat <- tmp[,!colnames(tmp) %in% c("diskID","Zone")]
        tmpMod <-slot(m,measure[i])
        if(ncol(tmpMod)==16)next()
        loc <- apply(tmpMod,1,grep,"fitted")
        tmpMod <- tmpMod[,seq(loc+1, ncol(tmpMod))]
        for(r in seq_len(nrow(tmp))){
          plot(seq_len(ncol(tmp)), tmp[r,],xlim=c(1,ncol(tmp)), ylim=c(0,1),
               xlab ="Time Points", ylab = measure[i], main = tag[r],
               pch =19 , col = "#03a5fc")
          mt <-as.numeric(as.vector(as.matrix(tmpMod[j,])))
          lines(seq_len(ncol(tmpMod)),mt, col ="#d65011")
        }
      }
    } else if(is.origin.empty & is.retain.empty & !is.traits.empty & !is.models.empty){
      d <- seed@measures
      m <- seed@traits
      sl <- slotNames(class(d))
      for(i in seq_along(measure)){
        tmpDat <- slot(d,measure[i])
        if(length(tmpDat) ==0) next()
        tagDat <- tmp$Zone
        tmpDat <- tmp[,!colnames(tmp) %in% c("diskID","Zone")]
        tmpMod <-slot(m,measure[i])
        if(length(tmpMod) ==0) next()
        if(ncol(tmpMod)==16)next()
        loc <- apply(tmpMod,1,grep,"fitted")
        tmpMod <- tmpMod[,seq(loc+1, ncol(tmpMod))]
        for(r in seq_len(nrow(tmp))){
          plot(seq_len(ncol(tmp)), tmp[r,],xlim=c(1,ncol(tmp)), ylim=c(0,1),
               xlab ="Time Points", ylab = measure[i], main = tag[r],
               pch =19 , col = "#03a5fc")
          mt <-as.numeric(as.vector(as.matrix(tmpMod[j,])))
          lines(seq_len(ncol(tmpMod)),mt, col ="#d65011")
        }
      }
    } else if(!is.origin.empty & is.retain.empty & is.traits.empty & is.models.empty){
          d <- seed@origin
          sl <- slotNames(class(d))
          for(i in seq_along(measure)){
            tmp <- slot(d,measure[i])
            if(length(tmp) ==0) next()
            for(ori in seq_along(tmp)){
                datInt <- tmp[[ori]][,!colnames(tmp[[ori]]) %in% c("diskID","pedigree","line","stem")]
                x <- rep(seq(1,ncol(datInt)),nrow(datInt))
                y <- as.vector(as.matrix(t(datInt)))
                plot(x,y,xlim=c(1,ncol(datInt)), ylim=c(0,1),
                     xlab ="Time Points", ylab = measure[i], main = names(tmp)[ori],
                     pch =19 , col = "#03a5fc")
            }
          }
    } else if(!is.origin.empty & is.retain.empty & !is.traits.empty & !is.models.empty){
      d <- seed@origin
      m <- seed@traits
      sl <- slotNames(class(d))
      for(i in seq_along(measure)){
        tmp <- slot(d,measure[i])
        if(length(tmp) ==0) next()

        tmpMod <- slot(m,measure[i])
        if(length(tmpMod) ==0) next()


        for(ori in seq_along(tmp)){
          datInt <- tmp[[ori]][,!colnames(tmp[[ori]]) %in% c("diskID","pedigree","line","stem")]
          tag <- tmp[[ori]][,colnames(tmp[[ori]]) %in% c("diskID","pedigree","line","stem")]
          tag <- apply(tag,1,paste, collapse ="")
          tag <- gsub(" ","",tag)
          tag <- gsub("missing","", tag )
          x <- rep(seq(1,ncol(datInt)),nrow(datInt))
          y <- as.vector(as.matrix(t(datInt)))
          plot(x,y,xlim=c(1,ncol(datInt)), ylim=c(0,1),
               xlab ="Time Points", ylab = measure[i], main = names(tmp)[ori],
               pch =19 , col = "#03a5fc")

          traittag <- tmpMod[,colnames(tmpMod) %in% c("diskID","Zone")]
          traittag <-apply(traittag,1,paste, collapse="")
          traittag <- gsub(" ","",traittag)
          names(traittag)<- NULL
          modReorder <- tmpMod[match(tag,traittag),]
          loc <- unique(apply(modReorder,1,function(x)grep("fitted",x)))




          modInt <- tmpMod[,seq(loc+1,ncol(tmpMod))]
          for(row in seq_len(nrow(modInt))){
            mx <-seq(1,ncol(modInt))
            my <-as.numeric(as.vector(as.matrix(modInt[row,])))
            lines(mx,my,col="#d65011",lwd=1.6)
          }

        }
      }
    }else if(!is.origin.empty & !is.retain.empty & !is.traits.empty & !is.retained.models.empty){
      d <- seed@origin
      m <- seed@traits
      sl <- slotNames(class(d))
      for(i in seq_along(measure)){
        tmp <- slot(d,measure[i])
        if(length(tmp) ==0) next()
        tmpMod <- slot(m,measure[i])
        if(length(tmpMod) ==0) next()

        for(ori in seq_along(tmp)){
          datInt <- tmp[[ori]][,!colnames(tmp[[ori]]) %in% c("diskID","plot","pedigree","line","stem")]
          tag <- tmp[[ori]][,colnames(tmp[[ori]]) %in% c("diskID","plot","pedigree","line","stem")]
          tag <- apply(tag,1,paste, collapse ="")
          tag <- gsub(" ","",tag)
          tag <- gsub("missing","", tag )
          names(tag) <- NULL
          x <- rep(seq(1,ncol(datInt)),nrow(datInt))
          y <- as.vector(as.matrix(t(datInt)))
          plot(x,y,xlim=c(1,ncol(datInt)), ylim=c(0,1),
               xlab ="Time Points", ylab = measure[i], main = names(tmp)[ori],
               pch =19 , col = "#03a5fc")

          traittag <- tmpMod[,colnames(tmpMod) %in% c("diskID","Zone")]
          traittag <-apply(traittag,1,paste, collapse="")
          traittag <- gsub(" ","",traittag)
          names(traittag)<- NULL
          modReorder <- tmpMod[match(tag,traittag),]
          loc <- unique(apply(modReorder,1,function(x)grep("fitted",x)))


          modInt <- modReorder[,seq(loc+1,ncol(modReorder))]

          for(row in seq_len(nrow(modInt))){
              mx <-seq(1,ncol(modInt))
              my <-as.numeric(as.vector(as.matrix(modInt[row,])))
              lines(mx,my,col="#d65011",lwd=1.6)
          }



        }
      }
    }
}

