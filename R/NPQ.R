################################################################################
############################ NPQ analysis ######################################
################################################################################


runNPQAnalysis <-function(map_ID,input_data,output_data,plots,
                          LoadingType,threshold,splitby,combine,norm,
                          filterStep,filterMethod,filterthreshold,
                          modelType,time,individual,convertTime,
                          summaryOutPut,cores){

      ### Loading data
      message("Loading Data")
      #browser()
      loadedData <- batchLoading(input_data,map_ID,threshold=threshold,type=LoadingType)

      ## splitting by measure
      ## For every file -> split data by each measure type (NPQ,XE,EF,OE)
      message("Split by Measure")
      byMeasure <- batchExtractMeasure(loadedData)

      ## splitting by type
      ## In this case we are combining all data sets/files together
      ## Combined by measure as they have been split abive
      ## Then you split again for each split by (pedigree, line etc )
      ## This returns a list for all measures types and within each list
      ## You have the pedigree/ line split
      ## This takes into account multiple occurance of the same pedigree/line
      message("Split by ID")
      byID <- batchExtractByID(byMeasure, splitby=splitby, combine=combine,norm=norm)

      if(filterStep=="both"){
          message("2 step Filtering")
          ## This tell you which one to keep and which one to remove
          ## returns T/F for each disk
          byIDselected <- selectDisks(byID)

          ## Retains good disks
          byIDFiltered <- filterDisks(byID,byIDselected,keep="retain")

          ## Dropped disks after first step filtering
          preModDrop<-filterDisks(byID,byIDselected,keep="drop")

          message("step 1 : DONE")
          fitLocal <- batchModelFit(byIDFiltered,model=modelType,time=time,individual=F,cores=cores)
#browser()
          PostModSelect <- selectDisks(byIDFiltered,fitLocal,time=c(40,80),
                                       method=c(filterMethod,filterthreshold))

          postModDrop <- filterDisks(byIDFiltered,PostModSelect,keep="drop")
          byIDFiltered <-filterDisks(byIDFiltered,PostModSelect,keep="retain")

          message("step 2 : DONE")
          message("Merging dropped")
          dropped <-mergeFilteredDisks(preModDrop,postModDrop)

      }else if(filterStep=="preModel"){
          message("pre Model Filter")
          byIDselected <- selectDisks(byID)
          byIDFiltered <- filterDisks(byID,byIDselected,keep="retain")
          dropped<-filterDisks(byID,byIDselected,keep="drop")
          message("Filtering: DONE")
      } else if(filterStep=="postModel"){
        message("postModel Filtering")
        fitLocal <- batchModelFit(byID,model=modelType,time=time,individual=F,cores=cores)

        PostModSelect <- selectDisks(byID,fitLocal,time=c(40,80),
                                     method=c(filterMethod,filterthreshold))
        byIDFiltered <-filterDisks(byID,PostModSelect,keep="retain")
        dropped <- filterDisks(byID,PostModSelect,keep="drop")
        message("Filtering Done")

      } else if(filterStep=="none"){
           message("No Filter Applied")
           byIDFiltered <- byID
           dropped<-"NONE"
      } else{
        stop("Filter Step should be one of the following: preModel,postModel,both,none")
      }
browser()
      ### Model Fit on all data of filtered data
      message("Fitting Models to data")
      fit <-batchModelFit(byIDFiltered,model=modelType,time=time,individual=individual,cores=cores)

      ### Plotting

      ##setting up file names
      if(length(modelType)==2){
         modelType<-paste(modelType[1],"_",modelType[2])
      } else{
        modelType<-modelType[1]
      }

      if(combine){
         combine<-"combined"
      }else{
         combine <-"individual"
      }

      filename<-paste0(plots,modelType,"_",filterStep,"_",filterMethod,"_",combine,"_")
      message("Plotting Data and fitted curves")
      plotScoresOverTimePerZone(byIDFiltered,model=fit,file=filename, split="measure", time=time)
      if(dropped!="NONE"){
        message("Plotting Dropped Data")
        filename<-paste0(plots,modelType,"_",filterStep,"_",filterMethod,"_",combine,"_dropped_")
        plotScoresOverTimePerZone(dropped,model=droppedPostFilter,file=filename, split="measure", time=time)

      }

      ## preparinf for exporting
      if(dropped!="NONE"){
          message("Extracting and merging Parameters from dropped data")
          dropped<-.extractAndMergeParam(dropped)
      } else {
         dropped <-"NONE"
      }

      message("Extracting and merging Parameters")
      dataFiltered <- .extractAndMergeParam(byIDFiltered,fit)

      if(convertTime){
        message("converting Time")
         dataFiltered<-.convertTime(dataFiltered,loadedData)
      }


exportData(dropped,file=paste0(output_data,"DroppedDisks_"))

      if(summaryOutPut=="individualOnly"){
        message("Exporting Individual Files")
        fileRetDat<-paste0(output_data,"retainedData_")
        fileRetMod<-paste0(output_data,"retainedModel_")
        exportData(dataFiltered[[1]],file=fileRetDat)
        exportData(dataFiltered[[2]],file=fileRetMod)
        if(dropped!="NONE"){
            file <- paste0(output_data,"DroppedData_")
            exportData(dropped,file=file)
        }
      } else if(summaryOutPut=="summaryOnly"){
        message("Exporting Summary Files")
        fileRetDat<-paste0(output_data,"retainedData_")
        fileRetMod<-paste0(output_data,"retainedModel_")
        generateSummaryFiles(dataFiltered[[1]],file=fileRetDat)
        generateSummaryFiles(dataFiltered[[2]],file=fileRetMod)
        if(dropped!="NONE"){
          file <- paste0(output_data,"DroppedData_")
          generateSummaryFiles(dropped,file=file)
        }
      } else if(summaryOutPut=="both"){
        message("Exporting Individula and Summary Files")
        fileRetDat<-paste0(output_data,"retainedData_")
        fileRetMod<-paste0(output_data,"retainedModel_")
        exportData(dataFiltered[[1]],file=fileRetDat)
        exportData(dataFiltered[[2]],file=fileRetMod)
        generateSummaryFiles(dataFiltered[[1]],file=fileRetDat)
        generateSummaryFiles(dataFiltered[[2]],file=fileRetMod)
        if(dropped!="NONE"){
          file <- paste0(output_data,"DroppedData_")
          generateSummaryFiles(dropped,file=file)
          exportData(dropped,file=file)
        }
      }
      message("Analysis Completed")

}
