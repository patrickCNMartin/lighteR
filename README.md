# lighteR

The lighteR packages enables the analysis of Non-Photochemical quenching data as given by FluorImage. The analysis pipeline will go from loading data, cleaning data, filtering low quality disks, fitting models, plotting, extracting parameters and exporting results.

As of now, the analysis requires certain manual cleaning and formatting steps. These steps are described bellow.

## Installing lighteR

To install lighteR 


```

install.packages("devtools")
library(devtools)
install_github("patrickCNMartin/lighteR")


```
**IMPORTANT** lighteR is a "proto package" meaning that it has not been tested and will contain bugs when used by other users. Manual pages and user manual are sparse but will be spruced up in the future. See NEWS 

## Data Formating

### map ID
Map IDs contain the ID of each disk for a given plate. Ideally, mapID data should be stored in a .csv file 

The naming convention for map ID files in the following:

* date(ddmmyy)_plateNumber_mapID.csv

The purpose of this strict naming convention is to automatically assign the correct map to your input data.


### Input data
FluorImage exports many files but for this analysis, only two  are required: Zone data and Image Data.
The Zone data file contains all measures for each disk over all time points.
The Image data file contains a summary of the Zone data. Image data is optional
but required if you wish to convert time points to real time in seconds (TO BE ADDED)


As described above, a strict naming convention is required in order to map your ID's back to the data.

The naming convention for input_data is as follows:

* ZoneData = date(ddmmyy)_plateNumber_AllImages_ZoneData.csv
* ImageData = date(ddmmyy)_plateNumber_AllImages_ImageData.csv


## NEWS

* Bug Fixing 
* Increase robustness 
* Manual page re write 
* Vignette and user guide 
* Select and navigate functions 
* Goodness of fit threshold selection function 



## Example code 
```{r, eval =F}
dats <- "/Path/to/data/folder"
maps <- "/Path/to/mapID/"

seed <- sowSeed(dats, maps)
seed <- getMeasure(seed)
seed <- getOrigin(seed, splitby = c("pedigree","line"))
seed <- selectPlants(seed)
seed <- getTraits(seed)
seedMod <- modelPlants(seed, models = list("NPQ"=c("beta",3,1),"OE"=c(2,3)),
                    fit.to = "allPlants",cores=10)
seedSel <- selectPlants(seedMod,measure = c("NPQ","OE"), method ="RSD",threshold = 0.08)
seed <- getTraits(seedSel)

seedExp <- convertValues(seed)
seedExp <- convertTime(seedExp)

exportSeed(seedExp, file="../../NIAB",extension = ".csv")

pdf("../../NIAB_NPQ_NPQTriple_beta_3_1.pdf",width=12, height =12)
par(mfrow = c(3,3))
plotSeed(seed, measure = c("NPQ"))
dev.off()

pdf("../../NIAB_OE_2_3.pdf",width=12, height =12)
par(mfrow = c(3,3))
plotSeed(seed, measure = c("OE"))
dev.off()

pdf("../../NIAB_dropped.pdf",width=12, height =12)
par(mfrow = c(3,3))
plotSeed(seed, measure = c("NPQ","OE"), dropped=TRUE)
dev.off()




``
