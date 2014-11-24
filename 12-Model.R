## 12-Model.R ##
## R function to Use Training Results and Build Model ##
## Load Phenotype & Covariate Files
## Takes in SNP list (_ID.list), *compile.short, and Genotypes (*_012.raw) ##
## Puts together Predictive Model ##
## Tests it on "Test" set ##
## Calculate R2 after fitting model ## ...I guess.
## November 21, 2014 ##
## Kristopher Standish ##

## Usage ##
# 12-Model.R <New/Path/To/Pheno_> <New/Path/To/Cov.file> <Path/To/*compile.short> <Path/To/Genotypes/*_012.raw>

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

## Pull in Command Line Arguments
LINE <- commandArgs(trailingOnly = TRUE)

## Parse Command Line Arguments
PathToPheno <- LINE[1]
PathToCovs <- LINE[2]
PathToComp <- LINE[3]
PathToGeno <- LINE[4]

###############################################################
## LOAD DATA ##################################################
###############################################################

## Set Path to Output
PathToOut <- gsub("Cov_w_PCs.txt","MOD",PathToCov)

## Specify Paths to Training/Test File Outputs
PathToTrain <- paste(PathToPheno,"Train.txt",sep="_")
PathToTest <- paste(PathToPheno,"Test.txt",sep="_")

## Load Phenotype Files
PH.TR <- read.table( PathToTrain, sep="\t", header=T)
PH.TS <- read.table( PathToTest, sep="\t", header=T)

## Load Covariate File
COV <- read.table( PathToCov, sep="\t", header=T )

## Load Compiled File
CP <- read.table( PathToComp, sep="\t", header=T, quote="", comment.char="" )

## Load Genotype File
GT <- read.table( PathToGeno, header=T )

###############################################################
## PULL OUT BEST VARIANTS #####################################
###############################################################

## Remove HW Violations
THRSH <- 1e-8
RM.HWE <- which( CP$P_HW < THRSH )
CP.2 <- CP[ -RM.HWE, ]

## Pull out Best Remaining P-Values
BEST.COUNT <- 20
KP.P <- order(CP.2$P_Assoc)[1:BEST.COUNT]
KP.P.nm <- CP.2$SNP[KP.P]
CP.3 <- CP.2[CP.3,]

###############################################################
## MERGE FILES ################################################
###############################################################

## Merge Phenotype w/ Covs
MG.TR.1 <- merge( x=PH.TR, y=COV[c(1,3:ncol(COV))], by="FID" )
MG.TS.1 <- merge( x=PH.TS, y=COV[c(1,3:ncol(COV))], by="FID" )

## And Merge w/ Genotypes
MG.TR.2 <- merge( x=MG.TR.1, y=GT[c(1,7:ncol(GT))], by="FID" )
MG.TS.2 <- merge( x=MG.TS.1, y=GT[c(1,7:ncol(GT))], by="FID" )

## And get rid of IID & FID in table
MG.TR <- MG.TR.2[,3:ncol(MG.TR.2)]
MG.TS <- MG.TS.2[,3:ncol(MG.TS.2)]


### ***** CAN'T HAVE MORE PREDICTORS THAN SAMPLES ****** ###

###############################################################
## END OF DOC #################################################
###############################################################
