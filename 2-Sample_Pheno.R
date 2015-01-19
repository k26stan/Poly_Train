## 2-Sample_Pheno.R ##
## R function to Split Cohort into Training & Test Sets ##
## Takes in Phenotype File, spits out 2 Phenotype Files ##
## November 21, 2014 ##
## Kristopher Standish ##

## Usage ##
# 2-Sample_Pheno.R <Path/To/Pheno.txt> <New/Path/To/Pheno_>

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

## Pull in Command Line Arguments
LINE <- commandArgs(trailingOnly = TRUE)

## Parse Command Line Arguments
PathToPheno <- LINE[1]
PathToOut <- LINE[2]

###############################################################
## LOAD DATA ##################################################
###############################################################

## Specify Paths to Training/Test File Outputs
PathToTrain <- paste(PathToOut,"Train.txt",sep="_")
PathToTest <- paste(PathToOut,"Test.txt",sep="_")

## Load Table
PH <- read.table(PathToPheno,header=T)

###############################################################
## SAMPLE TEST/TRAINING SETS ##################################
###############################################################

## Sample Test/Training Sets
TRAIN_PRC <- 0.8
TRAIN <- sample( 1:nrow(PH), TRAIN_PRC*nrow(PH), replace=F )
PH.TRAIN <- PH[ TRAIN, ]
PH.TEST <- PH[ -TRAIN, ]

###############################################################
## WRITE OUTPUT TABLES ########################################
###############################################################

## Write Table
write.table(PH.TEST,PathToTest,sep="\t",row.names=F,col.names=T,quote=F)
write.table(PH.TRAIN,PathToTrain,sep="\t",row.names=F,col.names=T,quote=F)

###############################################################
## END OF DOC #################################################
###############################################################
