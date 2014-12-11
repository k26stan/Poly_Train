## 12-Model_B.R ##
## R function to Use Training Results and Build Model ##
## BINARY Phenotypes ##
## Load Phenotype & Covariate Files
## Takes in SNP list (_ID.list), *compile.short, and Genotypes (*_012.raw) ##
## Puts together Predictive Model ##
## Tests it on "Test" set ##
## Calculate R2 after fitting model ## ...I guess.
## November 21, 2014 ##
## Kristopher Standish ##

## Usage ##
# Rscript 12-Model.R <New/Path/To/Pheno_> <New/Path/To/Cov.file> <Path/To/*compile.short> <Path/To/Genotypes/*_012.raw> <COV_1,COV_2,COV_3>
# Rscript 12-Model.R ${PHENO_SETS} ${NEW_COV_FILE} ${CND_FILE%%txt}.compile.short ${CND_012}.raw ${COVS_COMMAND}

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

## Pull in Command Line Arguments
LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( paste( "/projects/janssen/Poly_Train/20141124test_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/", c("LT8_DEL_MNe_MN", "Cov_w_PCs.txt", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2_012.raw"), sep="" ), "DAS_BL_MN,PC1,PC2" )
# LINE <- c( paste( "/Users/kstandis/Downloads/20141124_Poly_Train/", c("LT8_DEL_MNe_MN", "Cov_w_PCs.txt", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2_012.raw"), sep="" ), "DAS_BL_MN,PC1,PC2" )

## Parse Command Line Arguments
PathToPheno <- LINE[1]
PathToCovs <- LINE[2]
PathToComp <- LINE[3]
PathToGeno <- LINE[4]
Covs_Command <- LINE[5]

###############################################################
## LOAD DATA ##################################################
###############################################################

## Separate Covariates into Vector
Cov_Cols <- strsplit( Covs_Command, "," )[[1]]

## Set Path to Output
PathToOut <- gsub("Cov_w_PCs.txt","MOD",PathToCovs)

## Specify Paths to Training/Test File Outputs
PathToTrain <- paste(PathToPheno,"Train.txt",sep="_")
PathToTest <- paste(PathToPheno,"Test.txt",sep="_")

## Load Phenotype Files
PH.TR <- read.table( PathToTrain, sep="\t", header=T)
PH.TS <- read.table( PathToTest, sep="\t", header=T)

## Load Covariate File
COV <- read.table( PathToCovs, sep="\t", header=T )

## Load Compiled File
CP <- read.table( PathToComp, sep="\t", header=T, quote="", comment.char="" )

## Load Genotype File
GT <- read.table( PathToGeno, header=T )

###############################################################
## PULL OUT BEST VARIANTS #####################################
###############################################################

## Remove Duplicate SNPs
RM.DUP <- which(duplicated( CP$SNP ))
CP.1 <- CP[ -RM.DUP, ]

## Remove HW Violations
THRSH.HWE <- 1e-8
RM.HWE <- which( CP.1$P_HW < THRSH.HWE )
CP.2 <- CP.1[ -RM.HWE, ]

## Pull out Variants w/ MAF > 10%
THRSH.MAF <- 0.05
RM.MAF <- which( CP.2$REF_ALL < THRSH.MAF | CP.2$ALT_ALL < THRSH.MAF )
CP.3 <- CP.2[ -RM.MAF, ]

## Pull out Best Remaining P-Values (from Compiled Data)
BEST.COUNT <- nrow(CP.3) # 500
KP.P <- order(CP.3$P_Assoc)[1:BEST.COUNT]
CP.4 <- CP.3[ KP.P, ]
CP.4 <- CP.4[which(!is.na(CP.4$P_Assoc)),]

## Pull out Best Variants (from Genotype File)
 # Rename Columns of GT array
GT.nms <- gsub( "_A","",colnames(GT) )# lapply( strsplit( colnames(GT), "_" ), "[", 1 )
GT.nms <- gsub( "_C","",GT.nms )
GT.nms <- gsub( "_G","",GT.nms )
GT.nms <- gsub( "_T","",GT.nms )
GT.nms <- gsub( "X","",GT.nms )
colnames(GT) <- GT.nms
 # Pull out Variants of Interest
KP.P.id <- gsub( ":",".", CP.4$SNP, fixed=T )
KP.P.id <- gsub( "X","",KP.P.id )
GT.which.a <- match( KP.P.id, GT.nms ) # which( GT.nms %in% KP.P.id )
# GT.which.2 <- match( GT.nms, KP.P.id ) # which( GT.nms %in% KP.P.id )
# GT.which.nms <- GT.nms[ GT.which.a ]
GT.cand.a <- data.matrix( GT[,GT.which.a] )
rownames(GT.cand.a) <- as.character( GT$FID )
CP.4 <- CP.4[ which(KP.P.id %in% GT.nms[GT.which.a] ), ]
KP.P.id <- gsub( ":",".", CP.4$SNP, fixed=T )
KP.P.id <- gsub( "X","",KP.P.id )

## Get Beta values for each Variant that made the cut
BETA.gt.a <- CP.4$BETA ; names(BETA.gt.a) <- KP.P.id
P.gt.a <- CP.4$P_Assoc ; names(P.gt.a) <- KP.P.id

## Remove Variants w/ Correlated Genotypes
 # Correlation > .9
GT.corr.a <- cor( GT.cand.a )
THRSH.corr <- 0.9
RM.corr <- c()
in_loop <- 0
for ( v in 1:length(P.gt.a) ) {
	var <- names(P.gt.a)[v]
	row <- which( rownames(GT.corr.a)==var )
	if ( !(var %in% RM.corr) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,":",v,"-",var,"######" ))
		print( range( (row+1):ncol(GT.corr.a) ) )
		RM.temp.a <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.a, RM.temp.2 ) )
		print(paste("",length(P.gt.a)-length(RM.corr),"Remaining" ))
	}
}
BETA.gt.9 <- BETA.gt.a[ -which( names(BETA.gt.a) %in% RM.corr ) ] #; names(BETA.gt) <- KP.P.id[ -which( names(BETA.gt.a) %in% RM.corr ) ]
P.gt.9 <- P.gt.a[ -which( names(P.gt.a) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.a) %in% RM.corr ) ]
GT.cand.9 <- GT.cand.a[ ,-which( colnames(GT.cand.a) %in% RM.corr ) ]
GT.corr.9 <- cor( GT.cand.9 )
 # Correlation > .1
GT.corr.a <- cor( GT.cand.a )
THRSH.corr <- 0.1
RM.corr <- c()
in_loop <- 0
for ( v in 1:length(P.gt.a) ) {
	var <- names(P.gt.a)[v]
	row <- which( rownames(GT.corr.a)==var )
	if ( !(var %in% RM.corr) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,":",v,"-",var,"######" ))
		print( range( (row+1):ncol(GT.corr.a) ) )
		RM.temp.a <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.a, RM.temp.2 ) )
		print(paste("",length(P.gt.a)-length(RM.corr),"Remaining" ))
	}
}
BETA.gt.1 <- BETA.gt.a[ -which( names(BETA.gt.a) %in% RM.corr ) ] #; names(BETA.gt) <- KP.P.id[ -which( names(BETA.gt.a) %in% RM.corr ) ]
P.gt.1 <- P.gt.a[ -which( names(P.gt.a) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.a) %in% RM.corr ) ]
GT.cand.1 <- GT.cand.a[ ,-which( colnames(GT.cand.a) %in% RM.corr ) ]
GT.corr.1 <- cor( GT.cand.1 )
## Compare
par(mfrow=c(1,3))
hist( GT.corr.a, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )
hist( GT.corr.9, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )
hist( GT.corr.1, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )

###############################################################
## MERGE FILES ################################################
###############################################################

## Merge Phenotype w/ Covs
MG.TR.1 <- merge( x=PH.TR, y=COV[,c(1,3:ncol(COV))], by="FID" )
MG.TS.1 <- merge( x=PH.TS, y=COV[,c(1,3:ncol(COV))], by="FID" )

## And Merge w/ Genotypes
MG.TR.2.9 <- merge( x=MG.TR.1, y=GT.cand.9, by.x="FID", by.y="row.names" )
MG.TS.2.9 <- merge( x=MG.TS.1, y=GT.cand.9, by.x="FID", by.y="row.names" )
MG.TR.2.1 <- merge( x=MG.TR.1, y=GT.cand.1, by.x="FID", by.y="row.names" )
MG.TS.2.1 <- merge( x=MG.TS.1, y=GT.cand.1, by.x="FID", by.y="row.names" )

## And get rid of IID & FID in table
MG.TR.9 <- MG.TR.2.9[,3:ncol(MG.TR.2.9)] ; rownames(MG.TR.9) <- as.character( MG.TR.2.9[,1] )
MG.TS.9 <- MG.TS.2.9[,3:ncol(MG.TS.2.9)] ; rownames(MG.TS.9) <- as.character( MG.TS.2.9[,1] )
TEMP.TR.9 <- cor( MG.TR.9[87:ncol(MG.TR.9)] )
TEMP.TS.9 <- cor( MG.TS.9[87:ncol(MG.TS.9)] )
heatmap.2( TEMP.TR.9, scale="none", trace="none")
MG.TR.1 <- MG.TR.2.1[,3:ncol(MG.TR.2.1)] ; rownames(MG.TR.1) <- as.character( MG.TR.2.1[,1] )
MG.TS.1 <- MG.TS.2.1[,3:ncol(MG.TS.2.1)] ; rownames(MG.TS.1) <- as.character( MG.TS.2.1[,1] )
TEMP.TR.1 <- cor( MG.TR.1[87:ncol(MG.TR.1)] )
TEMP.TS.1 <- cor( MG.TS.1[87:ncol(MG.TS.1)] )
heatmap.2( TEMP.TR.1, scale="none", trace="none")

