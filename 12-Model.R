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
THRSH.MAF <- 0.1
RM.MAF <- which( CP.2$REF_ALL < THRSH.MAF | CP.2$ALT_ALL < THRSH.MAF )
CP.3 <- CP.2[ -RM.MAF, ]

## Pull out Best Remaining P-Values (from Compiled Data)
BEST.COUNT <- 500
KP.P <- order(CP.3$P_Assoc)[1:BEST.COUNT]
CP.4 <- CP.3[ KP.P, ]

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
GT.which.1 <- which( GT.nms %in% KP.P.id )
GT.cand.1 <- data.matrix( GT[,GT.which.1] ) ; rownames(GT.cand.1) <- as.character( GT$FID )
CP.4 <- CP.4[ which(KP.P.id %in% GT.nms[GT.which.1] ), ]
KP.P.id <- gsub( ":",".", CP.4$SNP, fixed=T )

## Get Beta values for each Variant that made the cut
BETA.gt.1 <- CP.4$BETA ; names(BETA.gt.1) <- KP.P.id
P.gt.1 <- CP.4$P_Assoc ; names(P.gt.1) <- KP.P.id

## Remove Variants w/ Correlated Genotypes
GT.corr <- cor( GT.cand.1 )
THRSH.corr <- 0.2
RM.corr <- c()
in_loop <- 0
for ( v in 1:length(P.gt.1) ) {
	var <- names(P.gt.1)[v]
	row <- which( rownames(GT.corr)==var )
	if ( !(var %in% RM.corr) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,"-",var,"######" ))
		RM.temp.1 <- colnames(GT.corr)[row+which(GT.corr[var,(row+1):ncol(GT.corr)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr)[row+which(GT.corr[var,(row+1):ncol(GT.corr)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.1, RM.temp.2 ) )
		print(paste("      ",length(P.gt.1)-length(RM.corr),"Remaining" ))
	}
}
BETA.gt <- BETA.gt.1[ -which( names(BETA.gt.1) %in% RM.corr ) ] #; names(BETA.gt) <- KP.P.id[ -which( names(BETA.gt.1) %in% RM.corr ) ]
P.gt <- P.gt.1[ -which( names(P.gt.1) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.1) %in% RM.corr ) ]
GT.cand <- GT.cand.1[ ,-which( colnames(GT.cand.1) %in% RM.corr ) ]
# GT.which <- GT.which.1[ -which( names(GT.which.1) %in% RM.corr ) ]
TEMP <- cor( GT.cand )
hist( TEMP, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )

###############################################################
## MERGE FILES ################################################
###############################################################

## Merge Phenotype w/ Covs
MG.TR.1 <- merge( x=PH.TR, y=COV[,c(1,3:ncol(COV))], by="FID" )
MG.TS.1 <- merge( x=PH.TS, y=COV[,c(1,3:ncol(COV))], by="FID" )

## And Merge w/ Genotypes
MG.TR.2 <- merge( x=MG.TR.1, y=GT.cand, by.x="FID", by.y="row.names" )
MG.TS.2 <- merge( x=MG.TS.1, y=GT.cand, by.x="FID", by.y="row.names" )
# MG.TR.2 <- merge( x=MG.TR.1, y=GT[,c(1,7,GT.which)], by="FID" )
# MG.TS.2 <- merge( x=MG.TS.1, y=GT[,c(1,7,GT.which)], by="FID" )

## And get rid of IID & FID in table
MG.TR <- MG.TR.2[,3:ncol(MG.TR.2)] ; rownames(MG.TR) <- as.character( MG.TR.2[,1] )
MG.TS <- MG.TS.2[,3:ncol(MG.TS.2)] ; rownames(MG.TS) <- as.character( MG.TS.2[,1] )

TEMP.TR <- cor( MG.TR[87:ncol(MG.TR)] )
TEMP.TS <- cor( MG.TS[87:ncol(MG.TS)] )
###############################################################
## MODEL w/ COVS + GTs ########################################
###############################################################

## To-Do List
#1) Rank Predictors
#2) Calculate Predicted Values by iteratively adding Predictors
#3) Calculate R2 & AdjR2 for Training and Test Sets
#4) Keep track of:
  # Number of Parameters
  # Which Parameters
  # R2/AdjR2 (Test & Train)
  # 

## Retreive Coefficients for Covariates Included in Model
TAB.TR.cov <- MG.TR[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
MOD.cov <- lm( Pheno ~ . , data=TAB.TR.cov )
BETA.cov <- coef( MOD.cov )
P.cov <- summary(MOD.cov)$coefficients[,"Pr(>|t|)"]

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt )
P.all <- c( P.cov[ORD.cov], P.gt[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt)[ORD.gt] )
BETA.all <- c( BETA.cov[ORD.cov], BETA.gt[ORD.gt] ) ; names(BETA.all) <- c( names(BETA.cov)[ORD.cov], names(BETA.gt)[ORD.gt] )
COEF.all <- data.frame( P.all, BETA.all ) ; rownames(COEF.all) <- names(BETA.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]

## Loop through Predictors and add each one iteratively to model
 # Then assess fit of model
R2 <- array( , dim=c(nrow(COEF.all),4) )
rownames(R2) <- rownames(COEF.all)
colnames(R2) <- c("TR","TS","Adj.TR","Adj.TS")
# for ( r in 1:nrow(COEF.all) ) {
for ( r in 1:20 ) {
	This_Pred <- rownames(COEF.all)[r]
	Which_Pred <- rownames(COEF.all)[r]
	Which_Pred <- rownames(COEF.all)[1:r]
	## Build Temporary Data Frame with These Covariates
	TEMP.TR <- data.matrix( MG.TR[, Which_Pred ] )
	TEMP.TS <- data.matrix( MG.TS[, Which_Pred ] )
	## Predict Values from Test Set based on Model with these Predictors
	PRED.TS <- COEF.int$BETA.all + TEMP.TS %*% COEF.all[Which_Pred,"BETA.all"]
	PRED.TR <- COEF.int$BETA.all + TEMP.TR %*% COEF.all[Which_Pred,"BETA.all"]
	## Calculate R2 for covariate model
	 # Calc Sum Squares
	DIFF.TS.mod <- PRED.TS - MG.TS$Pheno
	DIFF.TR.mod <- PRED.TR - MG.TR$Pheno
	SS.TS.mod <- sum(DIFF.TS.mod^2) - sum(DIFF.TS.mod)^2 / length(DIFF.TS.mod)
	SS.TR.mod <- sum(DIFF.TR.mod^2) - sum(DIFF.TR.mod)^2 / length(DIFF.TR.mod)
	SS.TS.val <- sum(MG.TS$Pheno^2) - sum(MG.TS$Pheno)^2 / nrow(MG.TS)
	SS.TR.val <- sum(MG.TR$Pheno^2) - sum(MG.TR$Pheno)^2 / nrow(MG.TR)
	 # R2
	R2.TS.all <- 1 - SS.TS.mod / SS.TS.val
	R2.TR.all <- 1 - SS.TR.mod / SS.TR.val
	 # Adjusted R2
	R2.adj.TS.all <- 1 - ( SS.TS.mod/(nrow(MG.TS)-r) / (SS.TS.val/(nrow(MG.TS)-1)) )
	R2.adj.TR.all <- 1 - ( SS.TR.mod/(nrow(MG.TR)-r) / (SS.TR.val/(nrow(MG.TR)-1)) )
	## Compile all the results for Plotting
	R2[r,"TR"] <- R2.TR.all
	R2[r,"TS"] <- R2.TS.all
	R2[r,"Adj.TR"] <- R2.adj.TR.all
	R2[r,"Adj.TS"] <- R2.adj.TS.all
	## Output Status
	if (r%%10==0) { 
		# print( tail(Which_Pred) )
		print(paste( "Finished",r,"of",nrow(COEF.all) ))
	}
}
head(R2,20)

## Predict values of phenotype for Test Group
 # With all Covariates
TAB.TS.all.cov <- data.matrix( MG.TS[, Cov_Cols ] ) 
PRED.TS.all <- BETA.cov["(Intercept)"] + TAB.TS.all.cov %*% BETA.cov[2:length(BETA.cov)]

## Calculate R2 for covariate model
 # Calc Sum Squares
DIFF.TS.all <- PRED.TS.all - MG.TS$Pheno
SS.TS.all <- sum(DIFF.TS.all^2) - sum(DIFF.TS.all)^2 / length(DIFF.TS.all)
SS.TS.val <- sum(MG.TS$Pheno^2) - sum(MG.TS$Pheno)^2 / nrow(MG.TS)
 # R2
R2.all <- 1 - SS.TS.all / SS.TS.val
R2.sig <- 1 - SS.TS.sig / SS.TS.val
 # Adjusted R2
R2.adj.all <- 1 - ( SS.TS.all/(nrow(MG.TS)-length(BETA.cov)) / (SS.TS.val/(nrow(MG.TS)-1)) )
R2.adj.sig <- 1 - ( SS.TS.sig/(nrow(MG.TS)-length(BETA.cov)) / (SS.TS.val/(nrow(MG.TS)-1)) )



###############################################################
## MODEL w/ COVARIATES ########################################
###############################################################

## Retreive Coefficients for Covariates Included in Model
TAB.TR.cov <- MG.TR[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
MOD.cov <- lm( Pheno ~ . , data=TAB.TR.cov )
BETA.cov <- coef( MOD.cov )

## Predict values of phenotype for Test Group
 # With all Covariates
TAB.TS.all.cov <- data.matrix( MG.TS[, Cov_Cols ] ) 
PRED.TS.all <- BETA.cov["(Intercept)"] + TAB.TS.all.cov %*% BETA.cov[2:length(BETA.cov)]
 # With only Significant Covariates
SIG.cov <- rownames(anova(MOD.cov))[ which( anova(MOD.cov)[,"Pr(>F)"]<.05 ) ]
TAB.TS.sig.cov <- data.matrix( MG.TS[, SIG.cov ] ) 
PRED.TS.sig <- BETA.cov["(Intercept)"] + TAB.TS.sig.cov %*% BETA.cov[SIG.cov]
cor( data.frame(PRED.TS.all, PRED.TS.sig, MG.TS$Pheno) )
pairs( data.frame(PRED.TS.all, PRED.TS.sig, MG.TS$Pheno) )
pairs( data.frame(PRED.TS.all, PRED.TS.sig, MG.TS$Pheno) )

## Calculate R2 for covariate model
 # Calc Sum Squares
DIFF.TS.all <- PRED.TS.all - MG.TS$Pheno
DIFF.TS.sig <- PRED.TS.sig - MG.TS$Pheno
SS.TS.all <- sum(DIFF.TS.all^2) - sum(DIFF.TS.all)^2 / length(DIFF.TS.all)
SS.TS.sig <- sum(DIFF.TS.sig^2) - sum(DIFF.TS.sig)^2 / length(DIFF.TS.sig)
SS.TS.val <- sum(MG.TS$Pheno^2) - sum(MG.TS$Pheno)^2 / nrow(MG.TS)
 # R2
R2.all <- 1 - SS.TS.all / SS.TS.val
R2.sig <- 1 - SS.TS.sig / SS.TS.val
 # Adjusted R2
R2.adj.all <- 1 - ( SS.TS.all/(nrow(MG.TS)-length(BETA.cov)) / (SS.TS.val/(nrow(MG.TS)-1)) )
R2.adj.sig <- 1 - ( SS.TS.sig/(nrow(MG.TS)-length(BETA.cov)) / (SS.TS.val/(nrow(MG.TS)-1)) )

##

## Try a bunch of different models and go with the one w/ the best Adjusted R2
 # This will get out of hand if I use all variants/covariates
 # 

###############################################################
## END OF DOC #################################################
###############################################################
