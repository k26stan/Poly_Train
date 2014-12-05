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
GT.which.1 <- match( KP.P.id, GT.nms ) # which( GT.nms %in% KP.P.id )
# GT.which.2 <- match( GT.nms, KP.P.id ) # which( GT.nms %in% KP.P.id )
# GT.which.nms <- GT.nms[ GT.which.1 ]
GT.cand.1 <- data.matrix( GT[,GT.which.1] )
rownames(GT.cand.1) <- as.character( GT$FID )
CP.4 <- CP.4[ which(KP.P.id %in% GT.nms[GT.which.1] ), ]
KP.P.id <- gsub( ":",".", CP.4$SNP, fixed=T )
KP.P.id <- gsub( "X","",KP.P.id )

## Get Beta values for each Variant that made the cut
BETA.gt.1 <- CP.4$BETA ; names(BETA.gt.1) <- KP.P.id
P.gt.1 <- CP.4$P_Assoc ; names(P.gt.1) <- KP.P.id

## Remove Variants w/ Correlated Genotypes
GT.corr.1 <- cor( GT.cand.1 )
THRSH.corr <- 0.9
RM.corr <- c()
in_loop <- 0
for ( v in 1:length(P.gt.1) ) {
	var <- names(P.gt.1)[v]
	row <- which( rownames(GT.corr)==var )
	if ( !(var %in% RM.corr) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,":",v,"-",var,"######" ))
		print( range( (row+1):ncol(GT.corr) ) )
		RM.temp.1 <- colnames(GT.corr)[row+which(GT.corr[var,(row+1):ncol(GT.corr)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr)[row+which(GT.corr[var,(row+1):ncol(GT.corr)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.1, RM.temp.2 ) )
		print(paste("",length(P.gt.1)-length(RM.corr),"Remaining" ))
	}
}
BETA.gt <- BETA.gt.1[ -which( names(BETA.gt.1) %in% RM.corr ) ] #; names(BETA.gt) <- KP.P.id[ -which( names(BETA.gt.1) %in% RM.corr ) ]
P.gt <- P.gt.1[ -which( names(P.gt.1) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.1) %in% RM.corr ) ]
GT.cand <- GT.cand.1[ ,-which( colnames(GT.cand.1) %in% RM.corr ) ]
# GT.which <- GT.which.1[ -which( names(GT.which.1) %in% RM.corr ) ]
GT.corr <- cor( GT.cand )
par(mfrow=c(1,2))
hist( GT.corr.1, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )
hist( GT.corr, breaks=seq(-1,1,.01),col="blue" ) ; abline( v=seq(-1,1,.2),lty=2,col="red" )

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
#1) Try Covariate Predictors Individually
#2) Take best and include in model
#3) Try rest of Covariates after including best
#4) Repeat 2-3 until all Covariates are Included
#5) Try Variant Predictors Individually (including covariates)
#6) Take best Variant and include in model
#7) Try rest of Variants w/ all previous predictors included
#8) Repeat until the best model is formed

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
################# CAN NOT USE INDIVIDUAL BETA VALUES WITH MULTIPLE PREDICTORS #################
TO_TEST.cov <- 1:(length(BETA.cov)-1)
TO_TEST.gt <- length(TO_TEST.cov) + 1:length(BETA.gt)
TO_TEST <- TO_TEST.cov
TO_INCL <- c()
TO_OMIT <- c()
MOD_SIZE <- 150 # nrow(COEF.all) # 20
N_PREDS <- nrow(COEF.all)
R2.test <- list()
R2.test$TS <- R2.test$TR <- R2.test$TS.adj <- R2.test$TR.adj <- array( , c(N_PREDS,MOD_SIZE) )
rownames(R2.test$TS) <- rownames(R2.test$TR) <- rownames(R2.test$TS.adj) <- rownames(R2.test$TR.adj) <- rownames(COEF.all)[1:N_PREDS]
R2 <- array( , dim=c(MOD_SIZE,4) )
colnames(R2) <- c("TR","TS","Adj.TR","Adj.TS")
# for ( r in 1:nrow(COEF.all) ) {
start_time <- proc.time()
for ( i in 1:MOD_SIZE ) {
	## Switch to Genotypes
	if ( i == length(TO_TEST.cov)+1 ) { TO_TEST <- TO_TEST.gt[1:(N_PREDS-length(TO_TEST.cov))] }
	## Loop through all Predictors and Find Best on to ADD to model
	for ( r in TO_TEST ) {
		This_Pred <- rownames(COEF.all)[r]
		Which_Pred <- c( TO_INCL, This_Pred )
		## Build Temporary Data Frame with These Covariates
		TEMP.TR <- data.matrix( MG.TR[, Which_Pred ] )
		TEMP.TS <- data.matrix( MG.TS[, Which_Pred ] )
		## Predict Values from Test Set based on Model with these Predictors
		PRED.TS <- COEF.int$BETA.all + TEMP.TS %*% COEF.all[Which_Pred,"BETA.all"]
		PRED.TR <- COEF.int$BETA.all + TEMP.TR %*% COEF.all[Which_Pred,"BETA.all"]
		## Fit Model Using lm()
		MOD.TR <- lm( MG.TR[,"Pheno"] ~ TEMP.TR )
		COEF.mod <- coef(MOD.TR)
		PRED.TR <- COEF.mod[1] + TEMP.TR %*% matrix(COEF.mod[-1])
		PRED.TS <- COEF.mod[1] + TEMP.TS %*% matrix(COEF.mod[-1])
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
		R2.adj.TS.all <- 1 - ( SS.TS.mod/(nrow(MG.TS)-length(Which_Pred)) / (SS.TS.val/(nrow(MG.TS)-1)) )
		R2.adj.TR.all <- 1 - ( SS.TR.mod/(nrow(MG.TR)-length(Which_Pred)) / (SS.TR.val/(nrow(MG.TR)-1)) )
		 # From Mod
		# R2.TR.all <- summary(MOD.TR)$r.squared
		# R2.adj.TR.all <- summary(MOD.TR)$adj.r.squared
		 # Compile Results to Compare for this Iteration
		R2.test$TS[r,i] <- R2.TS.all
		R2.test$TR[r,i] <- R2.TR.all
		R2.test$TS.adj[r,i] <- R2.adj.TS.all
		R2.test$TR.adj[r,i] <- R2.adj.TR.all
	}
	## Decide Which Model is the Best
	WHICH_TO_INC <- which.max( R2.test$TR[,i] )
	# WHICH_TO_INC <- which.max( R2.test$TR.adj[,i] )
	Pred_To_Inc <- rownames(COEF.all)[WHICH_TO_INC]
	print(paste( "Include",WHICH_TO_INC,"-",Pred_To_Inc ))
	## Compile all the results for Plotting
	R2[i,"TR"] <- R2.test$TR[WHICH_TO_INC,i]
	R2[i,"TS"] <- R2.test$TS[WHICH_TO_INC,i]
	R2[i,"Adj.TR"] <- R2.test$TR.adj[WHICH_TO_INC,i]
	R2[i,"Adj.TS"] <- R2.test$TS.adj[WHICH_TO_INC,i]
	## Include, but don't Test, that Predictor from now on
	TO_INCL <- c( TO_INCL, Pred_To_Inc )
	TO_TEST <- setdiff( TO_TEST, WHICH_TO_INC )
	## Output Status
	if (i%%10==0) { 
		# print( tail(Which_Pred) )
		print(paste( "Finished",i,"of",MOD_SIZE,"-",round( (proc.time()-start_time)[3], 1) ))
	}
}
rownames(R2) <- TO_INCL
colnames(R2.test$TS) <- colnames(R2.test$TR) <- colnames(R2.test$TS.adj) <- colnames(R2.test$TR.adj) <- TO_INCL
head(R2,20)

## Plot R2 Values for Test/Training Sets
XLIM <- c(1,MOD_SIZE)
YLIM <- c(-1,1) # range(R2)
COLS <- c( rep(c("chartreuse2","dodgerblue2"),2) ) # , rep("dodgerblue2",2) )
LTYS <- c( 1,1,2,2 )
plot( 0,0, type="n", xlim=XLIM, ylim=YLIM, xaxt="n", xlab="Predictors Included (Cummulative)", ylab="R2 & Adj R2", main="Variance Explained by Model" )
axis( 1, at=1:MOD_SIZE, label=TO_INCL, las=2 )
abline( h=seq( -1,1,.1), lty=c(rep(2,10),1,rep(2,10)), col=c(rep("firebrick3",10),"black",rep("grey50",10)) ) 
abline( v=seq( 0,MOD_SIZE,5), lty=2, col="grey50" ) 
for ( i in 1:4 ) {
	points( 1:nrow(R2), R2[,i], col=COLS[i], lty=LTYS[i], type="l", lwd=3 )
}
legend( quantile(XLIM,.8), quantile(YLIM,.9), legend=colnames(R2), lty=LTYS, col=COLS, lwd=3 )
## Plot This Shiz
COLS <- c( "gold1","chocolate2","firebrick3","black","slateblue3","steelblue2","springgreen1" )
HT_COLS <- colorRampPalette(COLS)(50)
## Correlation b/n Predictors
COR.INCL <- cor( MG.TR[, TO_INCL], use="pairwise.complete.obs", method="spearman" ) 
BRKS <- seq(-1,1,length.out=51)
heatmap.2( COR.INCL, col=HT_COLS, breaks=BRKS, scale="none", Colv=F, Rowv=T, dendrogram="row", trace="none" )
## Heatmap of Predictor Tests (R2.test$TR)
BRKS <- c( seq(min(R2.test$TR,na.rm=T),0,length.out=26)[1:25], 0, seq(0,5,length.out=26)[2:26] )
heatmap.2( R2.test$TR, col=HT_COLS, breaks=BRKS, scale="column", Colv=F, Rowv=F, dendrogram="none", trace="none" )
heatmap.2( t(R2.test$TR), col=HT_COLS, scale="row", Colv=F, Rowv=F, dendrogram="none", trace="none" )

###############################################################
## END OF DOC #################################################
###############################################################








































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