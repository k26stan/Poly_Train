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

library(gplots)

## Pull in Command Line Arguments
LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( paste( "/Users/kstandis/Downloads/20141124_Poly_Train/ACR20/", c("ACR20_100wk", "Cov_w_PCs.txt", "CND_ACR20_100wk_AGE_SEX_PC1_PC2.compile.short", "CND_ACR20_100wk_AGE_SEX_PC1_PC2_012.raw"), sep="" ), "AGE,SEX,PC1,PC2" )
# LINE <- c( paste( "/projects/janssen/Poly_Train/20150105_1_ACR50_100wk_AGE_SEX/", c("ACR50_100wk", "Cov_w_PCs.txt", "CND_ACR50_100wk_AGE_SEX.compile.short", "CND_ACR50_100wk_AGE_SEX_012.raw"), sep="" ), "AGE,SEX" )
# LINE <- c( paste( "/projects/janssen/Poly_Train/20150105_2_ACR50_100wk_AGE_SEX/", c("ACR50_100wk", "Cov_w_PCs.txt", "CND_ACR50_100wk_AGE_SEX.compile.short", "CND_ACR50_100wk_AGE_SEX_012.raw"), sep="" ), "AGE,SEX" )

## Parse Command Line Arguments
PathToPheno <- LINE[1]
PathToCovs <- LINE[2]
PathToComp <- LINE[3]
PathToGeno <- LINE[4]
Covs_Command <- LINE[5]

## Set Path to Save
Split_1 <- strsplit(PathToGeno,"/")[[1]]
PathToSave <- paste( c(Split_1[1:(length(Split_1)-1)],""), collapse="/" )
File_Name <- gsub( "CND_","", gsub( "_012.raw","", Split_1[length(Split_1)] ) )

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
if (length(RM.DUP)>0) { CP.1 <- CP[ -RM.DUP, ] }else{ CP.1 <- CP }

## Remove HW Violations
THRSH.HWE <- 1e-8
RM.HWE <- which( CP.1$P_HW < THRSH.HWE )
if (length(RM.HWE)>0) { CP.2 <- CP.1[ -RM.HWE, ] }else{ CP.2 <- CP.1 }

## Pull out Variants w/ MAF > p%
THRSH.MAF <- 0.05
RM.MAF <- which( CP.2$REF_ALL < THRSH.MAF | CP.2$ALT_ALL < THRSH.MAF )
if (length(RM.MAF)>0) { CP.3 <- CP.2[ -RM.MAF, ] }else{ CP.3 <- CP.2 }

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
OR.gt.a <- CP.4$OR ; names(OR.gt.a) <- KP.P.id
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
	if ( !(var %in% RM.corr) & row<ncol(GT.corr.a) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,":",v,"-",var,"######" ))
		print( range( (row+1):ncol(GT.corr.a) ) )
		RM.temp.a <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.a, RM.temp.2 ) )
		print(paste("",length(P.gt.a)-length(RM.corr),"Remaining" ))
	}
}
OR.gt.9 <- OR.gt.a[ -which( names(OR.gt.a) %in% RM.corr ) ] #; names(OR.gt) <- KP.P.id[ -which( names(OR.gt.a) %in% RM.corr ) ]
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
	if ( !(var %in% RM.corr) & row<ncol(GT.corr.a) ) {
		in_loop <- in_loop + 1
		print(paste("######",in_loop,":",v,"-",var,"######" ))
		print( range( (row+1):ncol(GT.corr.a) ) )
		RM.temp.a <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] > THRSH.corr)]
		RM.temp.2 <- colnames(GT.corr.a)[row+which(GT.corr.a[var,(row+1):ncol(GT.corr.a)] < -THRSH.corr)]
		RM.corr <- Reduce( union, list( RM.corr, RM.temp.a, RM.temp.2 ) )
		print(paste("",length(P.gt.a)-length(RM.corr),"Remaining" ))
	}
}
OR.gt.1 <- OR.gt.a[ -which( names(OR.gt.a) %in% RM.corr ) ] #; names(OR.gt) <- KP.P.id[ -which( names(OR.gt.a) %in% RM.corr ) ]
P.gt.1 <- P.gt.a[ -which( names(P.gt.a) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.a) %in% RM.corr ) ]
GT.cand.1 <- GT.cand.a[ ,-which( colnames(GT.cand.a) %in% RM.corr ) ]
GT.corr.1 <- cor( GT.cand.1 )
## Compare
jpeg( paste(PathToSave,"MOD_SNP_Corr",File_Name,".jpeg",sep=""), height=1000, width=2000, pointsize=30 )
par(mfrow=c(1,3))
hist( GT.corr.a, breaks=seq(-1,1,.01),col="dodgerblue2", main="Variant Correlation: All", xlab="R" ) ; abline( v=seq(-1,1,.2),lty=2,col="firebrick3" )
hist( GT.corr.9, breaks=seq(-1,1,.01),col="dodgerblue2", main="Variant Correlation: R2 < .9", xlab="R" ) ; abline( v=seq(-1,1,.2),lty=2,col="firebrick3" )
hist( GT.corr.1, breaks=seq(-1,1,.01),col="dodgerblue2", main="Variant Correlation: R2 < .1", xlab="R" ) ; abline( v=seq(-1,1,.2),lty=2,col="firebrick3" )
dev.off()

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
# heatmap.2( TEMP.TR.9, scale="none", trace="none")
MG.TR.1 <- MG.TR.2.1[,3:ncol(MG.TR.2.1)] ; rownames(MG.TR.1) <- as.character( MG.TR.2.1[,1] )
MG.TS.1 <- MG.TS.2.1[,3:ncol(MG.TS.2.1)] ; rownames(MG.TS.1) <- as.character( MG.TS.2.1[,1] )
TEMP.TR.1 <- cor( MG.TR.1[87:ncol(MG.TR.1)] )
TEMP.TS.1 <- cor( MG.TS.1[87:ncol(MG.TS.1)] )
# heatmap.2( TEMP.TR.1, scale="none", trace="none")

###############################################################
## MODEL w/ REGRESSION FIT of CORRELATED GENOTYPES ############
###############################################################
## Using <.9 correlated genotypes:
 # throw each predictor into a multiple regression (including best GT each time)
 # and use OR values calculated each time on TEST set

## Retreive Coefficients for Covariates Included in Model
TAB.TR.9.cov <- MG.TR.9[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
# MOD.cov.lm <- lm( Pheno ~ . , data=TAB.TR.9.cov )
MOD.cov <- glm( factor(Pheno) ~ . , data=TAB.TR.9.cov, family=binomial(link=logit) )
OR.cov <- coef( MOD.cov )
P.cov <- summary(MOD.cov)$coefficients[,"Pr(>|z|)"]

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt.9 )
P.all <- c( P.cov[ORD.cov], P.gt.9[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt.9)[ORD.gt] )
OR.all <- c( OR.cov[ORD.cov], OR.gt.9[ORD.gt] ) ; names(OR.all) <- c( names(OR.cov)[ORD.cov], names(OR.gt.9)[ORD.gt] )
COEF.all <- data.frame( P.all, OR.all ) ; rownames(COEF.all) <- names(OR.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]

## Loop through Predictors and add each one iteratively to model
 # Then assess fit of model
TO_TEST.cov <- 1:(length(OR.cov)-1)
TO_TEST.gt <- length(TO_TEST.cov) + 1:length(OR.gt.9)
TO_TEST <- TO_TEST.cov
TO_INCL <- c()
TO_OMIT <- c()
MOD_SIZE <- 30 ; MOD_SIZE <- min( MOD_SIZE, nrow(COEF.all) )
N_PREDS <- nrow(COEF.all)
AUC.test <- list()
AUC.test$TS <- AUC.test$TR <- array( , c(N_PREDS,MOD_SIZE) )
rownames(AUC.test$TS) <- rownames(AUC.test$TR) <- rownames(COEF.all)[1:N_PREDS]
AUC <- array( , dim=c(MOD_SIZE,2) )
colnames(AUC) <- c("TR","TS")
# for ( r in 1:nrow(COEF.all) ) {
start_time <- proc.time()
for ( i in 1:MOD_SIZE ) {
	## Switch to Genotypes
	if ( i == length(TO_TEST.cov)+1 ) { TO_TEST <- TO_TEST.gt[1:(N_PREDS-length(TO_TEST.cov))] }
	## Loop through all Predictors and Find Best on to ADD to model
	for ( r in TO_TEST ) {
		This_Pred <- rownames(COEF.all)[r]
		Which_Pred <- c( TO_INCL, This_Pred )
		## Specify Phenotype Values
		PHENO.TR.9 <- factor(MG.TR.9[,"Pheno"])
		PHENO.TS.9 <- factor(MG.TS.9[,"Pheno"])
		## Build Temporary Data Frame with These Covariates
		TEMP.TR.9 <- data.matrix( MG.TR.9[, Which_Pred ] ) ; colnames(TEMP.TR.9) <- Which_Pred
		TEMP.TS.9 <- data.matrix( MG.TS.9[, Which_Pred ] ) ; colnames(TEMP.TS.9) <- paste("TEMP.TR.9",Which_Pred,sep="")
		## Fit Model Using lm()
		MOD.TR.9 <- glm( PHENO.TR.9 ~ TEMP.TR.9, family=binomial(link=logit) )
		COEF.mod <- coef(MOD.TR.9)
		PRED.TR.9 <- COEF.mod[1] + TEMP.TR.9 %*% matrix(COEF.mod[-1])
		PRED.TS.9 <- COEF.mod[1] + TEMP.TS.9 %*% matrix(COEF.mod[-1])
		## Calculate Probability from Link Function
		 # p <- exp( b0 + b1X1 + ... + bnXn ) / ( 1 + exp( b0 + b1X1 + ... + bnXn ) )
		PROB.TR.9 <- exp(PRED.TR.9) / ( 1 + exp(PRED.TR.9) )
		PROB.TS.9 <- exp(PRED.TS.9) / ( 1 + exp(PRED.TS.9) )
		## Loop through different thresholds for calling Case/Control
		THRESHOLDS <- seq( 0,1,.002 )
		N_THRSH <- length(THRESHOLDS)
		TPR.TR.9 <- TPR.TS.9 <- FPR.TR.9 <- FPR.TS.9 <- numeric(N_THRSH)
		for ( t in 1:N_THRSH ) {
			thresh <- THRESHOLDS[t]
			## Predict Cases/Controls
			BIN.TR.9 <- PROB.TR.9 > thresh
			BIN.TS.9 <- PROB.TS.9 > thresh
			## Calculate False Positive and True Postive Rates
			TPR.TR.9[t] <- length(which( PHENO.TR.9==2 & BIN.TR.9==T )) / length(which( PHENO.TR.9==2 ))
			TPR.TS.9[t] <- length(which( PHENO.TS.9==2 & BIN.TS.9==T )) / length(which( PHENO.TS.9==2 ))
			FPR.TR.9[t] <- length(which( PHENO.TR.9==1 & BIN.TR.9==T )) / length(which( PHENO.TR.9==1 ))
			FPR.TS.9[t] <- length(which( PHENO.TS.9==1 & BIN.TS.9==T )) / length(which( PHENO.TS.9==1 ))
		}
		# plot( 0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="FPR",ylab="TPR", main="ROC" )
		# abline( 0,1, lty=1, col="black", lwd=2 )
		# abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
		# abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
		# points( FPR.TR.9, TPR.TR.9, type="o", pch=20, col="chartreuse2", lwd=2 )
		# points( FPR.TS.9, TPR.TS.9, type="o", pch=20, col="dodgerblue2", lwd=2 )
		## Calculate Area Under Curve
		AUC.test$TR[r,i] <- sum( -( FPR.TR.9[2:N_THRSH]-FPR.TR.9[2:N_THRSH-1] ) * ( TPR.TR.9[2:N_THRSH]+TPR.TR.9[2:N_THRSH-1] ) /2 )
		AUC.test$TS[r,i] <- sum( -( FPR.TS.9[2:N_THRSH]-FPR.TS.9[2:N_THRSH-1] ) * ( TPR.TS.9[2:N_THRSH]+TPR.TS.9[2:N_THRSH-1] ) /2 )
	}
	## Decide Which Model is the Best
	WHICH_TO_INC <- which.max( AUC.test$TR[,i] )
	Pred_To_Inc <- rownames(COEF.all)[WHICH_TO_INC]
	print(paste( "Include",WHICH_TO_INC,"-",Pred_To_Inc ))
	## Compile all the results for Plotting
	AUC[i,"TR"] <- AUC.test$TR[WHICH_TO_INC,i]
	AUC[i,"TS"] <- AUC.test$TS[WHICH_TO_INC,i]
	## Include, but don't Test, that Predictor from now on
	TO_INCL <- c( TO_INCL, Pred_To_Inc )
	TO_TEST <- setdiff( TO_TEST, WHICH_TO_INC )
	## Output Status
	if (i%%5==0) { 
		# print( tail(Which_Pred) )
		print(paste( "Finished",i,"of",MOD_SIZE,"-",round( (proc.time()-start_time)[3], 1) ))
	}
}
rownames(AUC) <- TO_INCL
colnames(AUC.test$TS) <- colnames(AUC.test$TR) <- TO_INCL
head(AUC,20)
AUC.9 <- AUC
AUC.9.test <- AUC.test

## Plot AUC.9 Values for Test/Training Sets
XLIM <- c(1,nrow(AUC.9))
YLIM <- c(0,1) # range(AUC.9)
COLS <- c( rep(c("chartreuse2","dodgerblue2"),2) ) # , rep("dodgerblue2",2) )
LTYS <- c( 1,1,2,2 )
jpeg( paste(PathToSave,"MOD_Reg_",File_Name,".jpeg",sep=""), height=1000, width=2000, pointsize=30 )
plot( 0,0, type="n", xlim=XLIM, ylim=YLIM, xaxt="n", xlab="Predictors Included (Cummulative)", ylab="AUC", main="AUC for Training/Test Sets (Regression Model)" )
axis( 1, at=1:nrow(AUC.9), label=TO_INCL, las=2 )
abline( h=seq( 0,1,.1), lty=c(rep(2,5),1,rep(2,5)), col=c(rep("firebrick3",5),"black",rep("grey50",5)) ) 
abline( v=seq( 0,nrow(AUC.9),5), lty=2, col="grey50" ) 
for ( i in 1:2 ) {
	points( 1:nrow(AUC.9), AUC.9[,i], col=COLS[i], lty=LTYS[i], type="l", lwd=3 )
}
legend( quantile(XLIM,.8), quantile(YLIM,.9), legend=colnames(AUC.9), lty=LTYS, col=COLS, lwd=3 )
dev.off()

###############################################################
## MODEL w/ +/- SUM of GT COUNTS ##############################
###############################################################
## Using <.1 correlated genotypes:
 # calculate sum (incl. +/- for OR) of genotypes person to yield a score
 # fit that single score to response data
 # and use OR value for score calculated each time on TEST set

## Retreive Coefficients for Covariates Included in Model
TAB.TR.1.cov <- MG.TR.1[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
# MOD.cov.lm <- lm( Pheno ~ . , data=TAB.TR.1.cov )
MOD.cov <- glm( factor(Pheno) ~ . , data=TAB.TR.1.cov, family=binomial(link=logit) )
OR.cov <- coef( MOD.cov )
P.cov <- summary(MOD.cov)$coefficients[,"Pr(>|z|)"]

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt.1 )
P.all <- c( P.cov[ORD.cov], P.gt.1[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt.1)[ORD.gt] )
OR.all <- c( OR.cov[ORD.cov], OR.gt.1[ORD.gt] ) ; names(OR.all) <- c( names(OR.cov)[ORD.cov], names(OR.gt.1)[ORD.gt] )
COEF.all <- data.frame( P.all, OR.all ) ; rownames(COEF.all) <- names(OR.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all$OR.all[(length(OR.cov)):nrow(COEF.all)] <- c(-1,1)[factor(COEF.all$OR.all[(length(OR.cov)):nrow(COEF.all)]>0)]

## Loop through Predictors and add each one iteratively to model
 # Then assess fit of model
TO_TEST.cov <- 1:(length(OR.cov)-1)
TO_TEST.gt <- length(TO_TEST.cov) + 1:length(OR.gt.1)
TO_TEST <- TO_TEST.cov
TO_INCL <- TO_OMIT <- c()
GT_Flag <- 0
MOD_SIZE <- 40 ; MOD_SIZE <- min( MOD_SIZE, nrow(COEF.all) ) # 20
N_PREDS <- nrow(COEF.all)
AUC.test <- list()
AUC.test$TS <- AUC.test$TR <- array( , c(N_PREDS,MOD_SIZE) )
rownames(AUC.test$TS) <- rownames(AUC.test$TR) <- rownames(COEF.all)[1:N_PREDS]
SCOR.test <- list()
SCOR.test$P <- SCOR.test$OR <- array( , c(N_PREDS,MOD_SIZE) )
rownames(SCOR.test$P) <- rownames(SCOR.test$OR) <- rownames(COEF.all)[1:N_PREDS]
AUC <- array( , dim=c(MOD_SIZE,2) )
colnames(AUC) <- c("TR","TS")
start_time <- proc.time()
for ( i in 1:MOD_SIZE ) {
# for ( i in 1:8 ) {
	## After Covariates, Switch to Genotypes
	if ( i == length(TO_TEST.cov)+1 ) {
		TO_TEST <- TO_TEST.gt[1:(N_PREDS-length(TO_TEST.cov))]
		GT_Flag <- 1
	}
	## Loop through all Predictors and Find Best one to ADD to model
	for ( r in TO_TEST ) {
		This_Pred <- rownames(COEF.all)[r]
		Which_Pred <- c( TO_INCL, This_Pred )
		Which_GT <- setdiff( Which_Pred, names(OR.cov) )
		Which_Cov <- setdiff( Which_Pred, Which_GT )
		## Specify Phenotype Values
		PHENO.TR.1 <- factor(MG.TR.9[,"Pheno"])
		PHENO.TS.1 <- factor(MG.TS.9[,"Pheno"])
		## Create Binary Matrix for Genotypes
		BIN.TR.1 <- data.matrix( MG.TR.1[, Which_GT ] ) ; BIN.TR.1[which(BIN.TR.1==2,arr.ind=T)] <- 1
		BIN.TS.1 <- data.matrix( MG.TS.1[, Which_GT ] ) ; BIN.TS.1[which(BIN.TS.1==2,arr.ind=T)] <- 1
		## Calculate Score for each Person
		SCOR.TR.1 <- BIN.TR.1 %*% COEF.all[Which_GT,"OR.all"]
		SCOR.TS.1 <- BIN.TS.1 %*% COEF.all[Which_GT,"OR.all"]
		## Build Temporary Data Frame with These Covariates
		if (GT_Flag==1) {
			TEMP.TR.1 <- data.matrix(data.frame( MG.TR.1[, c(Which_Cov) ], SCOR.TR.1 ))
			TEMP.TS.1 <- data.matrix(data.frame( MG.TS.1[, c(Which_Cov) ], SCOR.TS.1 ))
			}else{
			TEMP.TR.1 <- data.matrix( MG.TR.1[, c(Which_Cov) ] )
			TEMP.TS.1 <- data.matrix( MG.TS.1[, c(Which_Cov) ] )
		}
		## Fit Model Using lm()
		MOD.TR.1 <- glm( PHENO.TR.1 ~ TEMP.TR.1, family=binomial(link=logit) )
		COEF.mod <- coef(MOD.TR.1)
		PRED.TR.1 <- COEF.mod[1] + TEMP.TR.1 %*% matrix(COEF.mod[-1])
		PRED.TS.1 <- COEF.mod[1] + TEMP.TS.1 %*% matrix(COEF.mod[-1])
		## Calculate Probability from Link Function
		 # p <- exp( b0 + b1X1 + ... + bnXn ) / ( 1 + exp( b0 + b1X1 + ... + bnXn ) )
		PROB.TR.1 <- exp(PRED.TR.1) / ( 1 + exp(PRED.TR.1) )
		PROB.TS.1 <- exp(PRED.TS.1) / ( 1 + exp(PRED.TS.1) )
		## Loop through different thresholds for calling Case/Control
		THRESHOLDS <- seq( 0,1,.002 )
		N_THRSH <- length(THRESHOLDS)
		TPR.TR.1 <- TPR.TS.1 <- FPR.TR.1 <- FPR.TS.1 <- numeric(N_THRSH)
		for ( t in 1:N_THRSH ) {
			thresh <- THRESHOLDS[t]
			## Predict Cases/Controls
			BIN.TR.1 <- PROB.TR.1 > thresh
			BIN.TS.1 <- PROB.TS.1 > thresh
			## Calculate False Positive and True Postive Rates
			TPR.TR.1[t] <- length(which( PHENO.TR.1==2 & BIN.TR.1==T )) / length(which( PHENO.TR.1==2 ))
			TPR.TS.1[t] <- length(which( PHENO.TS.1==2 & BIN.TS.1==T )) / length(which( PHENO.TS.1==2 ))
			FPR.TR.1[t] <- length(which( PHENO.TR.1==1 & BIN.TR.1==T )) / length(which( PHENO.TR.1==1 ))
			FPR.TS.1[t] <- length(which( PHENO.TS.1==1 & BIN.TS.1==T )) / length(which( PHENO.TS.1==1 ))
		}
		# plot( 0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="FPR",ylab="TPR", main="ROC" )
		# abline( 0,1, lty=1, col="black", lwd=2 )
		# abline( h=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
		# abline( v=seq(0,1,.1), lty=2, col="grey50", lwd=1 )
		# points( FPR.TR.1, TPR.TR.1, type="o", pch=20, col="chartreuse2", lwd=2 )
		# points( FPR.TS.1, TPR.TS.1, type="o", pch=20, col="dodgerblue2", lwd=2 )
		## Calculate Area Under Curve
		AUC.test$TR[r,i] <- sum( -( FPR.TR.1[2:N_THRSH]-FPR.TR.1[2:N_THRSH-1] ) * ( TPR.TR.1[2:N_THRSH]+TPR.TR.1[2:N_THRSH-1] ) /2 )
		AUC.test$TS[r,i] <- sum( -( FPR.TS.1[2:N_THRSH]-FPR.TS.1[2:N_THRSH-1] ) * ( TPR.TS.1[2:N_THRSH]+TPR.TS.1[2:N_THRSH-1] ) /2 )
	}
	## Decide Which Model is the Best
	WHICH_TO_INC <- which.max( AUC.test$TR[,i] )
	Pred_To_Inc <- rownames(COEF.all)[WHICH_TO_INC]
	print(paste( "Include",WHICH_TO_INC,"-",Pred_To_Inc ))
	## Compile all the results for Plotting
	AUC[i,"TR"] <- AUC.test$TR[WHICH_TO_INC,i]
	AUC[i,"TS"] <- AUC.test$TS[WHICH_TO_INC,i]
	## Include, but don't Test, that Predictor from now on
	TO_INCL <- c( TO_INCL, Pred_To_Inc )
	TO_TEST <- setdiff( TO_TEST, WHICH_TO_INC )
	## Output Status
	if (i%%5==0) { 
		# print( tail(Which_Pred) )
		print(paste( "Finished",i,"of",MOD_SIZE,"-",round( (proc.time()-start_time)[3], 1) ))
	}
}
rownames(AUC) <- TO_INCL
colnames(AUC.test$TS) <- colnames(AUC.test$TR) <- TO_INCL
head(AUC,20)
AUC.1 <- AUC
AUC.1.test <- AUC.test

## Plot AUC.1 Values for Test/Training Sets
XLIM <- c(1,nrow(AUC.1))
YLIM <- c(0,1) # range(AUC.1)
COLS <- c( rep(c("chartreuse2","dodgerblue2"),2) ) # , rep("dodgerblue2",2) )
LTYS <- c( 1,1,2,2 )
jpeg( paste(PathToSave,"MOD_Scr_",File_Name,".jpeg",sep=""), height=1000, width=2000, pointsize=30 )
plot( 0,0, type="n", xlim=XLIM, ylim=YLIM, xaxt="n", xlab="Predictors Included (Cummulative)", ylab="AUC", main="AUC for Training/Test Sets (Score Model)" )
axis( 1, at=1:nrow(AUC.1), label=TO_INCL, las=2 )
abline( h=seq( 0,1,.1), lty=c(rep(2,5),1,rep(2,5)), col=c(rep("firebrick3",5),"black",rep("grey50",5)) ) 
abline( v=seq( 0,nrow(AUC.1),5), lty=2, col="grey50" ) 
for ( i in 1:2 ) {
	points( 1:nrow(AUC.1), AUC.1[,i], col=COLS[i], lty=LTYS[i], type="l", lwd=3 )
}
legend( XLIM[1], YLIM[2], legend=colnames(AUC.1), lty=LTYS, col=COLS, lwd=3 )
dev.off()




###############################################################
## SAVE DATA ##################################################
###############################################################

## Compile all the data I want to save
COMPILE <- list( AUC.9, AUC.9.test, AUC.1, AUC.1.test )
names(COMPILE) <- c("AUC_9","AUC_9_test" ,"AUC_1","AUC_1_test")

## Save all the data
save( COMPILE, file=paste(PathToSave,"MOD_",File_Name,".Rdata",sep="") )
















###############################################################
## END OF DOC #################################################
###############################################################
