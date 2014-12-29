## 12-Model.R ##
## R function to Use Training Results and Build Model ##
## CONTINUOUS Phenotypes ##
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
# LINE <- c( paste( "/projects/janssen/Poly_Train/20141124test_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/", c("LT8_DEL_MNe_MN", "Cov_w_PCs.txt", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2_012.raw"), sep="" ), "DAS_BL_MN,PC1,PC2" )
# LINE <- c( paste( "/Users/kstandis/Downloads/20141124_Poly_Train/", c("LT8_DEL_MNe_MN", "Cov_w_PCs.txt", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2_012.raw"), sep="" ), "DAS_BL_MN,PC1,PC2" )
# LINE <- c( paste( "/projects/janssen/Poly_Train/20141222_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2/", c("LT8_DEL_MNe_MN", "Cov_w_PCs.txt", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2.compile.short", "CND_LT8_DEL_MNe_MN_DAS_BL_MN_PC1_PC2_012.raw"), sep="" ), "DAS_BL_MN,PC1,PC2" )

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
BETA.gt.1 <- BETA.gt.a[ -which( names(BETA.gt.a) %in% RM.corr ) ] #; names(BETA.gt) <- KP.P.id[ -which( names(BETA.gt.a) %in% RM.corr ) ]
P.gt.1 <- P.gt.a[ -which( names(P.gt.a) %in% RM.corr ) ] #; names(P.gt) <- KP.P.id[ -which( names(P.gt.a) %in% RM.corr ) ]
GT.cand.1 <- GT.cand.a[ ,-which( colnames(GT.cand.a) %in% RM.corr ) ]
GT.corr.1 <- cor( GT.cand.1 )
## Compare
jpeg( paste(PathToSave,"MOD_SNP_Corr",File_Name,".jpeg",sep=""), height=1000, width=2500, pointsize=30 )
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
 # and use BETA values calculated each time on TEST set

## Retreive Coefficients for Covariates Included in Model
TAB.TR.9.cov <- MG.TR.9[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
MOD.cov <- lm( Pheno ~ . , data=TAB.TR.9.cov )
BETA.cov <- coef( MOD.cov )
P.cov <- summary(MOD.cov)$coefficients[,"Pr(>|t|)"]

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt.9 )
P.all <- c( P.cov[ORD.cov], P.gt.9[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt.9)[ORD.gt] )
BETA.all <- c( BETA.cov[ORD.cov], BETA.gt.9[ORD.gt] ) ; names(BETA.all) <- c( names(BETA.cov)[ORD.cov], names(BETA.gt.9)[ORD.gt] )
COEF.all <- data.frame( P.all, BETA.all ) ; rownames(COEF.all) <- names(BETA.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]

## Loop through Predictors and add each one iteratively to model
 # Then assess fit of model
TO_TEST.cov <- 1:(length(BETA.cov)-1)
TO_TEST.gt <- length(TO_TEST.cov) + 1:length(BETA.gt.9)
TO_TEST <- TO_TEST.cov
TO_INCL <- c()
TO_OMIT <- c()
MOD_SIZE <- 60 ; MOD_SIZE <- min( MOD_SIZE, nrow(COEF.all) )
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
		TEMP.TR.9 <- data.matrix( MG.TR.9[, Which_Pred ] )
		TEMP.TS.9 <- data.matrix( MG.TS.9[, Which_Pred ] )
		## Fit Model Using lm()
		MOD.TR.9 <- lm( MG.TR.9[,"Pheno"] ~ TEMP.TR.9 )
		COEF.mod <- coef(MOD.TR.9)
		## Predict Values from Test Set based on Model with these Predictors
		PRED.TR.9 <- COEF.mod[1] + TEMP.TR.9 %*% matrix(COEF.mod[-1])
		PRED.TS.9 <- COEF.mod[1] + TEMP.TS.9 %*% matrix(COEF.mod[-1])
		## Calculate R2 for covariate model
		 # Calc Sum Squares
		DIFF.TS.9.mod <- PRED.TS.9 - MG.TS.9$Pheno
		DIFF.TR.9.mod <- PRED.TR.9 - MG.TR.9$Pheno
		SS.TS.9.mod <- sum(DIFF.TS.9.mod^2) - sum(DIFF.TS.9.mod)^2 / length(DIFF.TS.9.mod)
		SS.TR.9.mod <- sum(DIFF.TR.9.mod^2) - sum(DIFF.TR.9.mod)^2 / length(DIFF.TR.9.mod)
		SS.TS.9.val <- sum(MG.TS.9$Pheno^2) - sum(MG.TS.9$Pheno)^2 / nrow(MG.TS.9)
		SS.TR.9.val <- sum(MG.TR.9$Pheno^2) - sum(MG.TR.9$Pheno)^2 / nrow(MG.TR.9)
		 # R2
		R2.TS.9.all <- 1 - SS.TS.9.mod / SS.TS.9.val
		R2.TR.9.all <- 1 - SS.TR.9.mod / SS.TR.9.val
		 # Adjusted R2
		R2.adj.TS.9.all <- 1 - ( SS.TS.9.mod/(nrow(MG.TS.9)-length(Which_Pred)) / (SS.TS.9.val/(nrow(MG.TS.9)-1)) )
		R2.adj.TR.9.all <- 1 - ( SS.TR.9.mod/(nrow(MG.TR.9)-length(Which_Pred)) / (SS.TR.9.val/(nrow(MG.TR.9)-1)) )
		 # From Mod
		# R2.TR.9.all <- summary(MOD.TR.9)$r.squared
		# R2.adj.TR.9.all <- summary(MOD.TR.9)$adj.r.squared
		 # Compile Results to Compare for this Iteration
		R2.test$TS[r,i] <- R2.TS.9.all
		R2.test$TR[r,i] <- R2.TR.9.all
		R2.test$TS.adj[r,i] <- R2.adj.TS.9.all
		R2.test$TR.adj[r,i] <- R2.adj.TR.9.all
	}
	## Decide Which Model is the Best
	WHICH_TO_INC <- which.max( R2.test$TR[,i] )
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
R2.9 <- R2
R2.9.test <- R2.test

## Plot R2.9 Values for Test/Training Sets
XLIM <- c(1,nrow(R2.9))
YLIM <- c(-1,1) # range(R2.9)
COLS <- c( rep(c("chartreuse2","dodgerblue2"),2) ) # , rep("dodgerblue2",2) )
LTYS <- c( 1,1,2,2 )
jpeg( paste(PathToSave,"MOD_Reg_",File_Name,".jpeg",sep=""), height=1000, width=2000, pointsize=30 )
plot( 0,0, type="n", xlim=XLIM, ylim=YLIM, xaxt="n", xlab="Predictors Included (Cummulative)", ylab="R2.9 & Adj R2.9", main="Variance Explained by Model (Regression Model)" )
axis( 1, at=1:nrow(R2.9), label=TO_INCL, las=2 )
abline( h=seq( -1,1,.1), lty=c(rep(2,10),1,rep(2,10)), col=c(rep("firebrick3",10),"black",rep("grey50",10)) ) 
abline( v=seq( 0,nrow(R2.9),5), lty=2, col="grey50" ) 
for ( i in 1:4 ) {
	points( 1:nrow(R2.9), R2.9[,i], col=COLS[i], lty=LTYS[i], type="l", lwd=3 )
}
legend( quantile(XLIM,.8), quantile(YLIM,.9), legend=colnames(R2.9), lty=LTYS, col=COLS, lwd=3 )
dev.off()

# ## Heatmaps of Predictors
# COLS <- c( "gold1","chocolate2","firebrick3","black","slateblue3","steelblue2","springgreen1" )
# HT_COLS <- colorRampPalette(COLS)(50)
# ## Correlation b/n Predictors
# COR.INCL <- cor( MG.TR[, TO_INCL], use="pairwise.complete.obs", method="spearman" ) 
# BRKS <- seq(-1,1,length.out=51)
# heatmap.2( COR.INCL, col=HT_COLS, breaks=BRKS, scale="none", Colv=F, Rowv=T, dendrogram="row", trace="none" )
# ## Heatmap of Predictor Tests (R2.9.test$TR)
# BRKS <- c( seq(min(R2.9.test$TR,na.rm=T),0,length.out=26)[1:25], 0, seq(0,5,length.out=26)[2:26] )
# heatmap.2( R2.9.test$TR, col=HT_COLS, breaks=BRKS, scale="column", Colv=F, Rowv=F, dendrogram="none", trace="none" )
# heatmap.2( t(R2.9.test$TR), col=HT_COLS, scale="row", Colv=F, Rowv=F, dendrogram="none", trace="none" )


###############################################################
## MODEL w/ +/- SUM of GT COUNTS ##############################
###############################################################
## Using <.1 correlated genotypes:
 # calculate sum (incl. +/- for BETA) of genotypes person to yield a score
 # fit that single score to response data
 # and use BETA value for score calculated each time on TEST set

## Retreive Coefficients for Covariates Included in Model
TAB.TR.1.cov <- MG.TR.1[, c("Pheno",Cov_Cols) ]# gsub( ",","+",Covs_Command )
MOD.cov <- lm( Pheno ~ . , data=TAB.TR.1.cov )
BETA.cov <- coef( MOD.cov )
P.cov <- summary(MOD.cov)$coefficients[,"Pr(>|t|)"]

## Compile P-Values and Beta Values for Covariates and Genotypes
ORD.cov <- order( P.cov ) ; ORD.gt <- order( P.gt.1 )
P.all <- c( P.cov[ORD.cov], P.gt.1[ORD.gt] ) ; names(P.all) <- c( names(P.cov)[ORD.cov], names(P.gt.1)[ORD.gt] )
BETA.all <- c( BETA.cov[ORD.cov], BETA.gt.1[ORD.gt] ) ; names(BETA.all) <- c( names(BETA.cov)[ORD.cov], names(BETA.gt.1)[ORD.gt] )
COEF.all <- data.frame( P.all, BETA.all ) ; rownames(COEF.all) <- names(BETA.all)
COEF.int <- COEF.all[ which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all <- COEF.all[ -which(rownames(COEF.all)=="(Intercept)"), ]
COEF.all$BETA.all[(length(BETA.cov)):nrow(COEF.all)] <- c(-1,1)[factor(COEF.all$BETA.all[(length(BETA.cov)):nrow(COEF.all)]>0)]

## Loop through Predictors and add each one iteratively to model
 # Then assess fit of model
TO_TEST.cov <- 1:(length(BETA.cov)-1)
TO_TEST.gt <- length(TO_TEST.cov) + 1:length(BETA.gt.1)
TO_TEST <- TO_TEST.cov
TO_INCL <- TO_OMIT <- c()
GT_Flag <- 0
MOD_SIZE <- 30 ; MOD_SIZE <- min( MOD_SIZE, nrow(COEF.all) ) # 20
# MOD_SIZE <- nrow(COEF.all) # 20
N_PREDS <- nrow(COEF.all)
R2.test <- list()
R2.test$TS <- R2.test$TR <- R2.test$TS.adj <- R2.test$TR.adj <- array( , c(N_PREDS,MOD_SIZE) )
rownames(R2.test$TS) <- rownames(R2.test$TR) <- rownames(R2.test$TS.adj) <- rownames(R2.test$TR.adj) <- rownames(COEF.all)[1:N_PREDS]
SCOR.test <- list()
SCOR.test$P <- SCOR.test$BETA <- array( , c(N_PREDS,MOD_SIZE) )
rownames(SCOR.test$P) <- rownames(SCOR.test$BETA) <- rownames(COEF.all)[1:N_PREDS]
R2 <- array( , dim=c(MOD_SIZE,4) )
colnames(R2) <- c("TR","TS","Adj.TR","Adj.TS")
start_time <- proc.time()
for ( i in 1:MOD_SIZE ) {
# for ( i in 1:10 ) {
	## After Covariates, Switch to Genotypes
	if ( i == length(TO_TEST.cov)+1 ) {
		TO_TEST <- TO_TEST.gt[1:(N_PREDS-length(TO_TEST.cov))]
		GT_Flag <- 1
	}
	## Loop through all Predictors and Find Best one to ADD to model
	for ( r in TO_TEST ) {
		This_Pred <- rownames(COEF.all)[r]
		Which_Pred <- c( TO_INCL, This_Pred )
		Which_GT <- setdiff( Which_Pred, names(BETA.cov) )
		Which_Cov <- setdiff( Which_Pred, Which_GT )
		## Create Binary Matrix for Genotypes
		BIN.TR.1 <- data.matrix( MG.TR.1[, Which_GT ] ) ; BIN.TR.1[which(BIN.TR.1==2,arr.ind=T)] <- 1
		BIN.TS.1 <- data.matrix( MG.TS.1[, Which_GT ] ) ; BIN.TS.1[which(BIN.TS.1==2,arr.ind=T)] <- 1
		## Calculate Score for each Person
		SCOR.TR.1 <- BIN.TR.1 %*% COEF.all[Which_GT,"BETA.all"]
		SCOR.TS.1 <- BIN.TS.1 %*% COEF.all[Which_GT,"BETA.all"]
		## Build Temporary Data Frame with These Covariates
		if (GT_Flag==1) {
			TEMP.TR.1 <- data.matrix(data.frame( MG.TR.1[, c(Which_Cov) ], SCOR.TR.1 ))
			TEMP.TS.1 <- data.matrix(data.frame( MG.TS.1[, c(Which_Cov) ], SCOR.TS.1 ))
			}else{
			TEMP.TR.1 <- data.matrix( MG.TR.1[, c(Which_Cov) ] )
			TEMP.TS.1 <- data.matrix( MG.TS.1[, c(Which_Cov) ] )
		}
		## Fit Model Using lm()
		MOD.TR.1 <- lm( MG.TR.1[,"Pheno"] ~ TEMP.TR.1 )
		COEF.mod <- coef(MOD.TR.1)
		PRED.TR.1 <- COEF.mod[1] + TEMP.TR.1 %*% matrix(COEF.mod[-1])
		PRED.TS.1 <- COEF.mod[1] + TEMP.TS.1 %*% matrix(COEF.mod[-1])
		## Calculate R2 for covariate model
		 # Calc Sum Squares
		DIFF.TS.1.mod <- PRED.TS.1 - MG.TS.1$Pheno
		DIFF.TR.1.mod <- PRED.TR.1 - MG.TR.1$Pheno
		SS.TS.1.mod <- sum(DIFF.TS.1.mod^2) - sum(DIFF.TS.1.mod)^2 / length(DIFF.TS.1.mod)
		SS.TR.1.mod <- sum(DIFF.TR.1.mod^2) - sum(DIFF.TR.1.mod)^2 / length(DIFF.TR.1.mod)
		SS.TS.1.val <- sum(MG.TS.1$Pheno^2) - sum(MG.TS.1$Pheno)^2 / nrow(MG.TS.1)
		SS.TR.1.val <- sum(MG.TR.1$Pheno^2) - sum(MG.TR.1$Pheno)^2 / nrow(MG.TR.1)
		 # R2
		R2.TS.1.all <- 1 - SS.TS.1.mod / SS.TS.1.val
		R2.TR.1.all <- 1 - SS.TR.1.mod / SS.TR.1.val
		 # Adjusted R2
		R2.adj.TS.1.all <- 1 - ( SS.TS.1.mod/(nrow(MG.TS.1)-length(Which_Pred)) / (SS.TS.1.val/(nrow(MG.TS.1)-1)) )
		R2.adj.TR.1.all <- 1 - ( SS.TR.1.mod/(nrow(MG.TR.1)-length(Which_Pred)) / (SS.TR.1.val/(nrow(MG.TR.1)-1)) )
		 # From Mod
		# R2.TR.1.all <- summary(MOD.TR.1)$r.squared
		# R2.adj.TR.1.all <- summary(MOD.TR.1)$adj.r.squared
		 # Compile Results to Compare for this Iteration
		R2.test$TS[r,i] <- R2.TS.1.all
		R2.test$TR[r,i] <- R2.TR.1.all
		R2.test$TS.adj[r,i] <- R2.adj.TS.1.all
		R2.test$TR.adj[r,i] <- R2.adj.TR.1.all
		SCOR.test$P[r,i] <- summary(MOD.TR.1)$coefficients[length(COEF.mod),"Pr(>|t|)"]
		SCOR.test$BETA[r,i] <- COEF.mod[length(COEF.mod)]
	}
	## Decide Which Model is the Best
	WHICH_TO_INC <- which.max( R2.test$TR[,i] )
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
R2.1 <- R2
R2.1.test <- R2.test

## Plot R2.1 Values for Test/Training Sets
XLIM <- c(1,nrow(R2.1))
YLIM <- c(-1,1) # range(R2.1)
COLS <- c( rep(c("chartreuse2","dodgerblue2"),2) ) # , rep("dodgerblue2",2) )
LTYS <- c( 1,1,2,2 )
jpeg( paste(PathToSave,"MOD_Scr_",File_Name,".jpeg",sep=""), height=1000, width=2000, pointsize=30 )
plot( 0,0, type="n", xlim=XLIM, ylim=YLIM, xaxt="n", xlab="Predictors Included (Cummulative)", ylab="R2.1 & Adj R2.1", main="Variance Explained by Model (Score Model)" )
axis( 1, at=1:nrow(R2.1), label=TO_INCL, las=2 )
abline( h=seq( -1,1,.1), lty=c(rep(2,10),1,rep(2,10)), col=c(rep("firebrick3",10),"black",rep("grey50",10)) ) 
abline( v=seq( 0,nrow(R2.1),5), lty=2, col="grey50" ) 
for ( i in 1:4 ) {
	points( 1:nrow(R2.1), R2.1[,i], col=COLS[i], lty=LTYS[i], type="l", lwd=3 )
}
legend( quantile(XLIM,.01), quantile(YLIM,.9), legend=colnames(R2.1), lty=LTYS, col=COLS, lwd=3 )
dev.off()

###############################################################
## SAVE DATA ##################################################
###############################################################

## Compile all the data I want to save
COMPILE <- list( R2.9, R2.9.test, R2.1, R2.1.test )
names(COMPILE) <- c("R2_9","R2_9_test","R2_1","R2_1_test")

## Save all the data
save( COMPILE, file=paste(PathToSave,"MOD_",File_Name,sep="") )







###############################################################
## END OF DOC #################################################
###############################################################