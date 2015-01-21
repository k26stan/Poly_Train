## Compile Results from many Polygenic Models ##
## From Janssen RA Cohort - Q3 Deliverables ##
## January 12, 2015 ##
## Kristopher Standish ##

#########################################################################
## TRANSFER DATA ONTO MAC ###############################################
#########################################################################

## Specify inputs/names/etc...
DATE=20150119
N=20
# PHENO_NAME=ACR50_100wk_AGE_SEX
PHENO_NAME=LT8_DEL_MNe_MN_DAS_BL_MN

mkdir /Users/kstandis/Data/Burn/Poly_Train/${DATE}_Meta/
mkdir /Users/kstandis/Data/Burn/Poly_Train/${DATE}_Meta/Data/
mkdir /Users/kstandis/Data/Burn/Poly_Train/${DATE}_Meta/Plots/
# On Mac Terminal
for num in `seq 1 ${N}`
do
	# Change TAG Identifier
	TAG=${DATE}_${num}
	# Set Path To Mac Directory
	MAC_PATH=/Users/kstandis/Data/Burn/Poly_Train/${DATE}_Meta/Data/${TAG}_MOD_${PHENO_NAME}.Rdata
	# Import Data
	SCP_IN /projects/janssen/Poly_Train/${TAG}_${PHENO_NAME}/*Rdata ${MAC_PATH}
	# SCP_IN /projects/janssen/Deliver_3/Polygenic_Modeling/20150105_${num}_ACR50_100wk_AGE_SEX/MOD_ACR50_100wk_AGE_SEX.Rdata /Users/kstandis/Data/Burn/Poly_Train/20150112_Meta/Data/20150105_${num}_MOD_ACR50_100wk_AGE_SEX.Rdata
	# SCP_IN /projects/janssen/Deliver_3/Polygenic_Modeling/20150105_${num}_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata /Users/kstandis/Data/Burn/Poly_Train/20150112_Meta/Data/20150105_${num}_MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
done

## Create Directory for Results
mkdir Dropbox/Schork/JNJ11/Association/Poly_Train/${DATE}_${PHENO_NAME}/
## Loop through N times
for num in `seq 1 ${N}`
do
	# Change TAG Identifier
	TAG=${DATE}_${num}
	# Set Path To Mac Directory
	MAC_DIR=Dropbox/Schork/JNJ11/Association/Poly_Train/${DATE}_${PHENO_NAME}/${num}/
	# Create Mac Directory
	mkdir ${MAC_DIR}
	SCP_IN /projects/janssen/Poly_Train/${TAG}_${PHENO_NAME}/*jpeg ${MAC_DIR}
done

#########################################################################
## LOAD DATA ############################################################
#########################################################################

DATE <- "20150119"

## Set Up Paths (Mac)
PathToData12 <- "/Users/kstandis/Data/Burn/Poly_Train/20150112_Meta/Data/"
PathToData19 <- "/Users/kstandis/Data/Burn/Poly_Train/20150119_Meta/Data/"
PathToSave <- paste("/Users/kstandis/Data/Burn/Poly_Train/",DATE,"_Meta/Plots/",sep="")

## Load R Data (DEL - 1/12)
DEL.12 <- list()
for ( i in 1:10 ) {
	PATH_DEL <- paste(PathToData12,"20150105_",i,"_MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata",sep="")
	if ( file.exists(PATH_DEL) ) {
		load( PATH_DEL )
		DEL.12[[paste("Run",i,sep="_")]] <- COMPILE
	}
}
## Load R Data (ACR - 1/12)
ACR.12 <- list()
for ( i in 1:10 ) {
	PATH_ACR <- paste(PathToData12,"20150105_",i,"_MOD_ACR50_100wk_AGE_SEX.Rdata",sep="")
	if ( file.exists(PATH_ACR) ) {
		load( PATH_ACR )
		ACR.12[[paste("Run",i,sep="_")]] <- COMPILE
	}
}
## Load R Data (DEL - 1/12)
DEL.19 <- list()
for ( i in c(1,3:5,7,9:20) ) {
	PATH_DEL <- paste(PathToData19,"20150119_",i,"_MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata",sep="")
	if ( file.exists(PATH_DEL) ) {
		load( PATH_DEL )
		DEL.19[[paste("Run",i,sep="_")]] <- COMPILE
	}
}

## Specify Covariates for each Model
COVS.del <- c("DAS_BL_MN")
COVS.acr <- c("AGE","SEX")

## Compile Runs for Analysis
DEL <- DEL.19
ACR <- ACR.12

#########################################################################
## PULL OUT VARIANTS ####################################################
#########################################################################

##################################################
## CONTINUOUS: DEL_MNe_MN Results
RUNS.del <- names(DEL)

## Set up variables to compile
 # List of variants in each model
VARS.9.del <- VARS.1.del <- list()
 # TS/TR Improvement for inclusion of each variant
R2UP.9.del <- R2UP.1.del <- list()

## Loop through Models & Compile Stats
r <- 0
for ( run in RUNS.del ) {
	r <- r + 1
	## Get variant IDs
	VARS.9.del[[run]] <- rownames( DEL[[run]]$R2_9 )[ (length(COVS.del)+1) : nrow(DEL[[run]]$R2_9) ]
	VARS.1.del[[run]] <- rownames( DEL[[run]]$R2_1 )[ (length(COVS.del)+1) : nrow(DEL[[run]]$R2_1) ]
	## Get change in R2 from inclusion of each variant
	len.9 <- length(VARS.9.del[[run]])
	len.1 <- length(VARS.1.del[[run]])
	R2UP.9.del[[run]] <- array( , c(len.9,2) )
	R2UP.1.del[[run]] <- array( , c(len.1,2) )
	colnames(R2UP.9.del[[run]]) <- colnames(R2UP.1.del[[run]]) <- c("TR","TS")
	rownames(R2UP.9.del[[run]]) <- VARS.9.del[[run]]
	rownames(R2UP.1.del[[run]]) <- VARS.1.del[[run]]
	R2UP.9.del[[run]][,"TR"] <- DEL[[run]]$R2_9[ VARS.9.del[[run]], "TR" ] - DEL[[run]]$R2_9[ which(rownames(DEL[[run]]$R2_9) %in% VARS.9.del[[run]])-1, "TR" ]
	R2UP.9.del[[run]][,"TS"] <- DEL[[run]]$R2_9[ VARS.9.del[[run]], "TS" ] - DEL[[run]]$R2_9[ which(rownames(DEL[[run]]$R2_9) %in% VARS.9.del[[run]])-1, "TS" ]
	R2UP.1.del[[run]][,"TR"] <- DEL[[run]]$R2_1[ VARS.1.del[[run]], "TR" ] - DEL[[run]]$R2_1[ which(rownames(DEL[[run]]$R2_1) %in% VARS.1.del[[run]])-1, "TR" ]
	R2UP.1.del[[run]][,"TS"] <- DEL[[run]]$R2_1[ VARS.1.del[[run]], "TS" ] - DEL[[run]]$R2_1[ which(rownames(DEL[[run]]$R2_1) %in% VARS.1.del[[run]])-1, "TS" ]
}

##################################################
## BINARY: ACR50_100wk Results
RUNS.acr <- names(ACR)

## Set up variables to compile
 # List of variants in each model
VARS.9.acr <- VARS.1.acr <- list()
 # TS/TR Improvement for inclusion of each variant
AUCUP.9.acr <- AUCUP.1.acr <- list()

## Loop through Models & Compile Stats
r <- 0
for ( run in RUNS.acr ) {
	r <- r + 1
	## Get variant IDs
	VARS.9.acr[[run]] <- rownames( ACR[[run]]$AUC_9 )[ (length(COVS.acr)+1) : nrow(ACR[[run]]$AUC_9) ]
	# VARS.1.acr[[run]] <- rownames( ACR[[run]]$AUC_1 )[ (length(COVS.acr)+1) : nrow(ACR[[run]]$AUC_1) ]
	## Get change in AUC from inclusion of each variant
	len.9 <- length(VARS.9.acr[[run]])
	# len.1 <- length(VARS.1.acr[[run]])
	AUCUP.9.acr[[run]] <- array( , c(len.9,2) )
	AUCUP.1.acr[[run]] <- array( , c(len.1,2) )
	colnames(AUCUP.9.acr[[run]]) <- colnames(AUCUP.1.acr[[run]]) <- c("TR","TS")
	rownames(AUCUP.9.acr[[run]]) <- VARS.9.acr[[run]]
	# rownames(AUCUP.1.acr[[run]]) <- VARS.1.acr[[run]]
	AUCUP.9.acr[[run]][,"TR"] <- ACR[[run]]$AUC_9[ VARS.9.acr[[run]], "TR" ] - ACR[[run]]$AUC_9[ (length(COVS.acr)):(length(COVS.acr)+len.9-1), "TR" ]
	AUCUP.9.acr[[run]][,"TS"] <- ACR[[run]]$AUC_9[ VARS.9.acr[[run]], "TS" ] - ACR[[run]]$AUC_9[ (length(COVS.acr)):(length(COVS.acr)+len.9-1), "TS" ]
	# AUCUP.1.acr[[run]][,"TR"] <- ACR[[run]]$AUC_1[ VARS.1.acr[[run]], "TR" ] - ACR[[run]]$AUC_1[ (length(COVS.acr)):(len.1), "TR" ]
	# AUCUP.1.acr[[run]][,"TS"] <- ACR[[run]]$AUC_1[ VARS.1.acr[[run]], "TS" ] - ACR[[run]]$AUC_1[ (length(COVS.acr)):(len.1), "TS" ]
}

#########################################################################
## COMPILE RESULTS ######################################################
#########################################################################

##################################################
## CONTINUOUS: DEL_MNe_MN Results
 # Regression Model
ALL.VARS.9.del <- Reduce( union, VARS.9.del )
RANKS.9.del <- WTS.9.del <- array( , c(length(ALL.VARS.9.del),length(RUNS.del)) )
rownames(RANKS.9.del) <- rownames(WTS.9.del) <- ALL.VARS.9.del
colnames(RANKS.9.del) <- colnames(WTS.9.del) <- RUNS.del
for ( v in 1:length(ALL.VARS.9.del) ) {
	var <- ALL.VARS.9.del[v]
	for ( run in RUNS.del ) {
		if ( var %in% VARS.9.del[[run]] ) {
			RANKS.9.del[var,run] <- which( VARS.9.del[[run]] == var )
			WTS.9.del[var,run] <- R2UP.9.del[[run]][var,"TS"]
		}
	}
}
 # Score Model
ALL.VARS.1.del <- Reduce( union, VARS.1.del )
RANKS.1.del <- WTS.1.del <- array( , c(length(ALL.VARS.1.del),length(RUNS.del)) )
rownames(RANKS.1.del) <- rownames(WTS.1.del) <- ALL.VARS.1.del
colnames(RANKS.1.del) <- colnames(WTS.1.del) <- RUNS.del
for ( v in 1:length(ALL.VARS.1.del) ) {
	var <- ALL.VARS.1.del[v]
	for ( run in RUNS.del ) {
		if ( var %in% VARS.1.del[[run]] ) {
			RANKS.1.del[var,run] <- which( VARS.1.del[[run]] == var )
			WTS.1.del[var,run] <- R2UP.1.del[[run]][var,"TS"]
		}
	}
}

##################################################
## BINARY: ACR50_100wk Results
 # Regression Model
ALL.VARS.9.acr <- Reduce( union, VARS.9.acr )
RANKS.9.acr <- WTS.9.acr <- array( , c(length(ALL.VARS.9.acr),length(RUNS.acr)) )
rownames(RANKS.9.acr) <- rownames(WTS.9.acr) <- ALL.VARS.9.acr
colnames(RANKS.9.acr) <- colnames(WTS.9.acr) <- RUNS.acr
for ( v in 1:length(ALL.VARS.9.acr) ) {
	var <- ALL.VARS.9.acr[v]
	for ( run in RUNS.acr ) {
		if ( var %in% VARS.9.acr[[run]] ) {
			RANKS.9.acr[var,run] <- which( VARS.9.acr[[run]] == var )
			WTS.9.acr[var,run] <- AUCUP.9.acr[[run]][var,"TS"]
		}
	}
}
#  # Score Model
# ALL.VARS.1.acr <- Reduce( union, VARS.1.acr )
# RANKS.1.acr <- WTS.1.acr <- array( , c(length(ALL.VARS.1.acr),length(RUNS.acr)) )
# rownames(RANKS.1.acr) <- rownames(WTS.1.acr) <- ALL.VARS.1.acr
# colnames(RANKS.1.acr) <- colnames(WTS.1.acr) <- RUNS.acr
# for ( v in 1:length(ALL.VARS.1.acr) ) {
# 	var <- ALL.VARS.1.acr[v]
# 	for ( run in RUNS.acr ) {
# 		if ( var %in% VARS.1.acr[[run]] ) {
# 			RANKS.1.acr[var,run] <- which( VARS.1.acr[[run]] == var )
# 			WTS.1.acr[var,run] <- R2UP.1.acr[[run]][var,"TS"]
# 		}
# 	}
# }

#########################################################################
## PRIORITIZE RESULTS ###################################################
#########################################################################

##################################################
## CONTINUOUS: DEL_MNe_MN Results

 # Count how many models each variant appears in
COUNT.9.del <- apply( RANKS.9.del, 1, function(x) length(which(!is.na(x))) )
COUNT.1.del <- apply( RANKS.1.del, 1, function(x) length(which(!is.na(x))) )
 # Get weight based on R2 improvement
WT.9.del <- apply( WTS.9.del, 1, sum, na.rm=T )
WT.1.del <- apply( WTS.1.del, 1, sum, na.rm=T )
 # Get best rank across models for each variant
MIN.9.del <- apply( RANKS.9.del, 1, min, na.rm=T )
MIN.1.del <- apply( RANKS.1.del, 1, min, na.rm=T )
 # Calculate the MEAN rank of each variant
MEAN.9.del <- apply( RANKS.9.del, 1, mean, na.rm=T )
MEAN.1.del <- apply( RANKS.1.del, 1, mean, na.rm=T )
 # Take MEAN rank and divide by number of runs in which variant appears
 # This boosts variants that appear in multiple models by giving a lower rank STAT
STAT.9.del <- MIN.9.del / COUNT.9.del
STAT.1.del <- MIN.1.del / COUNT.1.del

## Compile into single table
COMP.9.del <- data.frame( COUNT=COUNT.9.del, MEAN_RNK=round(MEAN.9.del,2), STAT=round(STAT.9.del,2), R2_SUM=round(WT.9.del,2), RANKS.9.del )
COMP.9.del <- COMP.9.del[order(COMP.9.del[,"STAT"]),]
COMP.9.del.R2 <- COMP.9.del[order(COMP.9.del[,"R2_SUM"],decreasing=T),]
COMP.1.del <- data.frame( COUNT=COUNT.1.del, MEAN_RNK=round(MEAN.1.del,2), STAT=round(STAT.1.del,2), R2_SUM=round(WT.1.del,2), RANKS.1.del )
COMP.1.del <- COMP.1.del[order(COMP.1.del[,"STAT"]),]
COMP.1.del.R2 <- COMP.1.del[order(COMP.1.del[,"R2_SUM"],decreasing=T),]

## Heatmap
library(gplots)
COLS.list <- c("springgreen1","steelblue1","slateblue2","black")
COLS <- colorRampPalette(COLS.list)(100)
heatmap.2( t(data.matrix(COMP.1.del)), col=COLS, scale="none", dendrogram="none", Colv=F, Rowv=F, trace="none" )
heatmap.2( t(data.matrix(COMP.9.del)), col=COLS, scale="none", dendrogram="none", Colv=F, Rowv=F, trace="none" )

## Save Tables
write.table( COMP.1.del, paste(PathToSave,"TAB_1_DEL.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )
write.table( COMP.9.del, paste(PathToSave,"TAB_9_DEL.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F )

##################################################
## BINARY: ACR50_100wk Results

 # Count how many models each variant appears in
COUNT.9.acr <- apply( RANKS.9.acr, 1, function(x) length(which(!is.na(x))) )
# COUNT.1.acr <- apply( RANKS.1.acr, 1, function(x) length(which(!is.na(x))) )
 # Get weight based on R2 improvement
WT.9.acr <- apply( WTS.9.acr, 1, sum, na.rm=T )
# WT.1.acr <- apply( WTS.1.acr, 1, sum, na.rm=T )
 # Get best rank across models for each variant
MIN.9.acr <- apply( RANKS.9.acr, 1, min, na.rm=T )
# MIN.1.acr <- apply( RANKS.1.acr, 1, min, na.rm=T )
 # Calculate the MEAN rank of each variant
MEAN.9.acr <- apply( RANKS.9.acr, 1, mean, na.rm=T )
# MEAN.1.acr <- apply( RANKS.1.acr, 1, mean, na.rm=T )
 # Take MEAN rank and divide by number of runs in which variant appears
 # This boosts variants that appear in multiple models by giving a lower rank STAT
STAT.9.acr <- MIN.9.acr / COUNT.9.acr
# STAT.1.acr <- MEAN.1.acr / COUNT.1.acr

## Compile into single table
COMP.9.acr <- data.frame( COUNT=COUNT.9.acr, MEAN_RNK=round(MEAN.9.acr,2), STAT=round(STAT.9.acr,2), AUC_SUM=round(WT.9.acr,2), RANKS.9.acr )
COMP.9.acr <- COMP.9.acr[order(COMP.9.acr[,"STAT"]),]
COMP.9.acr.AUC <- COMP.9.acr[order(COMP.9.acr[,"AUC_SUM"],decreasing=T),]
# COMP.1.acr <- data.frame( COUNT=COUNT.1.acr, MEAN_RNK=round(MEAN.1.acr,2), STAT=round(STAT.1.acr,2), AUC_SUM=round(WT.1.acr,2), RANKS.1.acr )
# COMP.1.acr <- COMP.1.acr[order(COMP.1.acr[,"STAT"]),]
# COMP.1.acr.AUC <- COMP.1.acr[order(COMP.1.acr[,"AUC_SUM"],decreasing=T),]













#########################################################################
## END OF DOC ###########################################################
#########################################################################