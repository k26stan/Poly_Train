DATE=20150119_2
HOME_DIR=/projects/janssen/Poly_Train
cd ${HOME_DIR}

## Files
VAR_FILE=${HOME_DIR}/../scripts/AS-ASSOCIATION/BED_FILES/BED_HC_FULL_SNP.bed
ANNOTS=${HOME_DIR}/../ANNOTATE/JnJ_121613_all_annotations.txt.bgz
PHENO_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/ACR50_100wk.txt
PHENO_TYPE=B
COV_FILE=${HOME_DIR}/../ASSOCIATION/PH-PHENOTYPES/COV.txt
COVS=`echo AGE SEX`
EIG_VEC=${HOME_DIR}/../ASSOCIATION/EIGEN/HC_FULL.eigenvec
PC_COUNT=0
START_STEP=12

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

########################################
#### Pheno = LT8_DEL_MNe_MN ############
## Run The Script
/projects/janssen/Poly_Train/SCRIPTS/Poly_Train.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${ANNOTS} \
${PHENO_FILE} \
${PHENO_TYPE} \
${COV_FILE} \
${COVS} \
${EIG_VEC} \
${PC_COUNT} \
${START_STEP}





20150119_10_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_11_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_12_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_13_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_14_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_15_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_16_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_17_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_18_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_19_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_1_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_20_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_3_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_4_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_5_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_7_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata
20150119_9_LT8_DEL_MNe_MN_DAS_BL_MN/MOD_LT8_DEL_MNe_MN_DAS_BL_MN.Rdata





RUN_N

#################################################
## RUN ON TSCC ##################################
#################################################

## Specify inputs/names/etc...
DATE=20150119
N=20
# PHENO_NAME=ACR50_100wk_AGE_SEX
# PHENO_NAME=LT8_DEL_MNe_MN_DAS_BL_MN
PHENO_NAME=LT8_DEL_lCRP_MNe_MN_lCRP_BL_MN

## Loop through N times
for num in `seq 1 ${N}`
do
	# Change TAG Identifier
	TAG=20150119_${num}
	# Make new ".run" file
	cp ${DATE}_${PHENO_NAME}.run ${TAG}_${PHENO_NAME}.run
	# Replace placeholder with Identifier
	sed -i "s/INSERTDATEHERE/${TAG}/g" ${TAG}_${PHENO_NAME}.run
	# Submit job
	qsub ${TAG}_${PHENO_NAME}.run
done

#################################################
## IMPORT TO MAC ################################
#################################################

## Specify inputs/names/etc...
DATE=20150119
N=20
# PHENO_NAME=ACR50_100wk_AGE_SEX
PHENO_NAME=LT8_DEL_MNe_MN_DAS_BL_MN

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
