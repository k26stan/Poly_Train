## Write Script to Train/Test Models using GWAS ##
## November 21, 2014 ##
## Kristopher Standish ##

######################################################
## OVERVIEW/GOALS ####################################
######################################################

## Steps ##

# Randomly Sample Patients for Training/Test Sets
  # 70/30 -> Training/Test
# Run GWAS on Training Set
# Pull out all variants beyond a given threshold
  # p < 1e-5
# Collect Beta Values and Build Regression Model
# Throw Test Set into Regression Model
  # Calculate R2 & Adjusted R2
    # R2 = 1 - (  sum( yobs - yexp )^2 / sum( yobs - ymn )^2 )
    # adjR2 = 1 - (1-R2)(n-1)/(n-p)
      # n = number of patients (in Test set)
      # p = number of parameters
  # *Could try running with subset of SNPs and find optimal model...
    # May take a while
# Run this way a bunch of times
  # Compile best predictors
  # Determine optimal beta values (somehow?)
  # Calculate R2/adjR2 for entire cohort using compiled model/predictors

##########################################################################
## 1 ## Set up Paths #####################################################
##########################################################################
 # Use Bash
 # Take in arguments, set up directories/paths for files/tools
echo \### 1 - `date` \###
echo \### Define Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

## Names/Paths for Output
DATE=$1
HOME_DIR=$2

## Parameters and Files
VAR_FILE=$3
ANNOTS=$4 # Path to Annotation File
PHENO_FILE=$5 # Which Phenotype File are you using?
PHENO_TYPE=$6 # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=$7 # Path to Covariate File or "F"
COVS=$8 # Which Covariates to Include?
EIG_VEC=$9 # Output from Plink's --pca command (MAF>1%) or "F"
PC_COUNT=${10} # How many PCs to Include as Covariates?
START_STEP=${11} # Which Step do you want to start on?

###########################################################
## Constant Paths ##

## Public Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt

## Custom Scripts
POLY_DIR=/projects/janssen/Poly_Train/SCRIPTS
s2_SAMPLE_PHENO_R=${POLY_DIR}/2-Sample_Pheno.R
s5_MAKE_COV_TAB_R=${POLY_DIR}/5-Make_Cov_Table.R
s8_MANHAT_PLOT_R=${POLY_DIR}/8-Manhat_Plot.R
s10_GET_PLINK_DAT_PY=${POLY_DIR}/10-Get_Plink_Dat.py
s10_COMPILE_OUTS_R=${POLY_DIR}/10-Compile_Outs.R
s12_Model_C_R=${POLY_DIR}/12-Model_C.R
s12_Model_B_R=${POLY_DIR}/12-Model_B.R

###########################################################
## Pull some Info out of Parameters ##

## Get Names of Specific Files
DIRS=(${VAR_FILE//\// })
VAR_FILE_NAME=${DIRS[${#DIRS[@]} - 1]} # Get Name of Variant File
DIRS=(${PHENO_FILE//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

## Specify list of Covariates to include (for command and for filename)
if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

## Specify extensions for Continuous vs Binary Phenotype
if [ $PHENO_TYPE = "C" ]
then
SUFFIX=linear
else
SUFFIX=logistic
fi

## Set up Directory for today's adventure
echo \### Moving to Base Directory at `date`: \###
echo $ASSOC
ASSOC=${HOME_DIR}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME}
mkdir ${ASSOC}

## Specify Output File Paths & Names
PHENO_SETS=${ASSOC}/${PHENO%%.txt}
NEW_COV_FILE=${ASSOC}/Cov_w_PCs.txt
ASSOC_TR=${ASSOC}/${DATE}_TR_${PHENO%%.txt}_${COVS_FILENAME}
ASSOC_TS=${ASSOC}/${DATE}_TS_${PHENO%%.txt}_${COVS_FILENAME}
P_FILE=${ASSOC}/${DATE}_${PHENO%%.txt}_${COVS_FILENAME}.P
CND_FILE=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.txt
CND_ANNOTS=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.Annot.txt
CND_GENES=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}.Gene.txt
CND_012=${ASSOC}/CND_${PHENO%%.txt}_${COVS_FILENAME}_012
cd ${ASSOC}

## Specify a File to which to Write Updates
UPDATE_FILE=${ASSOC}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 ## Sample Training/Test Set #########################################
##########################################################################
 # Use R
 # Split Phenotype File into Training/Test Sets
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Sample Training/Test Sets \###
echo `date` "2 - Sample Training/Test Sets" >> ${UPDATE_FILE}

## Take in Phenotype File
 # Spit out 2 separate Phenotype Files (Training/Test Sets)
Rscript ${s2_SAMPLE_PHENO_R} ${PHENO_FILE} ${PHENO_SETS}

## Done
echo `date` "2 - Sample Training/Test Sets - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Determine Variant Format and Adjust ##############################
##########################################################################
 # Use Bash
 # Convert Variant File into BED format
   # If File is VCF, make a PED/MAP file, then BED file
   # If File is PED, make BED file
   # If File is BED, skip ahead
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Determine/Adjust Variant File Formats \###
echo `date` "3 - Determine/Adjust Variant File Formats" >> ${UPDATE_FILE}

## Determing File Type and Converty to .bed File
if [ ${VAR_FILE: -4} == ".vcf" ] ; then
        ${VCF_TOOLS} --plink --vcf ${VAR_FILE} --out ${ASSOC}_${VAR_FILE_NAME%%.vcf}
        VAR_FILE=${ASSOC}_${VAR_FILE_NAME%%.vcf}.ped
fi
if [ ${VAR_FILE: -4} == ".ped" ] ; then
        ${PLINK} --make-bed --file ${VAR_FILE%%.ped} --out ${VAR_FILE%%.ped}
        VAR_FILE=${VAR_FILE%%.ped}.bed
fi

## Done
echo `date` "3 - Determine/Adjust Variant File Formats - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 4 ## Principal Components #############################################
##########################################################################
 # Use Bash, Plink
 # If PCs are provided, move on
 # If not, Calculate them using Plink
if [ "$START_STEP" -le 4 ]; then
echo \### 4 - `date` \###
echo \### Calculate/Specify PCs \###
echo `date` "4 - Calculate/Specify PCs" >> ${UPDATE_FILE}

## If No Principal Components Exist, Make Them
if [ $EIG_VEC = "F" ] ; then
	# Use BED file to run PCA
	${PLINK} --bfile ${VAR_FILE%%.bed} \
	--pca header \
	--allow-no-sex \
	--out ${ASSOC}/${VAR_FILE%%.bed}
	EIG_VEC=${VAR_FILE%%.bed}.eigenvec
fi

## Done
echo `date` "4 - Calculate/Specify PCs - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 5 ## Make New Covariate File ##########################################
##########################################################################
 # Use R
 # Take in Covariate File & PCs
 # Spit out Single Covariate file containing PCs
if [ "$START_STEP" -le 5 ]; then
echo \### 5 - `date` \###
echo \### Create New Covariate File \###
echo `date` "5 - Create New Covariate File" >> ${UPDATE_FILE}

## If There is no Covariate File
if [ $COV_FILE = "F" ] ; then
	# Make PC file the Covariate File
	cp ${EIG_VEC} ${NEW_COV_FILE}
else
	# Make new Covariate File
	Rscript ${s5_MAKE_COV_TAB_R} ${EIG_VEC} ${COV_FILE} ${NEW_COV_FILE}
fi

## Done
echo `date` "5 - Create New Covariate File - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 6 ## Perform Single-Locus Association on Training Set #################
##########################################################################
 # Use Plink
 # Use Training Phenotype File
if [ "$START_STEP" -le 6 ]; then
echo \### 6 - `date` \###
echo \### Perform GWAS on Training Set \###
echo `date` "6 - Perform GWAS on Training Set" >> ${UPDATE_FILE}

## Perform Association on Training Set
${PLINK} --bfile ${VAR_FILE%%.bed} \
--pheno ${PHENO_SETS}_Train.txt \
--covar ${NEW_COV_FILE} --covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar --adjust qq-plot \
--allow-no-sex \
--hardy midp \
--maf 0.01 \
--keep-allele-order \
--freqx \
--out ${ASSOC_TR}

## Done
echo `date` "6 - Create New Covariate File - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 7 ## Pull out Results (for QQ/Manhat Plots) ###########################
##########################################################################
 # Use Bash
 # Pull out only pertinent columns (for now)
if [ "$START_STEP" -le 7 ]; then
echo \### 7 - `date` \###
echo \### Pull out Results for QQ/Manhat Plots \###
echo `date` "7 - Pull out Results for QQ/Manhat Plots" >> ${UPDATE_FILE}

## Pull out just the p-values & identifying information
cat ${ASSOC_TR}.assoc.${SUFFIX} | awk '{print $1"\t"$2"\t"$3"\t"$9}' > ${P_FILE}

## Done
echo `date` "7 - Pull out GWAS Results on Training Set - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 8 ## Manhattan/QQ Plots ###############################################
##########################################################################
 # Use R
 # Take in P-values and Positions
 # Spit out Manhattand and QQ Plots
 # Spit out list of SNPs beyond a given threshold
if [ "$START_STEP" -le 8 ]; then
echo \### 8 - `date` \###
echo \### Make Manhattan and QQ Plots \###
echo `date` "8 - Make Manhattan and QQ Plots" >> ${UPDATE_FILE}

## Make Manhattan/QQ plots & Candidate Table of SNPs beyond Threshold
Rscript ${s8_MANHAT_PLOT_R} ${P_FILE} ${PHENO} ${COVS_FILENAME}

## Done
echo `date` "8 - Make Manhattan and QQ Plots - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 9 ## Pull out Cypher Annotations ######################################
##########################################################################
 # Use Tabix
 # Take in Cypher Annots and List of Candidate SNPs
 # Spit out Annotations for all Candidate SNPs
if [ "$START_STEP" -le 9 ]; then
echo \### 9 - `date` \###
echo \### Pull Cypher Annotations \###
echo `date` "9 - Pull Cypher Annotations" >> ${UPDATE_FILE}

## Pull out Variants from Candidate SNP list into CHR:START-STOP format
tail -n +2 ${CND_FILE} | awk '{print "chr"$1":"$3"-"$3}' > ${CND_FILE%%txt}list

## Use Tabix to Pull out Cypher Annotations
zcat ${ANNOTS} | head -1 > ${CND_ANNOTS}
for i in `cat ${CND_FILE%%txt}list`
do
tabix ${ANNOTS} ${i} >> ${CND_ANNOTS}
done

## Done
echo `date` "9 - Pull Cypher Annotations - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 10 ## Compile Stats and Annotations ###################################
##########################################################################
 # Use Bash & Python & R
 # Bash: Pull out Important Annotations from Cypher file
 # Python: Pull Info form Multiple Files
 # R: Compile all the Data/Annotations into Tables/Plots
if [ "$START_STEP" -le 10 ]; then
echo \### 10 - `date` \###
echo \### Compile Stats and Annotations \###
echo `date` "10 - Compile Stats and Annotations" >> ${UPDATE_FILE}

## Pull out Subset of Cypher Annotations
cat ${CND_ANNOTS} | cut -d$'\t' -f2-7,20-22,24,46-50,54,64,68,70,74,78-80,100-101,105  > ${CND_GENES}

## Use Python to Pull Info from Multiple Files
python ${s10_GET_PLINK_DAT_PY} ${CND_FILE} \
${ASSOC_TR}.hwe ${ASSOC_TR}.assoc.${SUFFIX}.adjusted ${ASSOC_TR}.frqx ${ASSOC_TR}.assoc.${SUFFIX} \
${CND_FILE%%txt}hwe ${CND_FILE%%txt}adj ${CND_FILE%%txt}frqx ${CND_FILE%%txt}pv
echo \### Done with Python, move to R \###

## Use Rscript to Pull Together Desired Info
Rscript ${s10_COMPILE_OUTS_R} ${CND_FILE} ${CND_FILE%%txt}hwe ${CND_FILE%%txt}adj ${CND_FILE%%txt}frqx ${CND_FILE%%txt}pv ${CND_GENES} ${CND_ANNOTS} ${COVS_COMMAND} ${SUFFIX}

## Done
echo `date` "10 - Compile Stats and Annotations - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 11 ## Pull Genotypes for Candidate SNPs ###############################
##########################################################################
 # Use Bash & Plink
 # Pull out SNP id for Candidate SNPs
 # Take in Candidate SNP list
 # Make PED/012 file for Candidates
if [ "$START_STEP" -le 11 ]; then
echo \### 11 - `date` \###
echo \### Pull out Genotypes for Candidates \###
echo `date` "11 - Pull out Genotypes for Candidates" >> ${UPDATE_FILE}

## Pull out ID for Candidates
cat ${CND_FILE%%txt}compile.short | awk '{print $7}' > ${CND_FILE%%.txt}_ID.list

## Output to 012 Format
${PLINK} \
--bfile ${VAR_FILE%%.bed} \
--silent \
--recode A \
--extract ${CND_FILE%%.txt}_ID.list \
--out ${CND_012} # Save as ..._012.raw

## Done
echo `date` "11 - Pull out Genotypes for Candidates - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 12 ## Filter Candidate SNPs and Build Model ###########################
##########################################################################
 # Use R
 # Take in Compiled Candidate SNP file
 # Take in Genotypes for Candidate SNPs
 # Filter out Best SNPs (good hwe, higher MAF)
 # Spit out Model w/ relevant Coefficients
if [ "$START_STEP" -le 12 ]; then
echo \### 12 - `date` \###
echo \### Filter Candidates and Build Model \###
echo `date` "12 - Filter Candidates & Build Model" >> ${UPDATE_FILE}

if [ $PHENO_TYPE = "C" ]
then
Rscript ${s12_Model_C_R} ${PHENO_SETS} ${NEW_COV_FILE} ${CND_FILE%%txt}compile.short ${CND_012}.raw ${COVS_COMMAND}
else
Rscript ${s12_Model_B_R} ${PHENO_SETS} ${NEW_COV_FILE} ${CND_FILE%%txt}compile.short ${CND_012}.raw ${COVS_COMMAND}
fi

## Done
echo `date` "12 - Filter Candidates and Build Model - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
#####################################################################
## END OF DOC #######################################################
#####################################################################