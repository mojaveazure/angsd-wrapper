#!/usr/bin/bash

set -e
set -o pipefail

#   Load variables from supplied config file
source "$1"

#   Are we using Common_Config? If so, source it
if [[ -f "${COMMON}" ]]
then
    source "${COMMON}"
fi

#   Where is angsd-wrapper located?
SOURCE="$2"

#   Where is ANGSD located?
ANGSD_DIR="${SOURCE}/dependencies/angsd"

#   Create outdirectory
OUT="${SCRATCH}/${PROJECT}/SNP_Calls"
mkdir -p "${OUT}"

if [[ -f "${REGIONS}" ]]
then
    ${ANGSD_DIR}/angsd \
        -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -doMajorMinor "${DO_MAJOR_MINOR}" \
        -out "${OUT}/${PROJECT}" \
        -rf "${REGIONS}"
elif [[ -z "${REGIONS}" ]]
then
    ${ANGSD_DIR}/angsd \
        -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -doMajorMinor "${DO_MAJOR_MINOR}" \
        -out "${OUT}/${PROJECT}"
else
    ${ANGSD_DIR}/angsd \
        -bam "${SAMPLE_LIST}" \
        -GL "${GT_LIKELIHOOD}" \
        -doMaf "${DO_MAF}" \
        -SNP_pval "${SNP_PVAL}" \
        -doMajorMinor "${DO_MAJOR_MINOR}" \
        -out "${OUT}/${PROJECT}" \
        -r "${REGIONS}"
fi
