#!/bin/bash

set -e
set -o pipefail

#   A simple script to hold variables for SNP Calling
#   Are you using the Common_Config file?
#       If so, where is it?
COMMON=

##############################################################################################
#   If we aren't using the Common_Config file, specify these variables
#   If Common_Config is specified, leave these blank
#   Define a list of samples
SAMPLE_LIST=

#   Define a list of inbreeding coefficients
SAMPLE_INBREEDING=

#   Ancestral and Reference sequences
ANC_SEQ=
REF_SEQ=

#   Name the project
PROJECT=

#   Where do we put the outfiles?
    #   Note, the final outdirectory will be
    #   ${SCRATCH}/${PROJECT}/SFS
SCRATCH=

#   Define the region being looked at
#       Optional, but ANGSD is expensive to run without specifying regions to look at
REGIONS=

#   Parameters that are specified in Common_Config
#       Calculate genotype likelihoods
GT_LIKELIHOOD=2
#       Set the minimum mapping quality for a base to be used
N_CORES=32
#       Determine major and minor alleles
DO_MAJORMINOR=1
#       Calculate allele frequencies
DO_MAF=1

##############################################################################################

#   SNP Calling Parameters
#       Listed below are the defaults, please modify for your samples
#       Output genotype likelihood frequency file
DO_GLF=
#       Set the maximum p-value for polymorphic sites
SNP_PVAL=1e-6
