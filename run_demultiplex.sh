#!/usr/bin/bash
#SBATCH --time=1-00:00:00

set -o "errexit"

module load bbtools

FASTQ=$1
OUTPUTDIR=$2
OUTPUT=TMP.${FASTQ}

module load perl

perl ./extract_barcode_from_read.pl  -f ${FASTQ} | gzip > ${OUTPUT}

bbtools demuxbyname in=${OUTPUT} \
    out=${OUTPUTDIR}/out_%.fastq \
    outu=${OUTPUTDIR}/unmatched.fastq \
    suffixmode=t \
    length=8 \
    hdist=1 \
    names=AATAATGT,CAACACTT,ATAATTCT,GTCCATAT,CAAGTGAT,CGACTTGG,GCGAGTTG,AAGACGGG,CTGTAGGG,GCCGAGGG,GGTACCGG,GCCCTCCG,ACATGGGT,GTACCGAG,AAGTGCCG,AGTTAGAG,CGATCCAG,GTGGCTGC,ATACAGGC,ACAATCTC,GATGGCTC,GCAGTTAC,CTTGGGCC,CGCCAGAC


echo "#################################"
echo Success !!
echo "#################################"

