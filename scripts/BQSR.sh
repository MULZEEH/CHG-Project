#!/bin/bash
#name is going to be changed for sure
# Script cna both be runned with standard or input as such
# ========== GENERAL USE (SNAKE) ===========
# sample_script.sh -o Sample.sorted.recalibrated.bam -i data/Sample.sorted.bam -r annotations/human_g1k_v37.fasta -k annotations/hapmap_3.3.b37.vcf -l annotations/CancerGenesSel.bed


REFERENCE="annotations/human_g1k_v37.fasta"
INPUT_BAM="data/Sample.sorted.bam"
KNOWN_SITES="annotations/hapmap_3.3.b37.vcf"
INTERVALS="annotations/CancerGenesSel.bed"
OUTPUT_BAM="data/Sample.sorted.recalibrated.bam"
REPORT="results/recalibration_report.pdf"

# Same but with input parameters

#  stupid way
# if [ -n $1 ]; then
#     echo "using miaomiao"
#     REFERENCE=$1
# fi
# if [ -n $2 ]; then
#     INPUT_BAM=$2
# fi
# if [ -n $3 ]; then
#     KNOWN_SITES=$3
# fi
# if [ -n $4 ]; then
#     INTERVALS=$4
# fi
# if [ -n $5 ]; then
#     OUTPUT_BAM=$5
# fi

while getopts "i:o:r:k:l:" opt; do
  case $opt in
    i)
      echo "Using input BAM: $OPTARG"
      INPUT_BAM=$OPTARG
      ;;
    o)
      OUTPUT_BAM=$OPTARG
      ;;
    r)
      echo "Using reference FASTA: $OPTARG"
      REFERENCE=$OPTARG
      ;;
    k)
      KNOWN_SITES=$OPTARG
      ;;
    l)
      INTERVALS=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# check the file presence as the parameter:
if [ ! -f "$INPUT_BAM" ]; then
    echo "Input BAM file not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Reference genome file not found: $REFERENCE"
    exit 1
fi

if [ ! -f "$INTERVALS" ]; then
    echo "Intervals BED file not found: $INTERVALS"
    exit 1
fi

if [ ! -f "$KNOWN_SITES" ]; then
    echo "Known sites VCF file not found: $KNOWN_SITES"
    exit 1
fi

if [ -f "$OUTPUT_BAM" ]; then
    echo "Output BAM file already exists: $OUTPUT_BAM"
    exit 1
fi

#========== OLDER VERSION =============
# java -jar ../../Tools/GenomeAnalysisTK.jar 
#     -T BaseRecalibrator 
#     -R ../../Annotations/human_g1k_v37.fasta 
#     -I ../../02_Realignment/Data/Sample.sorted.realigned.bam 
#     -knownSites ../../Annotations/hapmap_3.3.b37.vcf 
#     -o recal.table 
#     -L ../../Annotations/CancerGenesSel.bed

echo "Step 1: BaseRecalibrator Generate recalibration table..."
gatk BaseRecalibrator \
    -R $REFERENCE \
    -I $INPUT_BAM \
    --known-sites $KNOWN_SITES \
    -O .tmp/recal.table\
    -L $INTERVALS 

#========== OLDER VERSION =============
# java -jar ../../Tools/GenomeAnalysisTK.jar 
#     -T PrintReads 
#     -R ../../Annotations/human_g1k_v37.fasta 
#     -I ../../02_Realignment/Data/Sample.sorted.realigned.bam 
#     -BQSR recal.table 
#     -o Sample.sorted.realigned.recalibrated.bam 
#     -L ../../Annotations/CancerGenesSel.bed 
#     --emit_original_quals

echo "Step 2: Apply BQSR recalibration..."
gatk ApplyBQSR \
    -R $REFERENCE \
    -I $INPUT_BAM \
    -bqsr .tmp/recal.table \
    --emit-original-quals true \
    -O $OUTPUT_BAM 
    # -L $INTERVALS 


#========== OLDER VERSION =============
# java -jar ../../Tools/GenomeAnalysisTK.jar
#      -T BaseRecalibrator 
#      -R ../../Annotations/human_g1k_v37.fasta
#      -I ../../02_Realignment/Data/Sample.sorted.realigned.bam 
#      -knownSites ../../Annotations/hapmap_3.3.b37.vcf 
#      -BQSR recal.table 
#      -o after_recal.table 
#      -L ../../Annotations/CancerGenesSel.bed

echo "Step 3: Generate second recalibration table for validation..."
gatk BaseRecalibrator \
    -R $REFERENCE \
    -I $INPUT_BAM \
    --known-sites $KNOWN_SITES \
    -O .tmp/after_recal.table \
    # -L $INTERVALS


echo "Step 4: Create recalibration plots..."
gatk AnalyzeCovariates \
    -before .tmp/recal.table \
    -after .tmp/after_recal.table \
    -csv results/recalibration_covariates.csv \
    -plots $REPORT

echo "Step 5: Verify original qualities preserved..."
echo "Original qualities (OQ) count:"
samtools view $OUTPUT_BAM | grep OQ | wc -l

echo "Done!"
