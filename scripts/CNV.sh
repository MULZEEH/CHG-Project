#!/bin/bash

REFERENCE="annotations/human_g1k_v37.fasta"
CONTROL_BAM="data/Control.sorted.bam"
TUMOR_BAM="data/Tumor.sorted.bam"
OUTPUT_COPYCALLER="data/Control.sorted.recalibrated.copycaller.bam"
OUTPUT_COPYNUMBER="data/SCNA_C"

while getopts "i:t:o:r:c:" opt; do
  case $opt in
    i)
      echo "Using input BAM: $OPTARG"
      CONTROL_BAM=$OPTARG
      ;;
    t)
      TUMOR_BAM=$OPTARG
      ;;
    o)
      OUTPUT_COPYCALLER=$OPTARG
      ;;
    r)
      echo "Using reference FASTA: $OPTARG"
      REFERENCE=$OPTARG
      ;;
    c)
      OUTPUT_COPYNUMBER=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# check the file presence as the parameter:
if [ ! -f "$CONTROL_BAM" ]; then
    echo "Control BAM file not found: $CONTROL_BAM"
    exit 1
fi

if [ ! -f "$TUMOR_BAM" ]; then
    echo "Tumor BAM file not found: $TUMOR_BAM"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Reference genome file not found: $REFERENCE"
    exit 1
fi

if [ -f "$OUTPUT_COPYNUMBER" ]; then
    echo "Copy number file already exists: $OUTPUT_COPYNUMBER"
    exit 1
fi

if [ -f "$OUTPUT_COPYCALLER" ]; then
    echo "Output BAM file already exists: $OUTPUT_COPYCALLER"
    exit 1
fi

# samtools mpileup -q 1 -f annotations/human_g1k_v37.fasta data/Control.sorted.recalibrated.dedup.bam data/Tumor.sorted.recalibrated.dedup.bam | varscan copynumber --output-file SCNA --mpileup 1
samtools mpileup -q 1 -f $REFERENCE $CONTROL_BAM $TUMOR_BAM | varscan copynumber --mpileup 1 --output-file $OUTPUT_COPYNUMBER 
# varscan copyCaller SCNA.copynumber --output-file SCNA.copynumber.called
varscan copyCaller $OUTPUT_COPYNUMBER.copynumber --output-file $OUTPUT_COPYCALLER
# Rscript scripts/CBS.R SCNA.copynumber.called results
Rscript --vanilla --quiet --args $OUTPUT_COPYCALLER results/ < scripts/CBS.R