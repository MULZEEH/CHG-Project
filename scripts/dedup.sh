
REFERENCE="annotations/human_g1k_v37.fasta"
INPUT_BAM="data/Control.sorted.recalibrated.bam"
KNOWN_SITES="annotations/hapmap_3.3.b37.vcf"
INTERVALS="data/Captured_Regions.bed"
OUTPUT_BAM="data/Control.sorted.recalibrated.dedup.bam"
METRICS="results/dedup_metrics.txt"

while getopts "i:o:r:k:l:m:" opt; do
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
    m)
      METRICS=$OPTARG
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

gatk MarkDuplicates \
        -I $INPUT_BAM \
        -O $OUTPUT_BAM \
        -M $METRICS \
        --REMOVE_DUPLICATES true \
        -AS