#!/bin/bash

REFERENCE="annotations/human_g1k_v37.fasta"
CONTROL_BAM="data/Control.sorted.recalibrated.dedup.bam"
TUMOR_BAM="data/Tumor.sorted.recalibrated.dedup.bam"
OUTPUT_VCF="results/Control_variants.vcf"
# PROCESSED_BAM="data/Sample.bam"

while getopts "i:t:o:r:v:" opt; do
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
    v)
      OUTPUT_VCF=$OPTARG
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

# code ASK MR YARI-CHAN: 
# Since HaplotypeCaller do ade Novo allignment that should compensate the missing MiaoShiftBaseStufff
# should i also apply it on the Tumor sample? even if snp on that is silly? maybe
# VARIANT CALLING
gatk HaplotypeCaller \
            -R $REFERENCE \
            -I $CONTROL_BAM \
            -O $OUTPUT_VCF \
            --bam-output data/Control.processed.bam

vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --vcf $CONTROL_BAM --out Control.with_indels --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf $CONTROL_BAM --out Control --recode --recode-INFO-all
$CRAZY_FILTERING=1
if [ $CRAZY_FILTERING ]; then
    gatk VariantRecalibrator \
        -R Homo_sapiens_assembly38.fasta \
        -V input.vcf.gz \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.sites.vcf.gz \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.sites.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O output.recal \
        --tranches-file output.tranches \
        --rscript-file output.plots.R

fi
gatk VariantRecalibrator \
        -R annotations/human_g1k_v37.fasta \
        -V results/Control_variants.vcf \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 annoations/hapmap_3.3.b37.vcf \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O output.recal \
        --tranches-file output.tranches \
        --rscript-file output.plots.R
# not sure
gatk HaplotypeCaller \
            -R $REFERENCE \
            -I $TUMOR_BAM \
            -O  results/Tumor_variants.vcf\
            --bam-output data/Tumor.processed.bam

vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --vcf $TUMOR_BAM --out Tumor.with_indels --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf $TUMOR_BAM --out Tumor --recode --recode-INFO-all


# VARIANT ANNOTATION
samtools mpileup -q 1 -f annotations/human_g1k_v37.fasta Control.BCF.recode.vcf > .tmp/Contorl_pileup
samtools mpileup -q 1 -f annotations/human_g1k_v37.fasta Tumor.BCF.recode.vcf   > .tmp/Tumor_pileup

snpEff -Xmx8g -v hg19kg Control.BCF.recode.vcf -s Sample.BCF.recode.ann.html > Sample.BCF.recode.ann.vcf
java -Xmx4g -jar ../../Tools/snpEff/SnpSift.jar Annotate ../../Annotations/hapmap_3.3.b37.vcf  Sample.BCF.recode.ann.vcf > Sample.BCF.recode.ann2.vcf
java -Xmx4g -jar ../../Tools/snpEff/SnpSift.jar Annotate ../../Annotations/clinvar_Pathogenic.vcf Sample.BCF.recode.ann2.vcf > Sample.BCF.recode.ann3.vcf

# snpEff -Xmx{resources.mem_gb}g \
#     -v {params.database} \
#     {input.vcf} \
#     -s {output.summary} \
#     > {output.annotated_vcf} \
#     2> {log}

# SnpSift -Xmx{resources.mem_gb}g Annotate \
#     {ANNOTATION_DIR}clinvar.vcf \
#     {output.annotated_vcf} \
#     >> {output.annotated_vcf} \
#     2>> {log}