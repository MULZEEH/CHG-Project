#Computational Human Genomics
This project born from CHG course held by Ciani and Demichelis...

## auth:
Project developed by:
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)

## Problem Context:

## Objectives:

***

## General Pipeline:
The general pipeline structured using Snakemake follows this structure:
1. Quality Control
2. Preprocessing (idk if 2. or 1.)
3. ...

### Quality Control
this step bla bla
### Preprocessing
this step bla bla

---

## Installation:

## Usage:


## Additional Notes:
for this Project we have adopted some State of Art tools such as:
- Samtools
- GATK (v.4)
- ...

############################################
# 3. Base Quality Score Recalibration (BQSR)
############################################

# Control
# Building recalibration model
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -o control.recal.table -L ../annotations/captured_regions.bed
# Creating a new bam file using the input table generated from the previous command which has accurate base substitution, insertion and deletion quality scores
java -jar ../tools/genome_analysis_TK.jar -T PrintReads -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -BQSR control.recal.table -o control.sorted.realigned.recalibrated.bam -L ../annotations/captured_regions.bed --emit_original_quals
# Second pass to evaluate what the data looks like after recalibration
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./control.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -BQSR control.recal.table -o control.after_recal.table -L ../annotations/captured_regions.bed
# Generating csv file based on before and after recalibration tables
java -jar ../tools/genome_analysis_TK.jar -T AnalyzeCovariates -R ../annotations/human_g1k_v37.fasta -before control.recal.table -after control.after_recal.table -csv control.recal.csv

# Tumor
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -o tumor.recal.table -L ../annotations/captured_regions.bed
java -jar ../tools/genome_analysis_TK.jar -T PrintReads -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -BQSR tumor.recal.table -o tumor.sorted.realigned.recalibrated.bam -L ../annotations/captured_regions.bed --emit_original_quals
java -jar ../tools/genome_analysis_TK.jar -T BaseRecalibrator -R ../annotations/human_g1k_v37.fasta -I ./tumor.sorted.realigned.bam -knownSites ../annotations/hapmap_3.3.b37.vcf -BQSR tumor.recal.table -o tumor.after_recal.table -L ../annotations/captured_regions.bed
java -jar ../tools/genome_analysis_TK.jar -T AnalyzeCovariates -R ../annotations/human_g1k_v37.fasta -before tumor.recal.table -after tumor.after_recal.table -csv tumor.recal.csv

# Executing R script to visualize the BQSR results
Rscript ../bqsr_plot_generator.R

######################################
# 4. Identifying & Removing Duplicates
######################################

# Control
# Removing duplicates
java -jar ../tools/picard.jar MarkDuplicates I=control.sorted.realigned.recalibrated.bam O=control.sorted.realigned.recalibrated.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=control.sorted.realigned.recalibrated.picard.log ASSUME_SORTED=true
# Reindexing the new deduplicated bam file
samtools index control.sorted.realigned.recalibrated.dedup.bam
# Generating flagstat report for the deduplicated bam file
samtools flagstat control.sorted.realigned.recalibrated.dedup.bam > control.sorted.realigned.recalibrated.dedup.flagstat.txt

# Tumor
java -jar ../tools/picard.jar MarkDuplicates I=tumor.sorted.realigned.recalibrated.bam O=tumor.sorted.realigned.recalibrated.dedup.bam REMOVE_DUPLICATES=true TMP_DIR=/tmp METRICS_FILE=tumor.sorted.realigned.recalibrated.picard.log ASSUME_SORTED=true
samtools index tumor.sorted.realigned.recalibrated.dedup.bam
samtools flagstat tumor.sorted.realigned.recalibrated.dedup.bam > tumor.sorted.realigned.recalibrated.dedup.flagstat.txt

################################
# 5. Somatic Copy Number Calling
################################

# Generating mpileup files for both control and tumor samples and piping them to VarScan for copy number analysis
samtools mpileup -q 1 -f ../annotations/human_g1k_v37.fasta control.sorted.realigned.recalibrated.dedup.bam tumor.sorted.realigned.recalibrated.dedup.bam | java -jar ../tools/var_scan.v2.3.9.jar copynumber --output-file SCNA --mpileup 1
# Calling somatic copy number variations
java -jar ../tools/var_scan.v2.3.9.jar copyCaller SCNA.copynumber --output-file SCNA.copynumber.called

# Executing R script for Circular Binary Segmentation (CBS) to analyze the somatic copy number variations
Rscript ../CBS.R

####################
# 6. Variant Calling
####################

# Control
# Using bcftools to call variants
bcftools mpileup -Ou -a DP -f ../annotations/human_g1k_v37.fasta control.sorted.realigned.recalibrated.dedup.bam | bcftools call -Ov -c -v > control.BCF.vcf
# Using GATK to call variants
java -jar ../tools/genome_analysis_TK.jar -T UnifiedGenotyper -R ../annotations/human_g1k_v37.fasta -I control.sorted.realigned.recalibrated.dedup.bam -o control.GATK.vcf -L ../annotations/captured_regions.bed
# Filtering the VCF files for the BCF and GATK outputs
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf control.BCF.vcf --out control.BCF --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf control.GATK.vcf --out control.GATK --recode --recode-INFO-all
# Comparing the two VCF files
vcftools --vcf control.BCF.recode.vcf --diff control.GATK.recode.vcf --diff-site --out control.BCF_vs_GATK

# Tumor
bcftools mpileup -Ou -a DP -f ../annotations/human_g1k_v37.fasta tumor.sorted.realigned.recalibrated.dedup.bam | bcftools call -Ov -c -v > tumor.BCF.vcf
java -jar ../tools/genome_analysis_TK.jar -T UnifiedGenotyper -R ../annotations/human_g1k_v37.fasta -I tumor.sorted.realigned.recalibrated.dedup.bam -o tumor.GATK.vcf -L ../annotations/captured_regions.bed
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf tumor.BCF.vcf --out tumor.BCF --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf tumor.GATK.vcf --out tumor.GATK --recode --recode-INFO-all
vcftools --vcf tumor.BCF.recode.vcf --diff tumor.GATK.recode.vcf --diff-site --out tumor.BCF_vs_GATK

#######################
# 7. Variant Annotation
#######################

# Control
# Annotating variants from bcf
java -Xmx4g -jar ../tools/snpEff/snpEff.jar -v hg19kg control.BCF.recode.vcf -s control.BCF.recode.ann.html > control.BCF.recode.ann.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/hapmap_3.3.b37.vcf control.BCF.recode.ann.vcf > control.BCF.recode.ann2.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/clinvar_pathogenic.vcf control.BCF.recode.ann2.vcf > control.BCF.recode.ann3.vcf
# Annotating variants from GATK
java -Xmx4g -jar ../tools/snpEff/snpEff.jar -v hg19kg control.GATK.recode.vcf -s control.GATK.recode.ann.html > control.GATK.recode.ann.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/hapmap_3.3.b37.vcf control.GATK.recode.ann.vcf > control.GATK.recode.ann2.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/clinvar_pathogenic.vcf control.GATK.recode.ann2.vcf > control.GATK.recode.ann3.vcf

# Tumor
java -Xmx4g -jar ../tools/snpEff/snpEff.jar -v hg19kg tumor.BCF.recode.vcf -s tumor.BCF.recode.ann.html > tumor.BCF.recode.ann.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/hapmap_3.3.b37.vcf tumor.BCF.recode.ann.vcf > tumor.BCF.recode.ann2.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/clinvar_pathogenic.vcf tumor.BCF.recode.ann2.vcf > tumor.BCF.recode.ann3.vcf
java -Xmx4g -jar ../tools/snpEff/snpEff.jar -v hg19kg tumor.GATK.recode.vcf -s tumor.GATK.recode.ann.html > tumor.GATK.recode.ann.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/hapmap_3.3.b37.vcf tumor.GATK.recode.ann.vcf > tumor.GATK.recode.ann2.vcf
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/clinvar_pathogenic.vcf tumor.GATK.recode.ann2.vcf > tumor.GATK.recode.ann3.vcf

############################
# 8. Somatic Variant Calling
############################

# Preparing control and tumor pileup files
samtools mpileup -q 1 -f ../annotations/human_g1k_v37.fasta control.sorted.realigned.recalibrated.dedup.bam > control.sorted.realigned.recalibrated.dedup.pileup
samtools mpileup -q 1 -f ../annotations/human_g1k_v37.fasta tumor.sorted.realigned.recalibrated.dedup.bam > tumor.sorted.realigned.recalibrated.dedup.pileup
# Running VarScan2 for somatic variant calling
java -jar ../tools/var_scan.v2.3.9.jar somatic control.sorted.realigned.recalibrated.dedup.pileup tumor.sorted.realigned.recalibrated.dedup.pileup --output-snp somatic.pm --output-indel somatic.indel --output-vcf 1
# Annotating the somatic point mutations
java -Xmx4g -jar ../tools/snpEff/snpSift.jar Annotate ../annotations/hapmap_3.3.b37.vcf somatic.pm.vcf > somatic.pm.vcf.hapmap_ann.vcf
# Filtering the annotated VCF file: one for SNPs only and one for non-SNPs
cat somatic.pm.vcf.hapmap_ann.vcf | java -Xmx4g -jar ../tools/snpEff/snpSift.jar filter "(exists ID) & ( ID =~ 'rs' )" > somatic.pm.onlySNPs.vcf
cat somatic.pm.vcf.hapmap_ann.vcf | java -Xmx4g -jar ../tools/snpEff/snpSift.jar filter "!(exists ID) & !( ID =~ 'rs' )" > somatic.pm.noSNPs.vcf

#################################
# 9. Purity and Ploidy Estimation
#################################

# Filtering the control VCF file for biallelic SNPs 
bcftools view -v snps -m2 -M2 control.BCF.vcf > control.BCF.biallelic_snps.vcf
# Extracting heterozygous SNPs from the biallelic SNPs VCF file
grep -E "(^#|0/1)" control.BCF.biallelic_snps.vcf > control.het.biallelic_snps.vcf
# Counting the allelic reads in the control and tumor samples
java -jar ../tools/genome_analysis_TK.jar -T ASEReadCounter -R ../annotations/human_g1k_v37.fasta -o control.csv -I control.sorted.realigned.recalibrated.dedup.bam -sites control.het.biallelic_snps.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20
java -jar ../tools/genome_analysis_TK.jar -T ASEReadCounter -R ../annotations/human_g1k_v37.fasta -o tumor.csv -I tumor.sorted.realigned.recalibrated.dedup.bam -sites control.het.biallelic_snps.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20
# Generating the somatic.pm and somatic.indel files without the .vcf extension for TPES
java -jar ../tools/var_scan.v2.3.9.jar somatic control.sorted.realigned.recalibrated.dedup.pileup tumor.sorted.realigned.recalibrated.dedup.pileup --output-snp somatic.pm --output-indel somatic.indel

# Executing the R script for purity and ploidy estimation
Rscript ../purity_and_ploidy_estimation.R