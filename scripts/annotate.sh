#!/bin/bash

REFERENCE="annotations/human_g1k_v37.fasta"
CONTROL_VCF="results/Control_variants.vcf"
TUMOR_VCF="results/Tumor_variants.vcf"

# CONTROL SAMPLE VARIANT
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --vcf $CONTROL_VCF --out variant_annotations/Control.with_indels --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf $CONTROL_VCF --out variant_annotations/Control --recode --recode-INFO-all

# Diff between with or without indels
vcftools --vcf variant_annotations/Control.with_indels.recode.vcf --diff variant_annotations/Control.recode.vcf --diff-site --out variant_annotations/Control.comparison.with_indels

# Annotate with snpEff the effect of that variants
snpEff -Xmx8g -v hg19kg variant_annotations/Control.recode.vcf -s variant_annotations/Control.ann.html > variant_annotations/Control.ann.vcf
snpEff -Xmx8g -v hg19kg variant_annotations/Control.with_indels.recode.vcf -s variant_annotations/Control.with_indels.ann.html > variant_annotations/Control.with_indels.ann.vcf

mkdir -p variant_annotations/anno

# Annotate with SnpSift the presence of that variants in the HapMap and ClinVar databases (Hapmap, ClinVar, Clinvar + Hapmap)
SnpSift -Xmx8g Annotate annotations/hapmap_3.3.b37.vcf variant_annotations/Control.ann.vcf > variant_annotations/anno/Control.hapmap.ann.vcf
SnpSift -Xmx8g Annotate annotations/clinvar_Pathogenic.vcf variant_annotations/Control.hapmap.ann.vcf > variant_annotations/anno/Control.hapmap_clinvar.ann.vcf
SnpSift -Xmx8g Annotate annotations/clinvar_Pathogenic.vcf variant_annotations/Control.ann.vcf > variant_annotations/anno/Control.clinvar.ann.vcf

# Filter out the variants
cat variant_annotations/anno/Control.hapmap.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" > variant_annotations/anno/Control.filtered.hapmap.ann.vc
cat variant_annotations/anno/Control.hapmap_clinvar.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" > variant_annotations/anno/Control.filtered.hapmap_clinvar.ann.vc
cat variant_annotations/anno/Control.clinvar.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists CLNSIG)" >variant_annotations/anno/Control.filtered.clinvar.ann.vc

# TUMOR SAMPLE VARIANT

vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --vcf $TUMOR_VCF --out variant_annotations/Tumor.with_indels --recode --recode-INFO-all
vcftools --minQ 20 --max-meanDP 200 --min-meanDP 5 --remove-indels --vcf $TUMOR_VCF --out variant_annotations/Tumor --recode --recode-INFO-all

# Diff between with or without indels
vcftools --vcf variant_annotations/Tumor.with_indels.recode.vcf --diff variant_annotations/Tumor.recode.vcf --diff-site --out variant_annotations/Tumor.comparison.with_indels


snpEff -Xmx8g -v hg19kg variant_annotations/Tumor.recode.vcf -s variant_annotations/Tumor.ann.html > variant_annotations/Tumor.ann.vcf
snpEff -Xmx8g -v hg19kg variant_annotations/Tumor.with_indels.recode.vcf -s variant_annotations/Tumor.with_indels.ann.html > variant_annotations/Tumor.with_indels.ann.vcf

SnpSift -Xmx8g Annotate annotations/hapmap_3.3.b37.vcf variant_annotations/Tumor.ann.vcf > variant_annotations/anno/Tumor.hapmap.ann.vcf
SnpSift -Xmx8g Annotate annotations/clinvar_Pathogenic.vcf variant_annotations/Tumor.hapmap.ann.vcf > variant_annotations/anno/Tumor.hapmap_clinvar.ann.vcf
SnpSift -Xmx8g Annotate annotations/clinvar_Pathogenic.vcf variant_annotations/Tumor.ann.vcf > variant_annotations/anno/Tumor.clinvar.ann.vcf

cat variant_annotations/anno/Tumor.hapmap.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" > variant_annotations/anno/Tumor.filtered.hapmap.ann.vc
cat variant_annotations/anno/Tumor.hapmap_clinvar.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists ID)" > variant_annotations/anno/Tumor.filtered.hapmap_clinvar.ann.vc
cat variant_annotations/anno/Tumor.clinvar.ann.vcf | SnpSift -Xmx4g filter "(ANN[ANY].IMPACT = 'HIGH') & (DP > 20) & (exists CLNSIG)" >variant_annotations/anno/Tumor.filtered.clinvar.ann.vc


if [ $REMOVE_LOG ]; then
    rm -rf variant_annotations/*.log
fi

if [ $MOVE_GROUP ]; then
    mkdir -p variant_annotations/Control
    mkdir -p variant_annotations/Tumor

    mv variant_annotations/Control* variant_annotations/Control/
    mv variant_annotations/Tumor* variant_annotations/Tumor/
fi