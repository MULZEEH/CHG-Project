# Snakemake (>=5.11) comes with a code quality checker (a so called linter), that analyzes your workflow and highlights issues that should be solved in order to follow best practices, achieve maximum readability, and reproducibility. The linter can be invoked with

# snakemake --lint

# Including results in a report

# For including results into the report, the Snakefile has to be annotated with additional information. Each output file that shall be part of the report has to be marked with the report flag, which optionally points to a caption in restructured text format and allows to define a category for grouping purposes. Moreover, a global workflow description can be defined via the report directive. Consider the following example:

# report: "report/workflow.rst"


# rule all:
#     input:
#         ["fig1.svg", "fig2.png", "testdir"]


# rule c:
#     output:
#         "test.{i}.out"
#     container:
#         "docker://continuumio/miniconda3:4.4.10"
#     conda:
#         "envs/test.yaml"
#     shell:
#         "sleep `shuf -i 1-3 -n 1`; touch {output}"


# rule a:
#     input:
#         expand("test.{i}.out", i=range(10))
#     output:
#         report("fig1.svg", caption="report/fig1.rst", category="Step 1")
#     shell:
#         "sleep `shuf -i 1-3 -n 1`; cp data/fig1.svg {output}"


# rule b:
#     input:
#         expand("{model}.{i}.out", i=range(10))
#     output:
#         report("fig2.png", caption="report/fig2.rst", category="Step 2", subcategory="{model}")
#     shell:
#         "sleep `shuf -i 1-3 -n 1`; cp data/fig2.png {output}"

# rule d:
#     output:
#         report(
#             directory("testdir"),
#             patterns=["{name}.txt"],
#             caption="report/somedata.rst",
#             category="Step 3")
#     shell:
#         "mkdir {output}; for i in 1 2 3; do echo $i > {output}/$i.txt; done"

# to avoid and log it correctly, you can redirect the output of snakemake to a file:
# &>file.txt

#============= PATHVAR =============
# pathvars:
#     data="data",
#     results="results"

# rule mytask:
#     input:
#         "<resources>/{dataset}.txt"
#     output:
#         "<results>/{dataset}.txt"
#     shell:
#         "some-tool {input} > {output}"

#============= WILDCARD STRATEGIES =============
# wildcard_constraints:
#     sample=r"[^_]+",
#     assembly=r"[^_]+",
#     sample_assembly=r"[^/]+",
#     num=r"[0-9]+",

# def get_ont_fq(wildcards):
#     if "filtlong" in wildcards.sample:
#         return "fastq-ont/" + wildcards.sample + ".fastq"
#     elif "rasusa" in wildcards.sample:
#         return "fastq-ont/" + wildcards.sample + ".fastq"
#     else:
#         return glob("fastq-ont/" + wildcards.sample + ".fastq*")


#============= DB-LOOKUP =============
# rule mytask:
#     input:
#         "path/to/{dataset}.txt"
#     output:
#         "result/{dataset}.txt"
#     params:
#         some_threshold=lookup(
#             dpath="some_tool/thresholds/{dataset}",
#             within=config,
#             default=0.1
#         )
#     shell:
#         "some-tool {input} > {output}"

# ============ OUTPUT PARAMETER FOR STRIPPING =============
# output:
#         foo="result/{dataset}.txt"
#     params:
#         prefix=subpath(output.foo, strip_suffix=".txt")
#     shell:
#         "some-tool {input} -o {params.prefix}"

#============= WRAPPERS =============
# might be intersesting to implement somehting like this
# rule samtools_sort:
#     input:
#         "mapped/{sample}.bam"
#     output:
#         "mapped/{sample}.sorted.bam"
#     params:
#         "-m 4G"
#     threads: 8
#     wrapper:
#         "0.2.0/bio/samtools/sort"

# might also include in the repor tosme igv retunrs:
#============= IGV NOTEBOOK INTEGRATION =============
# import igv_notebook
# igv_notebook.init()
# igv_browser= igv_notebook.Browser(
#     {
#         "genome": "hg19",
#         "locus": "chr22:24,376,166-24,376,456",
#         "tracks": [{
#             "name": "BAM",
#             "url": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam",
#             "indexURL": "https://s3.amazonaws.com/igv.org.demo/gstt1_sample.bam.bai",
#             "format": "bam",
#             "type": "alignment"
#         }],
#         "roi": [
#             {
#                 "name": "ROI set 1",
#                 "url": "https://s3.amazonaws.com/igv.org.test/data/roi/roi_bed_1.bed",
#                 "indexed": False,
#                 "color": "rgba(94,255,1,0.25)"
#             }
#         ]
#     }
# )
# SAMPLES = ["A", "B"]

#============= PATHVARS =============

# rule pathvars:
#     input: "<path>/{}"
#     output: "<path2>/{}"
#     run: 
## config
# pathvars:
#   path: "data"
#   path2: "results"

# # 1. The target rule
# rule all:
#     input:
#         expand("plots/{sample}_final_results.svg", sample=SAMPLES)

# # 2. A processing rule
# rule plot_results:
#     input:
#         "data/{sample}.vcf"
#     output:
#         # Tagged for the report with a category
#         plot = report("plots/{sample}_final_results.svg", 
#                       caption="report/vcf_plot.rst", 
#                       category="Variant Analysis")
#     shell:
#         "python scripts/my_plotter.py {input} {output.plot}"
# rule generate_plot:
#     input:
#         "data/results.txt"
#     output:
#         # Tagged for the report
#         fig = report("plots/summary_plot.png", 
#                      caption="report/plot_desc.rst", 
#                      category="Visualizations")
#     shell:
#         "python scripts/make_plot.py {input} {output.fig}"

# rule generate_table:
#     input:
#         "data/results.txt"
#     output:
#         # Tagged for the report - CSVs become interactive tables!
#         tbl = report("tables/metrics.csv", 
#                      caption="report/table_desc.rst", 
#                      category="Data Tables")
#     shell:
#         "python scripts/make_table.py {input} {output.tbl}"

# rule all:
#     input:
#         "plots/summary_plot.png",
#         "tables/metrics.csv"


# Snakefile
# usage: snakemake -c 1 --use-conda --report report.html
# reminders:
# - 

# Initial import statemenets
import os
import pandas as pd
# import igv_notebook
from glob import glob

from scripts.queries import get_somatic_query

config: "config.yaml"
# GLOBAL VARIABLES
KEEP_JUNK = config.get("keep_junk", False)
ANNOTATION_FOLDER = config.get("annotation_folder", "annotations/")
FULL_FASTA = config.get("full_fasta", "human_g1k_v37.fasta")
STARTING_DATA = config.get("bam_path", "data/")
JAVA_MEMORY = config.get("java_memory", 8)
SNP_REF_DB = config.get("snp_ref_db","hg19kg")

# No idea how to implement this without random errors, but idea is to ask the user/config 
# file the various path to the vcf of the convalidated snp/snv of clinical knowledge
# TESTING_VCF_LIBRARY = config.get("vcf_library", list())
# INPUT_FILE = config.get("input_bam", "data/tumor.bam")

# On error behavior
onsuccess:
    if os.path.exists(".tmp/pp.txt"):
        os.remove(".tmp/pp.txt")
    print("Workflow finished successfully")

onerror:
    if (not KEEP_JUNK):
        if os.path.exists(".tmp/pp.txt"):
            os.remove(".tmp/pp.txt")
    print("Workflow finished badly, go check log/ folder ")




# Wildcard bounds/regex to avoid unintentional hits
wildcard_constraints:
    sample=r"[^_]+",
    assembly=r"[^_]+",
    sample_assembly=r"[^/]+",
    num=r"[0-9]+", 

GROUPS=["control", "tumor"]

rule all:
    input: 
        expand("data/{group}.sorted.bam", group=GROUPS),
        expand("data/{group}.sorted.bai", group=GROUPS),
        # expand("{group}.intervals", group=GROUPS),
        # [list of things in the report]
        expand("data/{group}.sorted.recalibrated.bam", group=GROUPS)
        # ["results/*.svg", ""]

# first rule to always run, 
# to check the presence of all the correct
# Firstly check the presence of the control and tumor bam
# Currently using WRAPPER just out of curiosity
rule preprep:
    threads:
        1
    input: 
        "data/{group}.bam"
    output:
        "data/{group}.sorted.bam"
    log:
        "logs/samtools/sort/{group}.log"
    wrapper:
        "0.2.0/bio/samtools/sort"
    # run: 
    #     """
    #     echo "{wildcard.group} blabla"
    #     samtools sort {wildcard.group}
    #     echo "done"
    #     """
rule indexing:
    input: 
        "data/{group}.sorted.bam"
    output: 
        "data/{group}.sorted.bai"
    wrapper:
        "0.2.0/bio/samtools/index" 


# no idea in what execution i should use this
rule SPIA_check:
    input: 
        control="data/control.sorted.bam"
        tumor="data/tumor.sorted.bam"
    output: 
        test_result=report(
            "results/control/deossi.pdf", 
            caption="report/table_desc.rst", 
            category="Plot"
        )
    run: 
        """
        echo'Run of the SPIA script getting the quality of the SNP matching between the tumor and the contorl samples: if the distance is less then x caption will be OK, if intermediate value is likely to be a 1st degree relative, otherwise is NONO'
        """

#  gatk IndexFeatureFile -I annotations/*.vcf
def get_gatk_jar():
    return glob("tools/gatk/*.jar")[0]
# This will not be used since HaplotypeCaller of GATK4 already handles local de noveo assembly for indel realigniment
rule realignertargetcreator:
    input:
        bam="{sample}.sorted.bam",
        bai="{sample}.sorted.bai",
        ref=ANNOTATION_FOLDER + FULL_FASTA,
        fai=ANNOTATION_FOLDER + FULL_FASTA + ".fai",
        # dict="genome.dict",
        # known="dbsnp.vcf.gz",
        # known_idx="dbsnp.vcf.gz.tbi",
    output:
        intervals="{sample}.intervals",
    log:
        "logs/gatk/realignertargetcreator/{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    resources:
        mem_mb=1024,
    threads: 16
    wrapper:
        "master/bio/gatk3/realignertargetcreator"

rule Realignment_Maybe:
    input:
        bam="{sample}.sorted.bam",
        bai="{sample}.sorted.bai",
        ref=ANNOTATION_FOLDER + FULL_FASTA,
        fai=ANNOTATION_FOLDER + FULL_FASTA + ".fai",
        # dict="genome.dict",
        # known="dbsnp.vcf.gz",
        # known_idx="dbsnp.vcf.gz.tbi",
        target_intervals="{sample}.intervals",
    output:
        bam="{sample}.realigned.bam",
        bai="{sample}.realigned.bai",
    log:
        "logs/gatk3/indelrealigner/{sample}.log",
    params:
        extra="--defaultBaseQualities 20 --filter_reads_with_N_cigar",  # optional
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
        "0.2.0/bio/gatk3/indelrealigner"

# Step creating the jarjarbeans
rule BaseQualityScoreRecalibration:
    conda:
        "envs/general.yml"
    input: 
        input_bam=STARTING_DATA+"{group}.sorted.bam",
        ref=ANNOTATION_FOLDER + FULL_FASTA,
        known_sites=ANNOTATION_FOLDER + "hapmap_3.3.b37.vcf", #this can be standardized 
        intervals=ANNOTATION_FOLDER + "CancerGenesSel.bed"
    output: 
        output_bam="data/{group}.sorted.recalibrated.bam"
    script: # THIS RULE COULD RETURN ALREADY A PLOT OR TABLE OF THE INFORMAITION REGUADING THE ALLIGNMENT
        "scripts/BQSR.sh -r {input.ref} -i {input.input_bam} -o {output.output_bam} -k {input.known_sites} -l {input.intervals} "
        # I also need a better way to visualize the results

# output: 
#         output_vcf=report("results/{group}.vcf", category="Variants"),
#         bamout="results/{group}.bamout"
# Considering of using 
rule RemoveDuplicates:
    conda:
        "envs/general.yml"
    threads: 16
    input: 
        input bam=,
    output: 
        output_bam=,
        metrics=,
    run: """
    echo "Dedup step"
    gatk MarkDuplicates \
        -I {input.input_bam} \
        -O {output.output_bam} \
        -M {output.metrics} \
        -ASO=true

    samtools index {output.output_bam}
    samtools flagstat {output.output_bam} > {output.metrics}

    gatk MarkDuplicates \
        -I {input.input_bam} \
        -O {output.output_bam} \
        -M {output.metrics} \
        --REMOVING-DUPLICATES true \
        -AS
    """

rule CopyNumberVariation:
    input: 
    output: 
    run: """"
    samtools mpileup -q 1 -f {input.ref} {input.input_bam} | \
        varscan copynumber --mpileup 1 --output-file {output.copynumber} --output-file-prefix {wildcards.group}
    
    varscan copyCaller {output.copynumber} --output-file {output.copycaller} --output-file-prefix {wildcards.group}
    
    Rscript scripts/plot_cnv.R {output.copycaller} {output.plot}
    """
    
rule AncestryAnalysis:
    input: 
    output: 
    run: """
    echo "Performing Ancestry Analisys using EthSEQ"""

rule VariantCalling:
    conda:
        "envs/general.yml"
    threads: 16
    input: 
        input_bam="{RemoveDuplicates.output.output_bam}",
        ref=
    output: 
        output_vcf="results/{group}.vcf",
        bamout="results/{group}.bamout" # diventa nuovo main bam
    run: """
    gatk HaplotypeCaller  \
        -R {input.ref} \
        -I {input.input_bam} \
        -O {output.output_vcf} \
        -bamout {output.bamout}
    echo "perfroming quality again"
   """

#simple function that fiven the name of the input(contorl/tumor) reutnr the filename of the pileup (.bam -> .pileup)
def getname(wildcard):
    todo()  

rule SomaticVariantCall:
    input: 
    output: 
        control=tmp(PATH)
        tumor=tmp(PATH)
    run:
        """
        samtools mpileup -q 1 -f {input.ref} {input.control} > {output.control}
        samtools mpileup -q 1 -f {input.ref} {input.tumor} > {output.tumor}
        
        varscan -Xmx{resources.mem_gb}g somatic {output.control} {output.tumor} --output-snp somatic.pm --output-indel somatic.indel --ouput-vcf 1
        """

# Rule that handles the Alpha Genome run
rule VariantPrediction:
    input: 
        vcf="results/intermediate.vcf" # change name 
    output: 
        report=report("results/variant_prediction_report.html", category="Variant Prediction with AlphaGenome"),
    script:
        "scripts/variant_prediction.py {input} {output}"

# might incorporate this with the above
rule IGVPredictionVisualization:
    input: 
    output: 
    script:
        "script/igv_visualization.py {input} {output}"

def get_annotation_db(wildcard):
    # this is a function that given the name of the input vcf, return the correct annotation database for snpeff
    # for example if the input is hg19.vcf it will return hg19kg, if the input is hg38.vcf it will return hg38kg, etc
    # in addition it should merge the various annotation database with the custom chosen by the user (if chosen)
    todo()

# for this rule for the annotation and filtering part i was thinking of using a file with the possible query or idk how
rule VariantAnnotation:
    input:
    output:
    resources:
        mem_gb=JAVA_MEMORY
    params:
        database=SNP_REF_DB,

    run:"""
     snpEff -Xmx{resources.mem_gb}g -v {params.database} {input} -s {output.reportino} > {output.output}
     SnpSift -Xmx{resources.mem_gb}g Annotate 
     """
    

# 
rule PurityPloidy:
    input: 
        control="bcf.vcf"
    output: 
    run: 
        """
        echo "Filtering the control VCF file for biallelic SNPs"
        bcftools view \
            -v snps \
            -m2 -M2 {input.control} > {output.control}
        echo "Extracting heterozygous SNPs from the biallelic SNPs VCF file"
        grep -E "(^#|0/1)" \
            control.BCF.biallelic_snps.vcf > control.het.biallelic_snps.vcf
        echo "Counting the allelic reads in the control and tumor samples"
        java -jar ../tools/genome_analysis_TK.jar -T ASEReadCounter -R ../annotations/human_g1k_v37.fasta -o control.csv -I control.sorted.realigned.recalibrated.dedup.bam -sites control.het.biallelic_snps.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20
        echo "Counting the allelic reads in the tumor sample"
        java -jar ../tools/genome_analysis_TK.jar -T ASEReadCounter -R ../annotations/human_g1k_v37.fasta -o tumor.csv -I tumor.sorted.realigned.recalibrated.dedup.bam -sites control.het.biallelic_snps.vcf -U ALLOW_N_CIGAR_READS -minDepth 20 --minMappingQuality 20 --minBaseQuality 20
        echo "Generating the somatic.pm and somatic.indel files without the .vcf extension for TPES"
        java -jar ../tools/var_scan.v2.3.9.jar somatic control.sorted.realigned.recalibrated.dedup.pileup tumor.sorted.realigned.recalibrated.dedup.pileup --output-snp somatic.pm --output-indel somatic.indel
        """