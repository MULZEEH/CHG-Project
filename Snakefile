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

# Snakefile
# usage: snakemake -c 1 --use-conda --report report.html
# reminders:
# - 

# Initial import statemenets
import os
import pandas as pd
# import igv_notebook
from glob import glob

config: "config.yaml"
# GLOBAL VARIABLES
KEEP_JUNK = config.get("keep_junk", False)
ANNOTATION_FOLDER = config.get("annotation_folder", "annotations/")
FULL_FASTA = config.get("full_fasta", "human_g1k_v37.fasta")

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
        input_bam="data/{group}.sorted.bam",
        ref=ANNOTATION_FOLDER + FULL_FASTA,
        known_sites=ANNOTATION_FOLDER + "hapmap_3.3.b37.vcf", #this can be standardized 
        intervals=ANNOTATION_FOLDER + "CancerGenesSel.bed"
    output: 
        output_bam="data/{group}.sorted.recalibrated.bam"
    script: # THIS RULE COULD RETURN ALREADY A PLOT OR TABLE OF THE INFORMAITION REGUADING THE ALLIGNMENT
        "scripts/BQSR.sh -r {input.ref} -i {input.input_bam} -o {output.output_bam} -k {input.known_sites} -l {input.intervals} "
        # I also need a better way to visualize the results


rule RemoveDuplicates:
    input: 
        input bam=,

    output: 
        output_bam=,
        metrics=,
    run: """
    echo "Dedup step"
    gatk MarkDuplicates \
        -I {input.input_bam} \
        -O {output.ouptut} \
        -M {output.metrics}
    """

# rule VariantCalling:
#     input: 
#     output: 
#     run: 

# rule SPIA:
#     input: 
#     output: 
#     run: 