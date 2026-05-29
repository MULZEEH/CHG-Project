#Computational Human Genomics
This project born from CHG course held by Ciani and Demichelis...

## auth:
Project developed by:
YO BRO UPDATE WITH UR STUFF
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)
* Marco Pulze [@MULZEEH](https://github.com/MULZEEH)

## Problem Context:

## Objectives:

***

## Activation
remember to install conda and create the environemnt:
(environemtn found in envs/general.yml)
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


GENERAL IDEA OF THE PIPELINE:
- using spia to validate the affinity of the 2 samples:
- indexing and sorting shit
- (skipping realignment arund indels)
- bqsr -> might plot something but i doubt its informative
- removing dup
- variant calling -> find the variants (vcf)
- filter the vcf and define the intersest genes (alpha genome prediction on tons of things)
- variant annotation 
- Ancestry Analysis (ethseq)
- somatic copy number calling -> plot some copy number variation
- somatic variant calling:
    - mutation signal (sigProfiler)
    - alpha gneome
    - (if time copy number changes using segment (signal analysis on this))
- purity and ploidy estimation (maybe log2R is here)

---

# YARI PROJECT INDICATION
### Introduction
Most genomic studies in the context of oncology require the characterization of somatic events
assessed through the sequencing of bulk DNA tumor samples from patients’ cohorts (i.e., multiple
patients must be studied to obtain generalizable results and answer specific biological/clinical
questions). However, the proper analysis of single patients can provide informative data.
### Aim of the project
To detect and interpret the potential clinical relevance of somatic genomic aberrations harboured in
the genome of an oncologic patient.
### Execution
- Starting from the two provided BAM files (i.e., tumor and control DNA seq from the same individual),
- identify all somatic events (e.g., SNV, CN), and define their clinical relevance. 
- Running quality checks, manual inspections, and tumor population and sub-population characterization is advisable.
- Please, report the number of somatic variants, their type, and clinical relevance. Show in tabular
format the most relevant events, indicating the involved genes, statistical significance, and the
variables you consider most relevant. Provide visual representation of the obtained
biological/technical results (min 2 figures, max 4 figures, max 2 tables). - ( Don’t feel limited by the visualizations (and code snippets) used during lessons. Integrate your results with external cBioportal
(https://www.cbioportal.org/) and COSMIC (https://cancer.sanger.ac.uk/cosmic) databases.)
---
### Each group MUST submit the following components as part of the assignment:
- Informative Table of somatic events: Create a clear, well-organized table summarizing
the identified somatic events. Include relevant details such as gene name, variant type,
predicted impact, and any annotation information.
- Mutational signature analysis: provide a summary of results and a brief interpretation
explaining the biological or clinical relevance of the observed signatures.
- Purity and ploidy estimation: present your findings along with graphical visualizations. <span style="color:blue">This for *SURE* will have to use 1 of the 4 plot for this to actually see admixture lvl and subclonal populations </span>
- Choose one relevant somatic/genomic event and visualize it using IGV: Make sure the
visualization: Clearly supports the identified event; is appropriately zoomed and colored to
highlight the features of interest.
- Biological relevance contextualization: Write a short text explaining the biological
significance of your findings. Use COSMIC and cBioportal repositories as support.
---
Links:
[COSMIC](https://cancer.sanger.ac.uk/cosmic)
[cBioPortal](https://www.cbioportal.org/) <span style="color:blue"> this might be really interesting to compare after looking at changes in expression using *AlphaGenome* to *SELECT* the *GENES* <span>
[CBioPortal-CLI-Utilities](https://github.com/cBioPortal/cbiohubpy)
Format of the required project report:
Please comply with the following:
- 4 pages max, font size Arial 11, single line spacing.
- Include figures and related legends (min 2, max 4);
- include the following information: names of the students, project rationale (max 10 lines),
summarized computational workflow (you don’t need to report each command, but please
---
report meaningful options, such as filtering thresholds), relevant results with related
interpretations, and if necessary, pitfalls and criticisms (max 10 lines);
- figure axes, labels, and legends must be correct and complete. Size of figures should be
appropriate (visible by human eye in printed format). Figure legends should be clear,
informative, and self-contained (the reader should be able to understand the figure just by
reading the figure legend).
Notes and Suggestions:
i. You are not required to use tools other than those utilized during the classes.
ii. Provide tables and/or visualizations describing the number and features (for instance AF and
clinical relevance for mutations, log2Ratio for CNA) of identified somatic events.
iii. Scripts provided during the lessons are a good starting point. Feel free to optimize the code
for the visualization and presentation of the results.
iv. If you need to determine the genotype, you can use ASEReadCounter, then compute allelic
fraction (AF) and assign a genotype to each site considering the AF value (use the following
thresholds: AA<0.2, 0.2£AB£0.8, BB>0.8).
v. Please note that the provided BAM files are limited in genomic size for lighter processing.
Still, some intermediary files may take up some GB of space. If you use the provided virtual
machine, you may run out of space on disk: remember that some tools allow for setting a
max depth parameter to avoid wasting too much space and time during the analysis.

If you have any question, write to yari.ciani@unitn.it
Have Fun!

# MANAGMENT FOR US

### initial idea + working splits:
plots: 
- Log2/beta (SARAS)
- some track maybe/possibly (lets see)
- barplot varianti (or table) (Riccardo)
- Big figure with Tumor info (sigProfiler + Log2/beta maybe + cBioPortal with pathways/"coexpression of mutation variability"+ idk)

table:
- variants: Informative Table of somatic events: Create a clear, well-organized table summarizing the identified somatic events. Include relevant details such as gene name, variant type, predicted impact, and any annotation information. ->
In this it would be interesting to see how each of this event is present in each of the chrosomosome we are studying. (at least we will end up with more rows)
- 

text:
Spia (Mar)
EthSEQ (Ricc)
sigProfiler (Filippo)
cBioPortal (PIPUS)
AlphaGenome (Sara (+Marcol))

# work managing:
- possibly test ur script before pushing in git
- possibly lets use the notebook.Rmd file with chunks
- the initial pipeline (up to vcf (and annotation), somatic variant calling and segmentation file) on a different file (possibly Snakefile ) -> the rest on notebook.Rmd
- Ask before touching the chunk of someone
- Let's use the MarkDown Capability to describe what the chunk does (and what step of the "snkaefile pipeline" needs)
- try avoid HardCoding shit or sharing private keys to the workd

THANKS AND ENJOY GUYS :)

# note:

Clinvar found nothing but this:

[BRCA1](https://www.ncbi.nlm.nih.gov/clinvar/variation/54108/?term=54108%5BVariation+ID%5D)