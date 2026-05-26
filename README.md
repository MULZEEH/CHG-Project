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
