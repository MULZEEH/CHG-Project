from io import StringIO
from alphagenome.data import genome, gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
from alphagenome.visualization import plot_components
import numpy as np
import pandas as pd
import plotnine as gg

# import polar as pl

from tqdm import tqdm


import os
from dotenv import load_dotenv
import pandas as pd

# better implementation for a gene by gene analysis than ane xplorative one

# load config file
# load_dotenv("config.yaml")
# Load gene annotations (from GENCODE).
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# Define an extractor that fetches only MANE_select transcripts per gene.
# Mane select transcripts consists of of one curated transcript per locus.
gtf_transcript = gene_annotation.filter_protein_coding(gtf)
gtf_transcript = gene_annotation.filter_to_mane_select_transcript(
    gtf_transcript
)
transcript_extractor = transcript_utils.TranscriptExtractor(
    gtf_transcript
)
print(f"Loaded {len(gtf_transcript)} MANE select transcripts from GENCODE.")   

print(gtf_transcript.head())
print(head(gtf_transcript))
exit(1)
#====================================================================
#================================= Load VCF file ====================
#====================================================================
# vcf = pl.read_csv("variants.vcf", comment="#", sep="\t", header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
vcf = pd.read_csv("data/testM.vcf", comment="#", sep="\t", header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])

required_columns = ['CHROM', 'POS', 'REF', 'ALT']
for column in required_columns:
  if column not in vcf.columns:
    raise ValueError(f'VCF file is missing required column: {column}.')

# rewrite the ID in a more elegant way, and correcting the crhomosme notation
vcf['variant_id'] = vcf['CHROM'].astype(str) + "_" + vcf['POS'].astype(str) + "_" + vcf['REF'].astype(str) + ">" + vcf['ALT'].astype(str).astype(str)
vcf['CHROM'] = "chr" + vcf['CHROM'].astype(str) 

print(f"Loaded {len(vcf)} variants from VCF file.")
print(vcf.head())
#====================================================================
#============================ Init param and settings ===============
#====================================================================
# some @param
organism = 'human'
sequence_length = '1MB'  # @param ["16KB", "100KB", "500KB", "1MB"] { type:"string" }
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
    f'SEQUENCE_LENGTH_{sequence_length}'
]

# Might want to use a config file for this
# @markdown Specify which scorers to use to score your variants:
score_rna_seq = os.environ.get('ag_rna_seq', 'True').lower() == 'true'  # @param { type: "boolean"}
score_cage = os.environ.get('ag_cage', 'True').lower() == 'true'  # @param { type: "boolean" }
score_procap = os.environ.get('ag_procap', 'True').lower() == 'true'  # @param { type: "boolean" }
score_atac = os.environ.get('ag_atac', 'True').lower() == 'true'  # @param { type: "boolean" }
score_dnase = os.environ.get('ag_dnase', 'True').lower() == 'true'  # @param { type: "boolean" }
score_chip_histone = os.environ.get('ag_chip_histone', 'True').lower() == 'true'  # @param { type: "boolean" }
score_chip_tf = os.environ.get('ag_chip_tf', 'True').lower() == 'true'  # @param { type: "boolean" }
score_polyadenylation = os.environ.get('ag_polyadenylation', 'True').lower() == 'true'  # @param { type: "boolean" }
score_splice_sites = os.environ.get('ag_splice_sites', 'True').lower() == 'true'  # @param { type: "boolean" }
score_splice_site_usage = os.environ.get('ag_splice_site_usage', 'True').lower() == 'true'  # @param { type: "boolean" }
score_splice_junctions = os.environ.get('ag_splice_junctions', 'True').lower() == 'true'  # @param { type: "boolean" }

# @markdown Other settings:
download_predictions = os.environ.get('ag_download_predictions', 'False').lower() == 'true'  # @param { type: "boolean" }

# Parse organism specification. -> clearly human
organism_map = {
    'human': dna_client.Organism.HOMO_SAPIENS,
    'mouse': dna_client.Organism.MUS_MUSCULUS,
}
organism = organism_map[organism]


# Parse scorer specification. (praticamente fo tutto DIOCA)
scorer_selections = {
    'rna_seq': score_rna_seq,
    'cage': score_cage,
    'procap': score_procap,
    'atac': score_atac,
    'dnase': score_dnase,
    'chip_histone': score_chip_histone,
    'chip_tf': score_chip_tf,
    'polyadenylation': score_polyadenylation,
    'splice_sites': score_splice_sites,
    'splice_site_usage': score_splice_site_usage,
    'splice_junctions': score_splice_junctions,
}

all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]

# Remove any scorers or output types that are not supported for the chosen organism.
unsupported_scorers = [
    scorer
    for scorer in selected_scorers
    if (
        organism.value
        not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
    )
    | (
        (scorer.requested_output == dna_client.OutputType.PROCAP)
        & (organism == dna_client.Organism.MUS_MUSCULUS)
    )
]
if len(unsupported_scorers) > 0:
  print(
      f'Excluding {unsupported_scorers} scorers as they are not supported for'
      f' {organism}.'
  )
  for unsupported_scorer in unsupported_scorers:
    selected_scorers.remove(unsupported_scorer)


# Score variants in the VCF file.
#=== MIAO ====
#LOAD API KEY
load_dotenv()

api_key = os.getenv("ALPHA_GENOME_API_KEY")
print(api_key)
if not api_key:
    raise ValueError("API Key not found! Did you set it in the .env file dummy?")
dna_model = dna_client.create(str(api_key))
results = []

# to wrapp(trycatch)
for i, vcf_row in tqdm(vcf.iterrows(), total=len(vcf)):
  variant = genome.Variant(
      chromosome=str(vcf_row.CHROM),
      position=int(vcf_row.POS),
      reference_bases=vcf_row.REF,
      alternate_bases=vcf_row.ALT,
      name=vcf_row.variant_id,
  )
  interval = variant.reference_interval.resize(sequence_length)

  variant_scores = dna_model.score_variant(
      interval=interval,
      variant=variant,
      variant_scorers=selected_scorers,
      organism=organism,
  )
  results.append(variant_scores)

df_scores = variant_scorers.tidy_scores(results)

if download_predictions:
  df_scores.to_csv('variant_scores.csv', index=False)
# df_scores.to_csv('variant_scores.csv', index=False)

df_scores

exit(0)

# Currently hardcoding chromosome sizes for hg38, but in a real implementation, you would want to load this from a file or API
HG38_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
    'chrM': 16569
}

# can take the whole VCF file as input:
# No idea how to plot all that 

# not usefull for an exploration analysis but for a more systematic one, we might want to implement a function that takes as input a list of variants and outputs a dataframe with the scores for each variant and modality, which we can then filter based on the quantiles to identify significant hits.

# print("Main: Program finished.")

load_dotenv()

# # 1. Setup your credentials
api_key = os.getenv("ALPHA_GENOME_API_KEY")
print(api_key)
# if not api_key:
#     raise ValueError("API Key not found! Did you set it in the .env file?")
model = dna_client.create(str(api_key))

def score_my_variant(chrom, pos, ref, alt, genome_build='hg19'):
    """
    Scores a single variant and returns the predicted impact.
    """
    print(f"Scoring {chrom}:{pos} {ref}>{alt}...")
    
    # The client handles the 1Mb windowing automatically
    # results = client.score_variant(
    #     chromosome=chrom,
    #     position=pos,
    #     reference=ref,
    #     alternative=alt,
    #     genome_build=genome_build
    # )
    # model = dna_client.create(API_KEY)

    variant_string = f"{chrom}:{pos} {ref}>{alt}"
    # variant = genome.Variant.from_str(variant_string)
    # variant
    # as in yea
    #try catch code for the query retrieve
    try:
        # should handle automatically the size of the window but for later implementatin lets keep like this
        start_pos = max(0, pos - 524288)  # 500kb upstream (number used in the documentation)
        end_pos = min(HG38_SIZES.get(chrom, float('inf')), pos + 524288)  # 500kb downstream (still number used in the documentation)

        # might want to move these outside the try catch statement
        if start_pos >= end_pos:
            raise ValueError("Start position must be less than end position.")
        if start_pos < 0 or end_pos < 0:
            raise ValueError("Start and end positions must be non-negative.")
        
        interval = genome.Interval(chromosome=chrom, start=start_pos, end=end_pos)
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
        )

        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=['UBERON:0001157'],
            requested_outputs=[dna_client.OutputType.RNA_SEQ],
        )
        # model = dna_client.create(API_KEY)
    except Exception as e:
        print(f"Error scoring variant: {e}, might be passing wrong chromosome or chromosome build")
        return None
    
    plot_components.plot(
        [
            plot_components.OverlaidTracks(
                tdata={
                    'REF': outputs.reference.rna_seq,
                    'ALT': outputs.alternate.rna_seq,
                },
                colors={'REF': 'dimgrey', 'ALT': 'red'},
            ),
        ],
        interval=outputs.reference.rna_seq.interval.resize(2**15),
        # Annotate the location of the variant as a vertical line.
        annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
    )
    plt.show()
        
    return outputs

# 2. Example: Scoring a mutation in a known enhancer region
# Let's say we found a mutation in a regulatory region near a cancer gene
variant_impact = score_my_variant("chr22", 36201698, "A", "C")

# 3. Analyze the output
# The results usually come back as a dictionary of 'modalities' 
# (e.g., RNA_SEQ, ATAC_SEQ, CHIP_SEQ)
for modality, data in variant_impact.items():
    # We look for the 'quantile' to see if the impact is in the top 1%
    max_impact = data['quantile'].max()
    min_impact = data['quantile'].min()
    
    if max_impact > 0.99 or min_impact < 0.01:
        print(f"Significant hit found in {modality}!")
        print(f"Peak Quantile: {max_impact:.4f}")

# 4. Save to CSV for your report
# Converting the RNA_SEQ predictions to a dataframe
df_rna = pd.DataFrame(variant_impact['RNA_SEQ'])
df_rna.to_csv("alphagenome_rna_impact.csv", index=False)