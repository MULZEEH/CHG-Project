from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import matplotlib as plt
import os
from dotenv import load_dotenv
import pandas as pd

# 1. Setup your credentials
API_KEY = 
client = dna_client.create(API_KEY)

def score_my_variant(chrom, pos, ref, alt, genome_build='hg19'):
    """
    Scores a single variant and returns the predicted impact.
    """
    print(f"Scoring {chrom}:{pos} {ref}>{alt}...")
    
    # The client handles the 1Mb windowing automatically
    results = client.score_variant(
        chromosome=chrom,
        position=pos,
        reference=ref,
        alternative=alt,
        genome_build=genome_build
    )
    model = dna_client.create(API_KEY)

    # as in yea
    interval = genome.Interval(chromosome='chr22', start=35677410, end=36725986)
    variant = genome.Variant(
        chromosome='chr22',
        position=36201698,
        reference_bases='A',
        alternate_bases='C',
    )

    outputs = model.predict_variant(
        interval=interval,
        variant=variant,
        ontology_terms=['UBERON:0001157'],
        requested_outputs=[dna_client.OutputType.RNA_SEQ],
    )

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
        
    return results

# 2. Example: Scoring a mutation in a known enhancer region
# Let's say we found a mutation in a regulatory region near a cancer gene
variant_impact = score_my_variant("chr17", 7577121, "G", "A")

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