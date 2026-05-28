from io import StringIO
from alphagenome.data import genome, gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
from alphagenome.visualization import plot_components
import numpy as np
import pandas as pd
import plotnine as gg

# @title Utilities
def generate_background_variants(
    variant: genome.Variant, max_number: int = 100
) -> pd.DataFrame:
  """Generates a dataframe of background variants for a given variant.

  This is done by creating new sequences of the same length as the alternate
  allele.

  This allows us to test if the specific sequence of the oncogenic variant has a
  greater effect than a random sequence of the same length at the same location.

  Args:
    variant: The variant to generate ism variants for.
    max_number: The maximum number of ism variants to generate.

  Returns:
    A dataframe of variants.
  """
  nucleotides = np.array(list('ACGT'), dtype='<U1')

  def generate_unique_strings(n, max_number, random_seed=42):
    """Generates unique random strings of length n."""
    rng = np.random.default_rng(random_seed)

    if 4**n < max_number:
      raise ValueError(
          'Cannot generate that many unique strings for the given length.'
      )

    generated_strings = set()
    while len(generated_strings) < max_number:
      indices = rng.integers(0, 4, size=n)
      new_string = ''.join(nucleotides[indices])
      if new_string != variant.alternate_bases:
        generated_strings.add(new_string)
    return list(generated_strings)

  permutations = []
  if 4 ** len(variant.alternate_bases) < max_number:
    # Get all
    for p in itertools.product(
        nucleotides, repeat=len(variant.alternate_bases)
    ):
      permutations.append(''.join(p))
  else:
    # Sample some
    permutations = generate_unique_strings(
        len(variant.alternate_bases), max_number
    )
  ism_candidates = pd.DataFrame({
      'ID': ['mut_' + str(variant.position) + '_' + x for x in permutations],
      'CHROM': variant.chromosome,
      'POS': variant.position,
      'REF': variant.reference_bases,
      'ALT': permutations,
      'output': 0.0,
      'original_variant': variant.name,
  })
  return ism_candidates


def oncogenic_and_background_variants(
    input_sequence_length: int, number_of_background_variants: int = 20
) -> pd.DataFrame:
  """Generates a dataframe of all variants for this evaluation."""
  oncogenic_variants = oncogenic_tal1_variants()

  variants = []
  for vcf_row in oncogenic_variants.itertuples():
    variants.append(
        genome.Variant(
            chromosome=str(vcf_row.CHROM),
            position=int(vcf_row.POS),
            reference_bases=vcf_row.REF,
            alternate_bases=vcf_row.ALT,
            name=vcf_row.ID,
        )
    )

  background_variants = pd.concat([
      generate_background_variants(variant, number_of_background_variants)
      for variant in variants
  ])
  all_variants = pd.concat([oncogenic_variants, background_variants])
  return inference_df(all_variants, input_sequence_length=input_sequence_length)


def vcf_row_to_variant(vcf_row: pd.Series) -> genome.Variant:
  """Parse a row of a vcf df into a genome.Variant."""
  variant = genome.Variant(
      chromosome=str(vcf_row.CHROM),
      position=int(vcf_row.POS),
      reference_bases=vcf_row.REF,
      alternate_bases=vcf_row.ALT,
      name=vcf_row.ID,
  )
  return variant


def inference_df(
    qtl_df: pd.DataFrame,
    input_sequence_length: int,
) -> pd.DataFrame:
  """Returns a pd.DataFrame with variants and intervals ready for inference."""
  df = []
  for _, row in qtl_df.iterrows():
    variant = vcf_row_to_variant(row)

    interval = genome.Interval(
        chromosome=row['CHROM'], start=row['POS'], end=row['POS']
    ).resize(input_sequence_length)

    df.append({
        'interval': interval,
        'variant': variant,
        'output': row['output'],
        'variant_id': row['ID'],
        'POS': row['POS'],
        'REF': row['REF'],
        'ALT': row['ALT'],
        'CHROM': row['CHROM'],
    })
  return pd.DataFrame(df)


def coarse_grained_mute_groups(eval_df):
  grp = []
  for row in eval_df.itertuples():
    if row.POS >= 47239290:  # MUTE site.
      if row.ALT_len > 4:
        grp.append('MUTE' + '_other')
      else:
        grp.append('MUTE' + '_' + str(row.ALT_len))
    else:
      grp.append(str(row.POS) + '_' + str(row.ALT_len))

  grp = pd.Series(grp)
  return pd.Categorical(grp, categories=sorted(grp.unique()), ordered=True)


def find_genes_in_interval(
    bed_file_path, chromosome, position, size_kb, output_file=None
):
    """Finds genes from a BED file that overlap with a specified genomic interval.

    Parameters:
    - bed_file_path (str): Path to the BED file (tab-separated, no header).
    - chromosome (str): e.g., 'chr1' or '1' (must match BED file format).
    - position (int): The center or start position in base pairs.
    - size_kb (float): The size of the window in kilobases.
    - output_file (str): Optional path to save results to a CSV.
    """
    # 1. Calculate interval boundaries in base pairs
    size_bp = int(size_kb * 1000)
    half_window = size_bp // 2

    interval_start = max(0, position - half_window)
    interval_end = position + half_window

    print(f"Searching interval: {chromosome}:{interval_start}-{interval_end}")

    # 2. Load BED file
    # BED format standard: chrom, chromStart, chromEnd, name, ...
    # We only read the first 4 columns to save memory
    try:
        bed_df = pd.read_csv(
            bed_file_path,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "gene_name"],
            usecols=[0, 1, 2, 3],
            dtype={"chrom": str, "start": int, "end": int, "gene_name": str},
        )
    except Exception as e:
        print(f"Error reading BED file: {e}")
        return []

    # 3. Filter by chromosome first (massive reduction in data size)
    chrom_df = bed_df[bed_df["chrom"] == str(chromosome)]

    if chrom_df.empty:
        print(f"No entries found for chromosome {chromosome}.")
        return []

    # 4. Find overlaps
    # A gene overlaps if its start is before our end AND its end is after our start
    overlap_condition = (chrom_df["start"] < interval_end) & (
        chrom_df["end"] > interval_start
    )
    matching_genes_df = chrom_df[overlap_condition]

    # 5. Extract unique gene names
    # Drop duplicates in case a gene has multiple transcript isoforms listed
    unique_genes = matching_genes_df["gene_name"].dropna().unique().tolist()

    print(f"Found {len(unique_genes)} gene(s).")

    # Optional: Save to file
    if output_file and unique_genes:
        matching_genes_df.to_csv(output_file, index=False, sep="\t")
        print(f"Detailed overlap data saved to {output_file}")

    return unique_genes

def match_vcf_intervals(vcf_file_path, bed_file_path):
  """Returns a pd.Dataframe contianing the vcf position and 

    Parameters:
    - bed_file_path (str): Path to the BED file (tab-separated, no header).
    - chromosome (str): e.g., 'chr1' or '1' (must match BED file format).
    - position (int): The center or start position in base pairs.
    - size_kb (float): The size of the window in kilobases.
    - output_file (str): Optional path to save results to a CSV.
    """

# --- Example Usage ---
if __name__ == "__main__":
    # Example setup (Replace with your actual file path)
    bed_path = "hg38_genes.bed"

    # Define your region of interest
    chrom = "chr17"
    pos = 43044295  # Example: Near the BRCA1 region
    window_kb = 500  # 500 KB window around the position

    # Run the function
    genes = find_genes_in_interval(
        bed_file_path=bed_path,
        chromosome=chrom,
        position=pos,
        size_kb=window_kb,
    )

    print("\nGenes found in interval:")
    print(genes)