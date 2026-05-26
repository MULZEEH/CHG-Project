# find_gene.py
# this script holds some utility function to find the gene name 
# given the chromosome and the position of the variant. It uses the GTF file to find the gene name.

# the function
import pandas as pd

INTORNO_SIZE = 1,000,000

def get_gene_name(chromosome: str, position: int, gtf_file: str) -> str:
    """Returns the gene name given the chromosome and the position of the variant."""
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    # need to check the parsing of the right gtf information in the miao miao
    gtf.columns = ['CHROM', 'SOURCE', 'FEATURE', 'START', 'END', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTES']
    gene_name = gtf[(gtf['CHROM'] == chromosome) & (gtf['START'] <= position) & (gtf['END'] >= position) & (gtf['FEATURE'] == 'gene')]['ATTRIBUTES'].str.extract('gene_name "([^"]+)"')[0]
    if gene_name.empty:
        return None
    return gene_name.iloc[0]

# function that wraps the get_gene_name to be applied to a vcf file
def get_gene_name_vcf(vcf: pd.DataFrame, gtf_file: str) -> pd.DataFrame:
    """Returns a dataframe with the gene names for each variant in the VCF file."""
    vcf['GENE_NAME'] = vcf.apply(lambda row: get_gene_name(row['CHROM'], row['POS'], gtf_file), axis=1)
    return vcf

# review this ufnction: IT WAS 100% VIBECODED
def peak_finding(chromosome: str, position: int, gtf_file: str) -> pd.DataFrame:
    """Returns a dataframe with the gene names for each variant in the VCF file."""
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    gtf.columns = ['CHROM', 'SOURCE', 'FEATURE', 'START', 'END', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTES']
    peaks = gtf[(gtf['CHROM'] == chromosome) & (gtf['START'] <= position) & (gtf['END'] >= position) & (gtf['FEATURE'] == 'peak')]
    peaks['GENE_NAME'] = peaks['ATTRIBUTES'].str.extract('gene_name "([^"]+)"')[0]
    return peaks[['CHROM', 'POS', 'GENE_NAME']]

# need to handle the fact that most of the peaks are not in coding regions but still affects the transcription rate (that is measureble using alpha genome)
# idea is to use a backup function inc ase a vcf record is not inserted inside a coding region and define peaks 

# this function takes in input a chromosomic postion and a gtf file:
# the function returns a tuple of up to 2 gene name, with these codntion:
# are the 2 gene closest to the 
def close_gene_name(chromosome: str, position: int, gtf_file: str) -> tuple[str, str]:

    gtf[DISTANCE] = min(gtf['END']+position, gtf['START']-position)

    gene_name_upstream = gtf[(gtf['CHROM'] == chromosome) & (gtf['END']+INTORNO_SIZE >= position) & (gtf['FEATURE'] == 'gene')]['ATTRIBUTES'].str.extract('gene_name "([^"]+)"')[0]
    gene_name_downstream = gtf[(gtf['CHROM'] == chromosome) & (gtf['START'] > INTORNO_SIZE) (gtf['START']-INTORNO_SIZE <= position) & (gtf['FEATURE'] == 'gene')]['ATTRIBUTES'].str.extract('gene_name "([^"]+)"')[0]
    

def test():
    