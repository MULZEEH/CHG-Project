# Python scripts file with queries definitions for SnpSift filtering

# first query
def get_somatic_query(wildcards):
    # You can even make this dynamic based on the sample group!
    quality = 30
    impacts = "('HIGH', 'MODERATE')"
    
    query = (
        f"(QUAL > {quality}) & "
        f"(ANN[*].IMPACT has {impacts}) & "
        "(isHomRef(gen[0])) & (isHet(gen[1]))"
    )
    return query
