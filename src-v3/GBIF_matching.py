from pygbif import species
import pandas as pd
import os

def parse_taxonomy_string(tax_string):
    """Parse a taxonomy string with rank prefixes into a dict"""
    ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']  # domain to species
    tax_dict = {}
    
    parts = tax_string.split(';')
    for part in parts:
        part = part.strip()
        if not part:  # Handle empty fields
            continue
            
        # Extract rank and name
        if '__' in part:
            rank, name = part.split('__')
            if rank in ranks:
                tax_dict[rank] = name

    return tax_dict

def clean_species_name(name):
    """Remove common prefixes/suffixes that won't match in GBIF"""
    removals = [
        'uncultured_',
        'uncultured ',
        '_bacterium',
        'metagenome',
        'environmental'
    ]
    name = name.strip()
    for r in removals:
        name = name.replace(r, '')
    return name

def match_to_gbif(taxonomy_dict):
    """Match a taxonomy dict against GBIF"""
    # Start with most specific rank available
    for rank in ['s', 'g', 'f', 'o', 'c', 'p', 'd']:
        if rank in taxonomy_dict:
            name = taxonomy_dict[rank]
            name = clean_species_name(name)  # Removing species name- I've read GBIF species may be more useful than WoRMS?
            
            if not name:  # Skip if name is empty after cleaning
                continue
                
            try:
                # Try exact match first
                results = species.name_backbone(
                    name=name,
                    rank=rank,
                    strict=True  # Only match at this rank
                )
                
                # If no exact match, try fuzzy match
                if not results or results.get('matchType') == 'NONE':
                    results = species.name_backbone(
                        name=name,
                        rank=rank,
                        strict=False  # Allow matching to higher ranks
                    )
                
                # If still no match, try name_suggest as last resort
                if not results or results.get('matchType') == 'NONE':
                    suggestions = species.name_suggest(q=name)
                    if suggestions:
                        # Use name_backbone to get full classification
                        results = species.name_backbone(
                            name=suggestions[0]['scientificName']
                        )
                
                if results and results.get('matchType') != 'NONE':
                    classification = {
                        'kingdom': results.get('kingdom'),
                        'phylum': results.get('phylum'),
                        'class': results.get('class'),
                        'order': results.get('order'),
                        'family': results.get('family'),
                        'genus': results.get('genus')
                    }
                    
                    return {
                        'old_taxonRank': rank,
                        'old name': name,
                        'scientificName': results['scientificName'],
                        'scientificNameID': str(results['usageKey']),
                        'taxonRank': results['rank'],
                        'matchType': results['matchType'],  # Added to track match quality
                        **classification
                    }
            except Exception as e:
                print(f"Error matching {name}: {e}")
                
    return None

# Read test data
df = pd.read_csv('edna2obis/raw/test_data-taxa_sample_table_l7.tsv', 
                 sep='\t', 
                 skiprows=[1],  # Skip the "Constructed from biom file" line
                 comment='#',   # This will skip the "#OTU ID" line
                 nrows=25,   
                 index_col=0)   # First column is index

# Process first 5 rows
results = []
for tax_string in df.index:  # Iterate over index which contains our taxonomy strings
    tax_dict = parse_taxonomy_string(tax_string)
    match = match_to_gbif(tax_dict)
    if match:
        results.append({
            'full_tax': tax_string,
            'verbatimIdentification': tax_string,
            **match  # Unpack the match results
        })
    else:
        # Handle no match case
        results.append({
            'full_tax': tax_string,
            'verbatimIdentification': tax_string,
            'old_taxonRank': None,
            'old name': None,
            'scientificName': None,
            'scientificNameID': None,
            'kingdom': None,
            'phylum': None,
            'class': None,
            'order': None,
            'family': None,
            'genus': None,
            'taxonRank': None
        })

# Create DataFrame with columns in same order as WoRMS output
output_df = pd.DataFrame(results)[['full_tax', 'verbatimIdentification', 'old_taxonRank', 
                                 'old name', 'scientificName', 'scientificNameID',
                                 'kingdom', 'phylum', 'class', 'order', 'family', 
                                 'genus', 'taxonRank']]

# Add this before saving to check what's happening
print(f"Saving {len(output_df)} matched records to edna2obis/processed/gbif_test_matches_full.tsv")
print(f"First few results: {output_df.head(2)}")

# Create directory if it doesn't exist
os.makedirs('edna2obis/processed/', exist_ok=True)

output_df.to_csv('edna2obis/processed/gbif_test_matches_full.tsv', sep='\t', index=False)
print("File saved successfully!")