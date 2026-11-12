# Taxonomic mapping

**FAQ: Why “mapping” and not “assignment”?**  
Taxonomic assignment (matching ASVs/OTUs to taxa) happens before edna2obis. Taxonomic mapping to WoRMS (or GBIF) takes those names and maps them to backbone records so the data can be submitted to OBIS/GBIF.

## WoRMS

- Clean `verbatimIdentification`: split by `;`, replace `_ - /` with spaces, trim, remove trailing “sp.”/“spp.” and numbers, drop empty/“unassigned”.
- If it exists as a species by AphiaID in the local reference database file provided (optionally), look up AphiaID directly and fill fields.
- Otherwise, use WoRMS “Match Names”:
  - Start from the most specific level (ex: species) and move broader only if needed.
  - Keep only records with status “accepted” and add every accepted alternative to the review list.
  - Use the most specific level that returns an accepted record for the mapping (ex: favors a genus match over a class match).
- If an assay is set to skip species, drop the last level before searching.
- If nothing accepted, set `scientificName=incertae sedis` and clear ranks.
- Fill: `scientificName`, `scientificNameID` (LSID), `taxonRank`, ranks (`kingdom … species`), `nameAccordingTo=WoRMS`, plus `cleanedTaxonomy` and `match_type_debug`.
- Endpoints: AphiaID record; Match Names.

## GBIF

- Clean `verbatimIdentification` as above and remove rank prefixes like `d__`, `p__`, etc.
- Try the most specific level first; move broader only if needed.
- For each candidate from the initial lookup, confirm against the GBIF backbone to get a confidence score.
- Choose the highest-confidence accepted match by default; optionally choose by lineage `assignment_score` first when `gbif_use_assignment_score_for_selection` is true. Include accepted candidates and unaccepted/synonym review rows in the review list. If none, set `scientificName=incertae sedis`.
- If an assay is set to skip species, drop the last level before searching.
- Cache results to speed up future runs.
- Fill: `scientificName`, `taxonID` (e.g., `gbif:<usageKey>`), `taxonRank`, ranks (`kingdom … genus`), `confidence`, optional `higherClassification`, `environment` when available, `nameAccordingTo=GBIF`, plus `cleanedTaxonomy` and `match_type_debug`.
- Endpoints: name_lookup; name_backbone; name_usage parents.

## Settings that affect matching

- `taxonomic_api_source`: `"WoRMS"` or `"GBIF"`
- `assays_to_skip_species_match`: assays where species level is ignored
- Per-assay maximum depth: inferred from your inputs; used for species skipping
- GBIF candidate limit and caching
- `gbif_use_assignment_score_for_selection`: GBIF only. If true, lineage `assignment_score` is used before GBIF confidence when choosing `selected_match`
- `gbif_return_higher_classification`: GBIF only. If true, calls the parents endpoint and includes `higherClassification`
- Local reference database (WoRMS only) for fast AphiaID lookups

Also writes `taxa_mapping_INFO_<API>.xlsx` (sheet `taxa_mapping_INFO`). One row can appear per candidate match, not only the match used in the occurrence file.

## `taxa_mapping_INFO` columns

**Important column!** `selected_match` is True on the row the pipeline treats as the official mapping for each `verbatimIdentification` (the one reflected in the occurrence file). For `taxamapper.py`-only runs, it is True on the row that *would* be chosen in the full pipeline. When you summarize results, **filter on** `selected_match` **= TRUE**; other mapping rows exist for transparency and QA.

Columns are ordered with `selected_match` immediately after `cleanedTaxonomy` in the `.xlsx` file.

| Column                 | Meaning                                                                                                                                                                                                                               |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| verbatimIdentification | Original taxonomic string from your data before matching.                                                                                                                                                                             |
| cleanedTaxonomy        | Normalized / cleaned version of verbatimIdentification (what is actually used for API lookup).                                                                                                                                         |
| selected_match         | **Indicates the taxonomic mapping!** True only on the chosen mapping row for that `verbatimIdentification` (occurrence file in full pipeline; same rule for taxamapper-only output). False on alternate candidate rows kept for review. |
| ambiguous              | WoRMS and GBIF. True when more than one candidate mapping was found for this verbatim string.                                                                                                                                         |
| replaced_unaccepted    | WoRMS and GBIF. True when an unaccepted/synonym name was resolved to an accepted valid name.                                                                                                                                          |
| unaccepted_match       | WoRMS and GBIF. True when the mapping in that row is an unaccepted/synonym name shown for comparison.                                                                                                                                 |
| ranks_matched          | WoRMS and GBIF. Number of taxonomic ranks given in your verbatimIdentification that matched the candidate classification.                                                                                                             |
| ranks_provided         | WoRMS and GBIF. Number of taxonomic ranks given in your verbatimIdentification for that taxonomy.                                                                                                                                     |
| assignment_score       | WoRMS and GBIF. Ratio of ranks_matched / ranks_provided.                                                                                                                                                                              |
| environment            | WoRMS and GBIF habitat labels when available. GBIF habitat values come from GBIF usage data and may be coarser than WoRMS.                                                                                                            |
| scientificName         | Scientific name for this candidate row from the backbone.                                                                                                                                                                             |
| confidence             | GBIF only. Assignment confidence score from 0 to 100.                                                                                                                                                                                 |
| taxonRank              | Rank of scientificName assigned.                                                                                                                                                                                                      |
| scientificNameID       | WoRMS only. LSID for the taxon.                                                                                                                                                                                                       |
| taxonID                | GBIF only. Backbone taxonomic identifier.                                                                                                                                                                                             |
| higherClassification   | WoRMS and optionally GBIF. Pipe-separated lineage. GBIF uses the name_usage parents endpoint when `gbif_return_higher_classification` is true.                                                                                        |
| kingdom                | Darwin Core kingdom.                                                                                                                                                                                                                  |
| phylum                 | Darwin Core phylum.                                                                                                                                                                                                                   |
| class                  | Darwin Core class.                                                                                                                                                                                                                    |
| order                  | Darwin Core order.                                                                                                                                                                                                                    |
| family                 | Darwin Core family.                                                                                                                                                                                                                   |
| genus                  | Darwin Core genus.                                                                                                                                                                                                                    |
| match_type_debug       | Short code for how this row was matched via the API.                                                                                                                                                                                  |
| nameAccordingTo        | Which API you used to assign taxonomy.                                                                                                                                                                                                |
