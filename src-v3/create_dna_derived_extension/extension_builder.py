"""
DNA Derived Extension Builder for edna2obis
Contains the complex logic for building DNA derived extensions from ASV data
Create DNA derived extension file
"""

import pandas as pd
import numpy as np
import os


def create_dna_derived_extension(params, data, raw_data_tables, dwc_data, occurrence_core, all_processed_occurrence_dfs, checklist_df, reporter):
    """Create DNA derived extension file"""
    try:
        reporter.add_section("Creating DNA Derived Extension")

        # --- Load Mapper Config (MOVED TO TOP) ---
        # This must be loaded first to drive all subsequent logic
        import yaml
        with open('data_mapper.yaml', 'r', encoding='utf-8') as f:
            full_mapper = yaml.safe_load(f)
        
        format_prefix = "generic_" if params.get('metadata_format') == 'GENERIC' else ""
        dna_derived_key = f"{format_prefix}dna_derived_extension"
        dna_derived_map_config = full_mapper.get(dna_derived_key, {})
        DESIRED_DNA_DERIVED_COLUMNS = list(dna_derived_map_config.keys())

        # --- START FROM OCCURRENCE CORE ---
        # DNA derived extension should have ONE ROW PER OCCURRENCE (one-to-one with occurrence core)
        # NOT one row per unique ASV!
        reporter.add_text("Building DNA derived extension from occurrence core (one row per occurrence)...")
        
        # Determine which ID field to use (GBIF uses taxonID, WoRMS uses scientificNameID)
        id_field = None
        if 'taxonID' in occurrence_core.columns and not occurrence_core['taxonID'].isna().all():
            id_field = 'taxonID'
            reporter.add_text("Using 'taxonID' field (GBIF workflow)")
        elif 'scientificNameID' in occurrence_core.columns and not occurrence_core['scientificNameID'].isna().all():
            id_field = 'scientificNameID'
            reporter.add_text("Using 'scientificNameID' field (WoRMS workflow)")
        else:
            reporter.add_error("Could not find 'taxonID' or 'scientificNameID' in occurrence core. Cannot build DNA derived extension.")
            return
        
        # Start with occurrence core and extract the featureID from the ID field
        dna_derived_base = occurrence_core[['occurrenceID', id_field, 'eventID', 'parentEventID', 'materialSampleID', 'assay_name']].copy()
        dna_derived_base['featureID'] = dna_derived_base[id_field].str.replace('ASV:', '')
        
        reporter.add_text(f"Starting with {len(dna_derived_base)} occurrence records.")
        
        # The dataframe is now passed in with all necessary columns from the occurrence_builder's
        # pass-through logic. We no longer need to perform complex merges here.
        dna_derived_df_final = dna_derived_base.copy()

        # Bring through pass-through fields preserved in the occurrence core
        try:
            passthrough_terms = []
            for dwc_term, mapping_info in dna_derived_map_config.items():
                if isinstance(mapping_info, dict) and mapping_info.get('source') in ['sampleMetadata', 'experimentRunMetadata']:
                    faire_term = mapping_info.get('faire_term')
                    if faire_term and isinstance(faire_term, str):
                        passthrough_terms.append(faire_term)
                        # include unit columns if present
                        passthrough_terms.append(f"{faire_term}_unit")

            # Get all available columns from occurrence core that match our passthrough terms
            available_passthrough = [c for c in set(passthrough_terms) if c in occurrence_core.columns]
            if available_passthrough:
                # Merge only the columns we need, avoiding duplicates
                cols_to_merge = ['occurrenceID'] + [c for c in available_passthrough if c not in dna_derived_df_final.columns]
                if len(cols_to_merge) > 1:  # More than just occurrenceID
                    occ_subset = occurrence_core[cols_to_merge].copy()
                    dna_derived_df_final = dna_derived_df_final.merge(occ_subset, on='occurrenceID', how='left')
                    reporter.add_text(f"Brought forward {len(available_passthrough)} pass-through field(s) from occurrence core.")
                else:
                    reporter.add_text("No additional pass-through fields needed from occurrence core.")
            else:
                reporter.add_text("No pass-through fields found in occurrence core.")
        except Exception as _e:
            reporter.add_text(f"Could not merge pass-through fields from occurrence core: {_e}")

        # Now merge DNA sequences from taxonomy files
        # Build a lookup table of (featureID, assay_name) -> DNA_sequence
        all_tax_data = []
        for run_name, tables in raw_data_tables.items():
            if 'taxonomy' in tables and not tables['taxonomy'].empty:
                tax_copy = tables['taxonomy'].copy()
                tax_copy['_assay_name'] = params['datafiles'][run_name].get('assay_name', 'unknown')
                all_tax_data.append(tax_copy)
        
        if not all_tax_data:
            reporter.add_error("No taxonomy data found to build DNA derived extension.")
            return
        
        combined_tax_data = pd.concat(all_tax_data, ignore_index=True)
        
        # Find DNA sequence column
        dna_col = None
        for col in combined_tax_data.columns:
            if str(col).lower().strip() == 'dna_sequence':
                dna_col = col
                break
        
        if not dna_col:
            reporter.add_error(f"Could not find dna_sequence column in taxonomy files. Found columns: {list(combined_tax_data.columns)}")
            return
        
        # Create lookup: keep featureid, dna_sequence, and assay
        feature_col = None
        for col in combined_tax_data.columns:
            if str(col).lower().strip() == 'featureid':
                feature_col = col
                break
        
        if not feature_col:
            reporter.add_error(f"Could not find featureid column in taxonomy files.")
            return
        
        dna_lookup = combined_tax_data[[feature_col, dna_col, '_assay_name']].copy()
        dna_lookup.rename(columns={feature_col: 'featureID', dna_col: 'DNA_sequence'}, inplace=True)
        dna_lookup.drop_duplicates(subset=['featureID', '_assay_name'], keep='first', inplace=True)
        
        # Merge DNA sequences into the base (one DNA sequence per occurrence)
        dna_derived_with_occurrence = dna_derived_df_final.merge(
            dna_lookup,
            left_on=['featureID', 'assay_name'],
            right_on=['featureID', '_assay_name'],
            how='left'
        )
        
        # Clean up internal tracking columns
        internal_cleanup_cols = ['_assay_name', 'taxonID', 'scientificNameID', 'featureID']
        for col in internal_cleanup_cols:
            if col in dna_derived_with_occurrence.columns:
                dna_derived_with_occurrence.drop(columns=[col], inplace=True)

        if 'featureID' in dna_derived_with_occurrence.columns:
            dna_derived_with_occurrence.drop(columns=['featureID'], inplace=True)

        dna_derived_df_final = dna_derived_with_occurrence.copy()

        # Rename materialSampleID to source_mat_id for final output.
        if 'materialSampleID' in dna_derived_df_final.columns:
            dna_derived_df_final.rename(columns={'materialSampleID': 'source_mat_id'}, inplace=True)
        
        # If source_mat_id is missing or empty, fill it from sampleMetadata via parentEventID -> samp_name
        if ('source_mat_id' not in dna_derived_df_final.columns) or (dna_derived_df_final['source_mat_id'].isna().all()):
            try:
                if 'sampleMetadata' in data and not data['sampleMetadata'].empty:
                    sm_ms = data['sampleMetadata'][['samp_name','materialSampleID']].copy()
                    sm_ms.rename(columns={'materialSampleID':'source_mat_id'}, inplace=True)
                    dna_derived_df_final = dna_derived_df_final.merge(sm_ms, left_on='parentEventID', right_on='samp_name', how='left')
                    # Prefer existing non-null values
                    if 'source_mat_id_x' in dna_derived_df_final.columns and 'source_mat_id_y' in dna_derived_df_final.columns:
                        dna_derived_df_final['source_mat_id'] = dna_derived_df_final['source_mat_id_x'].combine_first(dna_derived_df_final['source_mat_id_y'])
                        dna_derived_df_final.drop(columns=['source_mat_id_x','source_mat_id_y'], inplace=True)
                    # Drop helper join column
                    if 'samp_name' in dna_derived_df_final.columns:
                        # keep original parentEventID; samp_name was only for join
                        pass
            except Exception:
                pass
        
        # Ensure no duplicate column labels (keep first occurrence)
        dna_derived_df_final = dna_derived_df_final.loc[:, ~dna_derived_df_final.columns.duplicated(keep='first')]
        
        reporter.add_text(f"Shape after all merges: {dna_derived_df_final.shape}")
        
        # --- STEP W: Merge Metadata using eventID (REVISED LOGIC) ---
        reporter.add_section("Merging Sample and Experiment Metadata", level=3)
        try:
            # 1. Check for necessary dataframes
            if 'experimentRunMetadata' not in data or data['experimentRunMetadata'].empty:
                raise ValueError("experimentRunMetadata sheet is missing or empty.")
            if 'sampleMetadata' not in data or data['sampleMetadata'].empty:
                raise ValueError("sampleMetadata sheet is missing or empty.")
            
            erm_df = data['experimentRunMetadata'].copy()
            sm_df = data['sampleMetadata'].copy()

            # 2. Prepare for robust joins by creating clean, standardized join keys
            dna_derived_df_final['eventID_join_key'] = dna_derived_df_final['eventID'].astype(str).str.strip()
            erm_df['lib_id_join_key'] = erm_df['lib_id'].astype(str).str.strip()
            sm_df['samp_name_join_key'] = sm_df['samp_name'].astype(str).str.strip()
            
            # Diagnostics: eventID vs lib_id
            try:
                event_keys = set(dna_derived_df_final['eventID_join_key'].dropna().astype(str).unique())
                lib_keys = set(erm_df['lib_id_join_key'].dropna().astype(str).unique())
                overlap_evt_lib = event_keys & lib_keys
                reporter.add_text(f"eventID/lib_id diagnostics: DNA rows={len(dna_derived_df_final)}, unique eventID={len(event_keys)}, unique lib_id={len(lib_keys)}, overlaps={len(overlap_evt_lib)}")
                if len(overlap_evt_lib) == 0:
                    sample_evt = list(sorted(event_keys))[:10]
                    sample_lib = list(sorted(lib_keys))[:10]
                    reporter.add_warning("No overlap between eventID and lib_id values. Merge will produce nulls.")
                    reporter.add_list(sample_evt, title="Sample eventID values (first 10):")
                    reporter.add_list(sample_lib, title="Sample lib_id values (first 10):")
            except Exception as diag_e:
                reporter.add_warning(f"Could not compute eventID/lib_id diagnostics: {diag_e}")

            # 3. Add all fields from experimentRunMetadata via eventID -> lib_id
            merged_df = dna_derived_df_final.merge(
                erm_df,
                left_on='eventID_join_key',
                right_on='lib_id_join_key',
                how='left',
                suffixes=('', '_from_erm')
            )
            # Post-merge diagnostics
            try:
                matched_erm = merged_df['lib_id_join_key'].notna().sum()
                reporter.add_text(f"ERM merge: matched rows={matched_erm} / {len(merged_df)}")
            except Exception:
                pass

            # Determine which samp_name column to use for the next merge
            samp_col_to_use = 'samp_name_from_erm' if 'samp_name_from_erm' in merged_df.columns else ('samp_name' if 'samp_name' in merged_df.columns else None)
            if samp_col_to_use is None:
                reporter.add_warning("No 'samp_name' column found after ERM merge; cannot merge sampleMetadata.")
                samp_col_to_use = 'samp_name'
                merged_df[samp_col_to_use] = pd.NA

            # 4. Use samp_name from the first merge to bring in all of sampleMetadata
            merged_df['samp_name_join_key'] = merged_df[samp_col_to_use].astype(str).str.strip()

            # Diagnostics: samp_name vs sampleMetadata
            try:
                samp_keys = set(merged_df['samp_name_join_key'].dropna().astype(str).unique())
                sm_keys = set(sm_df['samp_name_join_key'].dropna().astype(str).unique())
                overlap_samp = samp_keys & sm_keys
                reporter.add_text(f"samp_name diagnostics: unique samp_name(after ERM)={len(samp_keys)}, unique samp_name(sampleMetadata)={len(sm_keys)}, overlaps={len(overlap_samp)}")
                if len(overlap_samp) == 0:
                    sample_samp = list(sorted(samp_keys))[:10]
                    sample_sm = list(sorted(sm_keys))[:10]
                    reporter.add_warning("No overlap between samp_name from ERM and sampleMetadata.samp_name. Merge will produce nulls.")
                    reporter.add_list(sample_samp, title="Sample samp_name from ERM (first 10):")
                    reporter.add_list(sample_sm, title="Sample samp_name from sampleMetadata (first 10):")
            except Exception as diag2_e:
                reporter.add_warning(f"Could not compute samp_name diagnostics: {diag2_e}")

            final_merged_df = merged_df.merge(
                sm_df,
                on='samp_name_join_key',
                how='left',
                suffixes=('', '_from_sm')
            )
            # Post-merge diagnostics for sampleMetadata
            try:
                matched_sm = final_merged_df['samp_name_join_key'].notna().sum()
                crit_cols = [
                    'materialSampleID', 'depth_category', 'decimalLatitude', 'decimalLongitude',
                    'temp', 'salinity', 'samp_vol_we_dna_ext', 'concentration', 'concentration_unit'
                ]
                coverage_msgs = []
                for cc in crit_cols:
                    if cc in final_merged_df.columns:
                        na_ratio = float(final_merged_df[cc].isna().mean())
                        coverage_msgs.append(f"{cc}: {(1.0 - na_ratio)*100:.1f}% filled")
                    else:
                        coverage_msgs.append(f"{cc}: MISSING COLUMN")
                reporter.add_text(f"SampleMetadata merge: used column '{samp_col_to_use}', matched rows (by key)={matched_sm} / {len(final_merged_df)}")
                reporter.add_list(coverage_msgs, title="Post-merge coverage for key sampleMetadata fields:")
            except Exception:
                pass

            reporter.add_text("Merged with experimentRunMetadata and sampleMetadata.")

            # 5. Assign the fully merged dataframe back
            dna_derived_df_final = final_merged_df.copy()

            # 5b. Alias suffixed columns back to their base faire_term names for mapper compatibility
            try:
                aliased_count_sm = 0
                aliased_count_erm = 0
                # SampleMetadata
                for dwc_term, mapping_info in dna_derived_map_config.items():
                    if not isinstance(mapping_info, dict): continue
                    if mapping_info.get('source') == 'sampleMetadata':
                        ft = mapping_info.get('faire_term')
                        if not isinstance(ft, str) or not ft.strip(): continue
                        src_col = f"{ft}_from_sm"
                        if src_col in dna_derived_df_final.columns:
                            if ft not in dna_derived_df_final.columns:
                                dna_derived_df_final[ft] = dna_derived_df_final[src_col]
                            else:
                                dna_derived_df_final[ft] = dna_derived_df_final[ft].combine_first(dna_derived_df_final[src_col])
                            aliased_count_sm += 1
                # ERM
                for dwc_term, mapping_info in dna_derived_map_config.items():
                    if not isinstance(mapping_info, dict): continue
                    if mapping_info.get('source') == 'experimentRunMetadata':
                        ft = mapping_info.get('faire_term')
                        if not isinstance(ft, str) or not ft.strip(): continue
                        src_col = f"{ft}_from_erm"
                        if src_col in dna_derived_df_final.columns:
                            if ft not in dna_derived_df_final.columns:
                                dna_derived_df_final[ft] = dna_derived_df_final[src_col]
                            else:
                                dna_derived_df_final[ft] = dna_derived_df_final[ft].combine_first(dna_derived_df_final[src_col])
                            aliased_count_erm += 1
                # Before dropping suffixed columns, alias any unit columns to their base names
                # Example: 'samp_vol_we_dna_ext_unit_from_sm' -> 'samp_vol_we_dna_ext_unit'
                try:
                    # Prefer values from sampleMetadata over ERM when both exist
                    unit_from_sm = [c for c in dna_derived_df_final.columns if c.endswith('_unit_from_sm')]
                    unit_from_erm = [c for c in dna_derived_df_final.columns if c.endswith('_unit_from_erm')]
                    for c in unit_from_sm:
                        base = c.replace('_from_sm', '')
                        if base not in dna_derived_df_final.columns:
                            dna_derived_df_final[base] = dna_derived_df_final[c]
                        else:
                            dna_derived_df_final[base] = dna_derived_df_final[base].combine_first(dna_derived_df_final[c])
                    for c in unit_from_erm:
                        base = c.replace('_from_erm', '')
                        if base not in dna_derived_df_final.columns:
                            dna_derived_df_final[base] = dna_derived_df_final[c]
                        else:
                            dna_derived_df_final[base] = dna_derived_df_final[base].combine_first(dna_derived_df_final[c])
                except Exception as unit_alias_e:
                    reporter.add_warning(f"Unit column aliasing failed: {unit_alias_e}")

                # Drop leftover suffixed columns to avoid confusion
                drop_cols = [c for c in dna_derived_df_final.columns if c.endswith('_from_sm') or c.endswith('_from_erm')]
                if drop_cols:
                    dna_derived_df_final.drop(columns=drop_cols, inplace=True, errors='ignore')
                reporter.add_text(f"Aliased fields: sampleMetadata={aliased_count_sm}, experimentRunMetadata={aliased_count_erm}")
            except Exception as alias_e:
                reporter.add_warning(f"Aliasing step failed: {alias_e}")
            
            # 6. Cleanup join keys
            dna_derived_df_final.drop(columns=['eventID_join_key', 'lib_id_join_key', 'samp_name_join_key'], inplace=True, errors='ignore')

        except Exception as e:
            reporter.add_warning(f"Metadata merge failed with the new logic: {e}")
            import traceback
            reporter.add_warning(f"Traceback: {traceback.format_exc()}")

        # --- STEP X: Handle Units (Safe and non-destructive) ---
        reporter.add_text("Applying specific unit logic...")
        
        # This logic runs after the main data merge. It proactively looks for original
        # unit columns (e.g., 'samp_size_unit') and prepares them for combination.
        # It does NOT rely on the data_mapper.yaml to load these, making it robust.
        
        # 1) concentrationUnit: Stays as a separate column.
        if 'concentration_unit' in dna_derived_df_final.columns:
            dna_derived_df_final['concentrationUnit'] = dna_derived_df_final['concentration_unit'].fillna('ng/µl')
            reporter.add_text("Populated 'concentrationUnit' from source 'concentration_unit' column.")
        else:
            dna_derived_df_final['concentrationUnit'] = 'ng/µl'
            reporter.add_text("Set 'concentrationUnit' to default 'ng/µl'.")
        
        # 2) size_fracUnit: To be combined. Check for 'size_frac_unit' first.
        if 'size_frac_unit' in dna_derived_df_final.columns:
            dna_derived_df_final['size_fracUnit'] = dna_derived_df_final['size_frac_unit'].fillna('micrometer')
            reporter.add_text("Populated 'size_fracUnit' from source 'size_frac_unit' column.")
        else:
            dna_derived_df_final['size_fracUnit'] = 'micrometer'
            reporter.add_text("Set 'size_fracUnit' to default 'micrometer' (source column not found).")
        
        # 3) samp_vol_we_dna_extUnit: To be combined. Check for 'samp_vol_we_dna_ext_unit' first.
        if 'samp_vol_we_dna_ext_unit' in dna_derived_df_final.columns:
            # Use source data, filling any blanks with the default 'mL'.
            dna_derived_df_final['samp_vol_we_dna_extUnit'] = dna_derived_df_final['samp_vol_we_dna_ext_unit'].fillna('mL').astype(str).str.replace('ml', 'mL', case=False)
            reporter.add_text("Populated 'samp_vol_we_dna_extUnit' from source 'samp_vol_we_dna_ext_unit' column.")
        else:
            dna_derived_df_final['samp_vol_we_dna_extUnit'] = 'mL'
            reporter.add_text("Set 'samp_vol_we_dna_extUnit' to default 'mL' (source column not found).")
        
        # 4) samp_size_unit: To be combined. Check for 'samp_size_unit' first. No default.
        if 'samp_size_unit' in dna_derived_df_final.columns:
            dna_derived_df_final['samp_size_unit'] = dna_derived_df_final['samp_size_unit'].fillna('').astype(str).replace({'nan': '', 'None': ''})
            reporter.add_text("Populated 'samp_size_unit' from source 'samp_size_unit' column.")
        else:
            dna_derived_df_final['samp_size_unit'] = ''
            reporter.add_warning("'samp_size_unit' column not found in metadata. Units for 'samp_size' will be empty.")
            
        # --- Combine value and unit columns ---
        reporter.add_text("Combining specified value and unit fields...")
        
        # Define which value/unit pairs to combine
        pairs_to_combine = {
            'size_frac': 'size_fracUnit',
            'samp_vol_we_dna_ext': 'samp_vol_we_dna_extUnit',
            'samp_size': 'samp_size_unit'
        }
        
        cols_to_drop_after_combine = []
        
        for value_col, unit_col in pairs_to_combine.items():
            if value_col in dna_derived_df_final.columns and unit_col in dna_derived_df_final.columns:
                # Convert both to string and handle missing values gracefully
                value_series = dna_derived_df_final[value_col].astype(str).replace({'nan': '', 'None': ''})
                unit_series = dna_derived_df_final[unit_col].astype(str).replace({'nan': '', 'None': ''})
                
                # Combine them, adding a space only if both exist
                combined_series = value_series.str.strip() + ' ' + unit_series.str.strip()
                dna_derived_df_final[value_col] = combined_series.str.strip()
                
                cols_to_drop_after_combine.append(unit_col)
                reporter.add_success(f"Combined '{value_col}' and '{unit_col}'.")
        
        # Drop the unit columns that have been combined
        if cols_to_drop_after_combine:
            dna_derived_df_final.drop(columns=cols_to_drop_after_combine, inplace=True)
            reporter.add_text(f"Dropped combined unit columns: {cols_to_drop_after_combine}")
    
        # --- STEP Y: Populate fields from projectMetadata & analysisMetadata (REWRITTEN) ---
        reporter.add_text("Processing project and analysis metadata fields...")
        
        # Helper: resolve actual column names in a source DataFrame for given faire_terms
        def _resolve_source_columns(source_df: pd.DataFrame, requested_terms: list[str]) -> dict:
            resolved: dict = {}
            if not isinstance(source_df, pd.DataFrame) or source_df.empty or not requested_terms:
                return resolved
            # Build normalization map of columns
            def _norm(s: str) -> str:
                return ''.join(ch for ch in str(s).lower().strip().replace(' ', '').replace('\t','').replace('\n','') if ch.isalnum() or ch == '_')
            col_norm_map = {_norm(c): c for c in source_df.columns}
            for term in requested_terms:
                if not isinstance(term, str):
                    continue
                # 1) direct exact
                if term in source_df.columns:
                    resolved[term] = term
                    continue
                # 2) case-insensitive exact
                lower_map = {str(c).lower(): c for c in source_df.columns}
                if term.lower() in lower_map:
                    resolved[term] = lower_map[term.lower()]
                    continue
                # 3) normalized match
                tnorm = _norm(term)
                if tnorm in col_norm_map:
                    resolved[term] = col_norm_map[tnorm]
            return resolved
        
        # Build assay_map for GENERIC format (maps assay_name to column name in projectMetadata)
        assay_to_column_map = {}
        if params.get('metadata_format') == 'GENERIC' and 'projectMetadata' in data:
            assay_name_row = data['projectMetadata'][data['projectMetadata']['term_name'] == 'assay_name']
            if not assay_name_row.empty:
                assay_name_series = assay_name_row.iloc[0]
                for col_name, assay_name_val in assay_name_series.items():
                    if col_name in ['term_name', 'project_level'] or pd.isna(assay_name_val):
                        continue
                    assay_to_column_map[assay_name_val] = col_name
                reporter.add_text(f"Built assay-to-column map for GENERIC format: {assay_to_column_map}")
        
        def get_meta_value_series(faire_term, source_df, assay_series):
            """Return per-assay value if present, otherwise project-level value (case-insensitive term match)."""
            term_mask = source_df['term_name'].astype(str).str.strip().str.lower() == str(faire_term).strip().lower()
            field_row_df = source_df[term_mask]
            if field_row_df.empty:
                return pd.Series(None, index=assay_series.index if isinstance(assay_series, pd.Series) else range(len(dna_derived_df_final)))
            
            field_row = field_row_df.iloc[0]
            project_level_val = field_row.get('project_level')
            
            # If no assay context, use project-level
            if not isinstance(assay_series, pd.Series) or ('assay_name' not in dna_derived_df_final.columns):
                return pd.Series(project_level_val if pd.notna(project_level_val) and str(project_level_val).strip() else None,
                                 index=(assay_series.index if isinstance(assay_series, pd.Series) else range(len(dna_derived_df_final))))
            
            # Prefer per-assay value if available; else project-level
            def map_assay_value(assay):
                if assay_to_column_map and assay in assay_to_column_map:
                    col_name = assay_to_column_map[assay]
                    val = field_row.get(col_name)
                    if pd.notna(val) and str(val).strip():
                        return val
                if assay in field_row.index:
                    val2 = field_row.get(assay)
                    if pd.notna(val2) and str(val2).strip():
                        return val2
                return project_level_val if pd.notna(project_level_val) and str(project_level_val).strip() else None
            
            return assay_series.map(map_assay_value)
        
        # Apply projectMetadata fields
        for dwc_term, mapping_info in dna_derived_map_config.items():
            if not isinstance(mapping_info, dict): continue
            
            source = mapping_info.get('source')
            faire_term = mapping_info.get('faire_term')
            
            if source == 'projectMetadata':
                # Only populate if the field isn't already present from pass-through
                if dwc_term not in dna_derived_df_final.columns or dna_derived_df_final[dwc_term].isna().all():
                    values = get_meta_value_series(faire_term, data['projectMetadata'], dna_derived_df_final['assay_name'])
                    dna_derived_df_final[dwc_term] = values
                    reporter.add_text(f"Populated '{dwc_term}' from projectMetadata (faire_term: '{faire_term}')")
                else:
                    reporter.add_text(f"Field '{dwc_term}' already populated from pass-through data")
            
            elif source == 'analysisMetadata':
                values_to_apply = pd.Series([None] * len(dna_derived_df_final), index=dna_derived_df_final.index)
                for run_name, run_config in params['datafiles'].items():
                    assay_name = run_config['assay_name']
                    analysis_df = data.get('analysis_data_by_assay', {}).get(assay_name, {}).get(run_name)
                    if analysis_df is not None:
                        # NOAA-specific robustness: resolve columns and match case/format-insensitively
                        def _resolve_ci(df, name):
                            lower_map = {str(c).lower(): c for c in df.columns}
                            return lower_map.get(str(name).lower())
                        term_col = _resolve_ci(analysis_df, 'term_name') or 'term_name'
                        val_col = _resolve_ci(analysis_df, 'values') or 'values'
                        if term_col not in analysis_df.columns or val_col not in analysis_df.columns:
                            reporter.add_warning(f"analysisMetadata missing expected columns for assay '{assay_name}'. Found columns: {list(analysis_df.columns)}")
                            continue
                        def _norm2(s):
                            return ''.join(ch for ch in str(s).lower().strip() if ch.isalnum())
                        # Try to match on faire_term first, then DwC term as fallback
                        target_keys = [_norm2(faire_term)]
                        if dwc_term and _norm2(dwc_term) != _norm2(faire_term):
                            target_keys.append(_norm2(dwc_term))
                        series_norm = analysis_df[term_col].astype(str).map(_norm2)
                        mask_term = series_norm.isin(target_keys)
                        field_row_df = analysis_df[mask_term]
                        if not field_row_df.empty:
                            value = field_row_df.iloc[0][val_col]
                            mask_assay = dna_derived_df_final['assay_name'].astype(str).str.strip() == str(assay_name).strip()
                            values_to_apply.loc[mask_assay] = value
                        else:
                            try:
                                sample_terms = analysis_df[term_col].astype(str).head(15).tolist()
                                reporter.add_warning(f"analysisMetadata term not found for assay '{assay_name}': looked for '{faire_term}' (or '{dwc_term}'). Sample terms: {sample_terms}")
                            except Exception:
                                reporter.add_warning(f"analysisMetadata term not found for assay '{assay_name}': '{faire_term}' (or '{dwc_term}')")
                dna_derived_df_final[dwc_term] = values_to_apply
                reporter.add_text(f"Populated '{dwc_term}' from analysisMetadata (faire_term: '{faire_term}')")
        
        # --- SPECIAL LOGIC: Construct pcr_primer_reference ---
        reporter.add_text("Constructing 'pcr_primer_reference' field...")
        fwd_series = get_meta_value_series('pcr_primer_reference_forward', data['projectMetadata'], dna_derived_df_final['assay_name'])
        rev_series = get_meta_value_series('pcr_primer_reference_reverse', data['projectMetadata'], dna_derived_df_final['assay_name'])
        
        def combine_references(fwd, rev):
            fwd_str = str(fwd).strip() if pd.notna(fwd) else ""
            rev_str = str(rev).strip() if pd.notna(rev) else ""
            if fwd_str and rev_str:
                return fwd_str if fwd_str == rev_str else f"{fwd_str}|{rev_str}"
            return fwd_str or rev_str or pd.NA
        
        dna_derived_df_final['pcr_primer_reference'] = [combine_references(f, r) for f, r in zip(fwd_series, rev_series)]
        
        # --- SPECIAL LOGIC: Concatenate otu_db and otu_db_custom ---
        reporter.add_text("Constructing 'otu_db' field from 'otu_db' and 'otu_db_custom'...")
        try:
            # First, ensure 'otu_db' is populated if it's in the mapper
            if 'otu_db' in dna_derived_df_final.columns:
                otu_db_series = dna_derived_df_final['otu_db']
            else:
                otu_db_series = pd.Series([None] * len(dna_derived_df_final), index=dna_derived_df_final.index)
            
            # Then, get the custom values, which could be in either metadata type
            custom_db_series = pd.Series([None] * len(dna_derived_df_final), index=dna_derived_df_final.index)
            if params.get('metadata_format') == 'GENERIC':
                custom_db_series = get_meta_value_series('otu_db_custom', data['projectMetadata'], dna_derived_df_final['assay_name'])
            else: # NOAA format
                for run_name, run_config in params['datafiles'].items():
                    assay_name = run_config['assay_name']
                    analysis_df = data.get('analysis_data_by_assay', {}).get(assay_name, {}).get(run_name)
                    if analysis_df is not None:
                        # Find the value for 'otu_db_custom'
                        term_mask = analysis_df['term_name'].astype(str).str.strip().str.lower() == 'otu_db_custom'
                        if term_mask.any():
                            value = analysis_df.loc[term_mask, 'values'].iloc[0]
                            # Apply this value to all rows matching the current assay
                            assay_mask = dna_derived_df_final['assay_name'] == assay_name
                            custom_db_series.loc[assay_mask] = value
            
            # Combine them
            def combine_db_fields(db, custom):
                db_str = str(db).strip() if pd.notna(db) and str(db).strip() not in ['nan', 'None'] else ""
                custom_str = str(custom).strip() if pd.notna(custom) and str(custom).strip() not in ['nan', 'None'] else ""
                if db_str and custom_str:
                    return f"{db_str};{custom_str}"
                return db_str or custom_str or pd.NA
            
            dna_derived_df_final['otu_db'] = [combine_db_fields(db, custom) for db, custom in zip(otu_db_series, custom_db_series)]
            reporter.add_success("Successfully combined 'otu_db' and 'otu_db_custom' fields.")
            
        except Exception as e:
            reporter.add_warning(f"Could not combine otu_db fields: {e}")
        
        # --- SPECIAL LOGIC: Construct adapters field ---
        reporter.add_text("Constructing 'adapters' field...")
        try:
            fwd_series = get_meta_value_series('adapter_forward', data['projectMetadata'], dna_derived_df_final['assay_name'])
            rev_series = get_meta_value_series('adapter_reverse', data['projectMetadata'], dna_derived_df_final['assay_name'])
            
            def combine_adapters(fwd, rev):
                fwd_str = str(fwd).strip().upper() if pd.notna(fwd) else ""
                rev_str = str(rev).strip().upper() if pd.notna(rev) else ""
                if fwd_str and rev_str:
                    return f"{fwd_str};{rev_str}"
                return fwd_str or rev_str or pd.NA
            
            dna_derived_df_final['adapters'] = [combine_adapters(f, r) for f, r in zip(fwd_series, rev_series)]
            reporter.add_success("Successfully constructed 'adapters' field.")
            
        except Exception as e:
            reporter.add_warning(f"Could not construct adapters field: {e}")
        
        # --- FINAL STEP: Select and Order Columns ---
        
        # Now apply the mapper to select and rename columns from the merged data
        reporter.add_text("Applying Darwin Core mappings...")
        final_df = dna_derived_df_final.copy()
        
        rename_map = {}
        for dwc_term, mapping_info in dna_derived_map_config.items():
            if not isinstance(mapping_info, dict): continue
            source = mapping_info.get('source')
            faire_term = mapping_info.get('faire_term')
            if source in ['sampleMetadata', 'experimentRunMetadata'] and faire_term in final_df.columns:
                if faire_term != dwc_term: # Avoid renaming if names are already the same
                    rename_map[faire_term] = dwc_term
        
        final_df.rename(columns=rename_map, inplace=True)
        
        # --- Ensure we have all required columns and order them correctly ---
        # The mapper is the AUTHORITATIVE source. Only include what's in the mapper, in that exact order.
        final_ordered_columns = []
        for col in DESIRED_DNA_DERIVED_COLUMNS:
            final_ordered_columns.append(col)
            # Add unit columns right after their parent, if they exist
            unit_col_map = {
                'concentration': 'concentrationUnit'
            }
            # Only auto-inject the unit column if it's NOT explicitly listed in the mapper
            if col in unit_col_map and unit_col_map[col] not in final_ordered_columns and unit_col_map[col] not in DESIRED_DNA_DERIVED_COLUMNS:
                final_ordered_columns.append(unit_col_map[col])
        
        # Ensure all columns exist in the dataframe, adding empty ones if needed.
        for col in final_ordered_columns:
            if col not in final_df.columns:
                final_df[col] = pd.NA
        
        # Before selecting, ensure unique columns to avoid duplicates in output
        final_df = final_df.loc[:, ~final_df.columns.duplicated(keep='first')]
        
        # Select ONLY the columns from the mapper, in the exact order specified
        # and filter out any that might not exist in the final dataframe
        final_columns_that_exist = [col for col in final_ordered_columns if col in final_df.columns]
        dna_derived_extension = final_df[final_columns_that_exist]
        
        # --- FINAL DE-DUPLICATION ---
        # Drop any fully duplicate rows that might have been created during merges
        initial_rows = len(dna_derived_extension)
        # Avoid SettingWithCopyWarning by assigning the result rather than inplace on a potential slice
        if 'occurrenceID' in dna_derived_extension.columns and not dna_derived_extension['occurrenceID'].isna().all():
            dna_derived_extension = dna_derived_extension.drop_duplicates(subset=['occurrenceID'], keep='first')
        else:
            dna_derived_extension = dna_derived_extension.drop_duplicates(keep='first')
        
        final_rows = len(dna_derived_extension)
        if initial_rows > final_rows:
            reporter.add_text(f"Removed {initial_rows - final_rows} duplicate rows after merges.")
        
        reporter.add_text(f"Final DNA derived extension shape: {dna_derived_extension.shape}")
        
        # Add a preview of the dataframe to the HTML report
        reporter.add_dataframe(dna_derived_extension, "DNA Derived Extension Preview", max_rows=10)
        
        # Save DNA derived extension to CSV file
        reporter.add_text("Saving DNA derived extension to CSV...")
        
        output_dir = params.get('output_dir', 'processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        
        output_path = os.path.join(output_dir, "dna_derived_extension.csv")
        dna_derived_extension.to_csv(output_path, index=False, encoding='utf-8-sig') # Encoding helps with special characters in units
        
        reporter.add_success("DNA derived extension created successfully")
        reporter.add_text(f"Saved DNA derived extension: {len(dna_derived_extension):,} records")
        reporter.add_text(f"Output file: dna_derived_extension.csv")
        
        # Verify the file was created
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024*1024)  # Size in MB
            reporter.add_text(f"File size: {file_size:.2f} MB")
            print(f"DNA derived extension created! Saved {len(dna_derived_extension):,} records")
        else:
            reporter.add_error("Error: File was not created")
        
    except Exception as e:
        reporter.add_error(f"DNA derived extension creation failed: {str(e)}")
        print(f"DNA derived extension creation failed: {str(e)}")
        import traceback
        traceback.print_exc() 