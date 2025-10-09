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
                    reporter.add_text(f"✓ Brought forward {len(available_passthrough)} pass-through field(s) from occurrence core.")
                else:
                    reporter.add_text("✓ No additional pass-through fields needed from occurrence core.")
            else:
                reporter.add_text("⚠️ No pass-through fields found in occurrence core.")
        except Exception as _e:
            reporter.add_text(f"⚠️ Could not merge pass-through fields from occurrence core: {_e}")

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
        
        # --- STEP W: Merge sampleMetadata and experimentRunMetadata directly to ensure fields exist ---
        try:
            # Collect sampleMetadata faire_terms from mapper (include potential _unit companions)
            sm_terms = []
            for _dwc, mapping in dna_derived_map_config.items():
                if isinstance(mapping, dict) and mapping.get('source') == 'sampleMetadata':
                    ft = mapping.get('faire_term')
                    if isinstance(ft, str) and ft.strip():
                        sm_terms.append(ft)
                        sm_terms.append(f"{ft}_unit")

            sm_resolved = _resolve_source_columns(data.get('sampleMetadata', pd.DataFrame()), sm_terms)
            sm_source_cols = sorted(set(['samp_name'] + list(sm_resolved.values()))) if sm_resolved else []
            if sm_source_cols:
                sm_merge_df = data['sampleMetadata'][sm_source_cols].copy()
                sm_merge_df['samp_name'] = sm_merge_df['samp_name'].astype(str).str.strip()
                dna_derived_df_final['parentEventID'] = dna_derived_df_final['parentEventID'].astype(str).str.strip()
                dna_derived_df_final = dna_derived_df_final.merge(sm_merge_df, left_on='parentEventID', right_on='samp_name', how='left')
                # Drop any duplicate samp_name columns
                if 'samp_name_y' in dna_derived_df_final.columns:
                    dna_derived_df_final.drop(columns=['samp_name_y'], inplace=True)
                if 'samp_name_x' in dna_derived_df_final.columns:
                    dna_derived_df_final.rename(columns={'samp_name_x': 'samp_name'}, inplace=True)
                # Alias resolved columns back to faire_term names so rename_map later can find them
                for faire_term, real_col in sm_resolved.items():
                    if faire_term not in dna_derived_df_final.columns and real_col in dna_derived_df_final.columns:
                        dna_derived_df_final[faire_term] = dna_derived_df_final[real_col]
                reporter.add_text(f"✓ Merged {len(sm_resolved)} sampleMetadata field(s) via parentEventID→samp_name (resolved & aliased).")
            else:
                reporter.add_text("⚠️ No sampleMetadata fields from mapper could be resolved in sampleMetadata sheet.")
        except Exception as _e:
            reporter.add_text(f"⚠️ Could not merge sampleMetadata directly: {_e}")

        # Fallback: derive samp_name from ERM by eventID and re-merge sample fields if still empty
        try:
            if sm_source_cols:
                # Check if most of the sample terms are still empty
                empty_ratio = []
                for t in sm_terms:
                    if isinstance(t, str) and t in dna_derived_df_final.columns:
                        empty_ratio.append(float(dna_derived_df_final[t].isna().mean()))
                needs_fb = empty_ratio and any(r > 0.95 for r in empty_ratio)
                if needs_fb and 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                    erm_map = data['experimentRunMetadata'][['lib_id', 'samp_name']].copy()
                    erm_map.rename(columns={'lib_id': 'eventID', 'samp_name': '_samp_from_erm'}, inplace=True)
                    erm_map['eventID'] = erm_map['eventID'].astype(str).str.strip()
                    erm_map['_samp_from_erm'] = erm_map['_samp_from_erm'].astype(str).str.strip()
                    dna_derived_df_final['eventID'] = dna_derived_df_final['eventID'].astype(str).str.strip()
                    dna_derived_df_final = dna_derived_df_final.merge(erm_map.drop_duplicates('eventID'), on='eventID', how='left')
                    # Re-merge sample fields using _samp_from_erm
                    sm_fb_df = data['sampleMetadata'][sm_source_cols].copy()
                    sm_fb_df.rename(columns={'samp_name': '_samp_from_erm'}, inplace=True)
                    sm_fb_df['_samp_from_erm'] = sm_fb_df['_samp_from_erm'].astype(str).str.strip()
                    dna_derived_df_final = dna_derived_df_final.merge(sm_fb_df, on='_samp_from_erm', how='left', suffixes=('', '_fb'))
                    # Fill from resolved and fallback columns into the target faire_term column
                    for faire_term, real_col in sm_resolved.items():
                        fb_col = f"{real_col}_fb"
                        # Ensure target column exists
                        if faire_term not in dna_derived_df_final.columns:
                            dna_derived_df_final[faire_term] = pd.NA
                        # Prefer existing target, then resolved col, then fallback column
                        if real_col in dna_derived_df_final.columns:
                            dna_derived_df_final[faire_term] = dna_derived_df_final[faire_term].combine_first(dna_derived_df_final[real_col])
                        if fb_col in dna_derived_df_final.columns:
                            dna_derived_df_final[faire_term] = dna_derived_df_final[faire_term].combine_first(dna_derived_df_final[fb_col])
                            dna_derived_df_final.drop(columns=[fb_col], inplace=True)
                    reporter.add_text("✓ Applied ERM→samp_name fallback for sample fields.")
        except Exception:
            pass

        # Merge experimentRunMetadata fields via resolved columns
        try:
            erm_terms = []
            for _dwc, mapping in dna_derived_map_config.items():
                if isinstance(mapping, dict) and mapping.get('source') == 'experimentRunMetadata':
                    ft = mapping.get('faire_term')
                    if isinstance(ft, str) and ft.strip() and ft != 'lib_id':
                        erm_terms.append(ft)
            # Ensure lib_id is resolvable
            erm_resolved = _resolve_source_columns(data.get('experimentRunMetadata', pd.DataFrame()), erm_terms + ['lib_id'])
            if erm_resolved and 'lib_id' in erm_resolved:
                erm_cols = sorted(set([erm_resolved['lib_id']] + [erm_resolved[t] for t in erm_terms if t in erm_resolved]))
                erm_merge_df = data['experimentRunMetadata'][erm_cols].copy()
                erm_merge_df.rename(columns={erm_resolved['lib_id']: 'eventID'}, inplace=True)
                erm_merge_df['eventID'] = erm_merge_df['eventID'].astype(str).str.strip()
                dna_derived_df_final['eventID'] = dna_derived_df_final['eventID'].astype(str).str.strip()
                dna_derived_df_final = dna_derived_df_final.merge(erm_merge_df, on='eventID', how='left')
                # Alias back to faire_term names
                for t in erm_terms:
                    if t in erm_resolved:
                        real_col = erm_resolved[t]
                        if t not in dna_derived_df_final.columns and real_col in dna_derived_df_final.columns:
                            dna_derived_df_final[t] = dna_derived_df_final[real_col]
                reporter.add_text(f"✓ Merged {len(erm_terms)} ERM field(s) via eventID with resolved columns.")
            else:
                reporter.add_text("⚠️ Could not resolve 'lib_id' in ERM for eventID join.")
        except Exception as _e:
            reporter.add_text(f"⚠️ Could not merge ERM fields: {_e}")

        # --- STEP V: FORCE-BACKFILL mapped fields by direct key lookup (robust) ---
        try:
            # Helper: resolve a column in a DF by case-insensitive / normalized name
            def _resolve_col(df: pd.DataFrame, name: str):
                if not isinstance(df, pd.DataFrame) or df.empty or not isinstance(name, str):
                    return None
                if name in df.columns:
                    return name
                lower_map = {str(c).lower(): c for c in df.columns}
                if name.lower() in lower_map:
                    return lower_map[name.lower()]
                def _n(s: str) -> str:
                    return ''.join(ch for ch in str(s).lower().strip().replace(' ', '').replace('\t','').replace('\n','') if ch.isalnum() or ch == '_')
                col_norm_map = {_n(c): c for c in df.columns}
                key = _n(name)
                return col_norm_map.get(key)

            # 1) Backfill from sampleMetadata via parentEventID → samp_name
            if 'sampleMetadata' in data and not data['sampleMetadata'].empty:
                sm_df = data['sampleMetadata'].copy()
                sm_df['samp_name'] = sm_df['samp_name'].astype(str).str.strip()
                par_series = dna_derived_df_final['parentEventID'].astype(str).str.strip()
                for dwc_term, mapping in dna_derived_map_config.items():
                    if not isinstance(mapping, dict):
                        continue
                    if mapping.get('source') != 'sampleMetadata':
                        continue
                    faire_term = mapping.get('faire_term')
                    if not isinstance(faire_term, str) or not faire_term.strip():
                        continue
                    real_col = _resolve_col(sm_df, faire_term)
                    if not real_col:
                        # try unit companion
                        real_col = _resolve_col(sm_df, f"{faire_term}_unit")
                        if not real_col:
                            continue
                    value_map = sm_df.set_index('samp_name')[real_col]
                    if faire_term not in dna_derived_df_final.columns:
                        dna_derived_df_final[faire_term] = pd.NA
                    dna_derived_df_final[faire_term] = dna_derived_df_final[faire_term].fillna(par_series.map(value_map))

                # If most sample terms are still empty (generic workflows), fallback via ERM → samp_name
                sample_terms = [m.get('faire_term') for m in dna_derived_map_config.values() if isinstance(m, dict) and m.get('source') == 'sampleMetadata']
                sample_terms = [t for t in sample_terms if isinstance(t, str)]
                if sample_terms:
                    empties = [float(dna_derived_df_final[t].isna().mean()) for t in sample_terms if t in dna_derived_df_final.columns]
                    # In GENERIC mode, always attempt ERM→samp_name supplement to be robust
                    if (params.get('metadata_format') == 'GENERIC' or (empties and any(r > 0.95 for r in empties))) and 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                        erm_df2 = data['experimentRunMetadata'].copy()
                        lib_col2 = _resolve_col(erm_df2, 'lib_id') or 'lib_id'
                        erm_df2[lib_col2] = erm_df2[lib_col2].astype(str).str.strip()
                        erm_df2['_samp_from_erm'] = erm_df2['samp_name'].astype(str).str.strip() if 'samp_name' in erm_df2.columns else pd.NA
                        evt_series2 = dna_derived_df_final['eventID'].astype(str).str.strip()
                        lib_to_samp = erm_df2.set_index(lib_col2)['_samp_from_erm'] if '_samp_from_erm' in erm_df2.columns else None
                        if lib_to_samp is not None:
                            samp_from_evt = evt_series2.map(lib_to_samp)
                            # Map each term from sampleMetadata using _samp_from_erm
                            for t in sample_terms:
                                real_col2 = _resolve_col(sm_df, t) or _resolve_col(sm_df, f"{t}_unit")
                                if real_col2 is None or t not in dna_derived_df_final.columns:
                                    continue
                                val_map2 = sm_df.set_index('samp_name')[real_col2]
                                dna_derived_df_final[t] = dna_derived_df_final[t].fillna(samp_from_evt.map(val_map2))

            # 2) Backfill from ERM via eventID ↔ lib_id
            if 'experimentRunMetadata' in data and not data['experimentRunMetadata'].empty:
                erm_df = data['experimentRunMetadata'].copy()
                lib_col = _resolve_col(erm_df, 'lib_id') or 'lib_id'
                erm_df[lib_col] = erm_df[lib_col].astype(str).str.strip()
                evt_series = dna_derived_df_final['eventID'].astype(str).str.strip()
                for dwc_term, mapping in dna_derived_map_config.items():
                    if not isinstance(mapping, dict):
                        continue
                    if mapping.get('source') != 'experimentRunMetadata':
                        continue
                    faire_term = mapping.get('faire_term')
                    if not isinstance(faire_term, str) or not faire_term.strip() or faire_term == 'lib_id':
                        continue
                    real_col = _resolve_col(erm_df, faire_term)
                    if not real_col:
                        continue
                    value_map = erm_df.set_index(lib_col)[real_col]
                    if faire_term not in dna_derived_df_final.columns:
                        dna_derived_df_final[faire_term] = pd.NA
                    dna_derived_df_final[faire_term] = dna_derived_df_final[faire_term].fillna(evt_series.map(value_map))

            # 3) GENERIC fallback: map via materialSampleID if present (source_mat_id → sampleMetadata.materialSampleID)
            try:
                if params.get('metadata_format') == 'GENERIC' and 'sampleMetadata' in data and not data['sampleMetadata'].empty:
                    sm_df2 = data['sampleMetadata'].copy()
                    mat_col = _resolve_col(sm_df2, 'materialSampleID') or 'materialSampleID'
                    if mat_col in sm_df2.columns and 'source_mat_id' in dna_derived_df_final.columns:
                        sm_df2[mat_col] = sm_df2[mat_col].astype(str).str.strip()
                        src_series = dna_derived_df_final['source_mat_id'].astype(str).str.strip()
                        for dwc_term, mapping in dna_derived_map_config.items():
                            if not isinstance(mapping, dict):
                                continue
                            if mapping.get('source') != 'sampleMetadata':
                                continue
                            faire_term = mapping.get('faire_term')
                            if not isinstance(faire_term, str) or not faire_term.strip():
                                continue
                            real_col = _resolve_col(sm_df2, faire_term) or _resolve_col(sm_df2, f"{faire_term}_unit")
                            if not real_col:
                                continue
                            value_map3 = sm_df2.set_index(mat_col)[real_col]
                            if faire_term not in dna_derived_df_final.columns:
                                dna_derived_df_final[faire_term] = pd.NA
                            dna_derived_df_final[faire_term] = dna_derived_df_final[faire_term].fillna(src_series.map(value_map3))
            except Exception:
                pass

            reporter.add_text("✓ Force-backfilled sample/ERM fields by direct key mapping.")
        except Exception as _e:
            reporter.add_text(f"⚠️ Force-backfill step skipped due to error: {_e}")

        # --- STEP X: Handle Units (Safe and non-destructive) ---
        reporter.add_text("Applying specific unit logic...")

        # Helper to set a unit column only when missing or empty
        def _set_unit_if_missing(df, unit_col, default_value):
            if unit_col not in df.columns:
                df[unit_col] = default_value
            else:
                # fill only NAs/empty strings
                df[unit_col] = df[unit_col].where(~(df[unit_col].isna() | (df[unit_col].astype(str).str.strip() == '')), other=default_value)

        # 1) size_fracUnit: default 'micrometer' if absent/empty
        _set_unit_if_missing(dna_derived_df_final, 'size_fracUnit', 'micrometer')
        reporter.add_text("✓ Ensured 'size_fracUnit' present (default 'micrometer' if missing).")

        # 2) concentrationUnit: prefer mapped 'concentration_unit' if present, else default 'ng/µl'
        # First, check if we have a source concentration_unit column
        if 'concentration_unit' in dna_derived_df_final.columns:
            # Use the source column to populate concentrationUnit
            dna_derived_df_final['concentrationUnit'] = dna_derived_df_final['concentration_unit'].fillna('ng/\u00b5l')
            # Drop the source column to avoid duplicates
            dna_derived_df_final.drop(columns=['concentration_unit'], inplace=True)
            reporter.add_text("✓ Populated 'concentrationUnit' from source 'concentration_unit' column.")
        else:
            # No source column, set default
            _set_unit_if_missing(dna_derived_df_final, 'concentrationUnit', 'ng/\u00b5l')
            reporter.add_text("✓ Set 'concentrationUnit' to default 'ng/µl'.")

        # 3) samp_vol_we_dna_extUnit: prefer mapped unit column if present, else default 'mL'
        if 'samp_vol_we_dna_ext_unit' in dna_derived_df_final.columns:
            dna_derived_df_final['samp_vol_we_dna_extUnit'] = dna_derived_df_final['samp_vol_we_dna_ext_unit'].fillna('mL')
            # normalize casing ml -> mL
            dna_derived_df_final['samp_vol_we_dna_extUnit'] = dna_derived_df_final['samp_vol_we_dna_extUnit'].astype(str).str.replace('ml', 'mL', case=False)
            # Drop source column to avoid duplicates
            dna_derived_df_final.drop(columns=['samp_vol_we_dna_ext_unit'], inplace=True)
            reporter.add_text("✓ Populated 'samp_vol_we_dna_extUnit' from source 'samp_vol_we_dna_ext_unit' column.")
        else:
            _set_unit_if_missing(dna_derived_df_final, 'samp_vol_we_dna_extUnit', 'mL')
            reporter.add_text("✓ Set 'samp_vol_we_dna_extUnit' to default 'mL'.")

        
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
                    reporter.add_text(f"✓ Populated '{dwc_term}' from projectMetadata (faire_term: '{faire_term}')")
                else:
                    reporter.add_text(f"✓ Field '{dwc_term}' already populated from pass-through data")
            
            elif source == 'analysisMetadata':
                values_to_apply = pd.Series([None] * len(dna_derived_df_final), index=dna_derived_df_final.index)
                for run_name, run_config in params['datafiles'].items():
                    assay_name = run_config['assay_name']
                    analysis_df = data.get('analysis_data_by_assay', {}).get(assay_name, {}).get(run_name)
                    if analysis_df is not None:
                        field_row = analysis_df[analysis_df['term_name'] == faire_term]
                        if not field_row.empty:
                            value = field_row.iloc[0]['values']
                            mask = dna_derived_df_final['assay_name'] == assay_name
                            values_to_apply.loc[mask] = value
                dna_derived_df_final[dwc_term] = values_to_apply
                reporter.add_text(f"✓ Populated '{dwc_term}' from analysisMetadata (faire_term: '{faire_term}')")

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
                'size_frac': 'size_fracUnit',
                'concentration': 'concentrationUnit',
                'samp_vol_we_dna_ext': 'samp_vol_we_dna_extUnit'
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
            print(f"✅ DNA derived extension created! Saved {len(dna_derived_extension):,} records")
        else:
            reporter.add_error("❌ Error: File was not created")
        
    except Exception as e:
        reporter.add_error(f"DNA derived extension creation failed: {str(e)}")
        print(f"❌ DNA derived extension creation failed: {str(e)}")
        import traceback
        traceback.print_exc() 