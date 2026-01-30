"""
EML (Ecological Metadata Language) Builder for edna2obis

Generates EML XML files conforming to the GBIF Metadata Profile (version 1.1) 
from configuration parameters. EML files are required for OBIS/GBIF dataset submission.

Key Features:
- Generates valid XML conforming to GBIF EML profile
- Supports all major EML sections: people, coverage, keywords, methods, etc.
- Auto-generates UUIDs for document identification
- Validates required fields and provides helpful error messages
- Outputs to processed-v3/ directory automatically
"""

import os
import uuid
import yaml
import xml.etree.ElementTree as ET
from xml.dom import minidom
from datetime import datetime
from typing import Dict, List, Optional, Any
import pandas as pd


def _unique_clean_values_from_series(series: pd.Series) -> List[str]:
    """Return sorted unique non-empty strings from a Series, removing NaN/None/empty and NaN-like tokens."""
    try:
        if series is None or series.empty:
            return []
        cleaned = series.dropna().astype(str).map(lambda s: s.strip())
        invalid_tokens = {"", "nan", "none", "null", "na", "n/a"}
        cleaned = cleaned[~cleaned.str.lower().isin(invalid_tokens)]
        return sorted(pd.unique(cleaned))
    except Exception:
        return []


def _get_project_value(project_df: pd.DataFrame, field_name: str) -> Optional[str]:
    """Fetch a value for a given field from projectMetadata regardless of layout.

    Supports two layouts:
      1) Columnar layout: a column literally named `field_name` exists; returns first non-empty value.
      2) Term-list layout: a row where term_name == field_name; prefer 'project_level' column if present,
         otherwise use the first non-empty value among remaining columns besides 'term_name'.
    """
    if not isinstance(project_df, pd.DataFrame) or project_df.empty or not field_name:
        return None

    # Case 1: direct column
    if field_name in project_df.columns:
        try:
            series = project_df[field_name]
            vals = _unique_clean_values_from_series(series)
            if vals:
                return str(vals[0])
        except Exception:
            pass

    # Case 2: term_name layout
    cols_lower = {str(c).strip().lower(): c for c in project_df.columns}
    term_col = cols_lower.get('term_name')
    if term_col and term_col in project_df.columns:
        try:
            row = project_df[project_df[term_col].astype(str).str.strip().str.lower() == str(field_name).strip().lower()]
            if not row.empty:
                row0 = row.iloc[0]
                # Prefer project_level if present
                if 'project_level' in project_df.columns:
                    val = row0.get('project_level')
                    val = None if pd.isna(val) else str(val).strip()
                    if val:
                        return val
                # Else, scan other columns for the first non-empty value
                for col in project_df.columns:
                    if str(col).strip().lower() in {'term_name'}:
                        continue
                    val = row0.get(col)
                    val = None if pd.isna(val) else str(val).strip()
                    if val:
                        return val
        except Exception:
            pass

    return None


def _auto_populate_coverage(eml_config: Dict, occurrence_df: pd.DataFrame, reporter) -> Dict:
    """Auto-populates coverage data from the occurrence DataFrame based on config comments."""
    if occurrence_df.empty:
        reporter.add_text("Occurrence data not available, skipping auto-population of coverage.")
        return eml_config

    reporter.add_text("Attempting to auto-populate EML coverage fields from data...")

    # Geographic Coverage - Bounding Box
    geo_cfg = eml_config.get('geographic_coverage', {})
    lat_field = str(geo_cfg.get('latitude_field', 'decimalLatitude'))
    lon_field = str(geo_cfg.get('longitude_field', 'decimalLongitude'))

    if lat_field in occurrence_df.columns and lon_field in occurrence_df.columns:
        latitudes = pd.to_numeric(occurrence_df[lat_field], errors='coerce').dropna()
        longitudes = pd.to_numeric(occurrence_df[lon_field], errors='coerce').dropna()

        if not latitudes.empty and not longitudes.empty:
            eml_config.setdefault('geographic_coverage', {}).setdefault('bounding_coordinates', {})
            bounds = eml_config['geographic_coverage']['bounding_coordinates']
            if not bounds.get('west'): bounds['west'] = longitudes.min()
            if not bounds.get('east'): bounds['east'] = longitudes.max()
            if not bounds.get('north'): bounds['north'] = latitudes.max()
            if not bounds.get('south'): bounds['south'] = latitudes.min()
            reporter.add_success("Calculated geographic bounding coordinates.")

    # Geographic Coverage - Description (configurable source fields)
    geo_desc_field = 'description'
    eml_config.setdefault('geographic_coverage', {})
    current_desc = eml_config['geographic_coverage'].get(geo_desc_field, "")
    if not current_desc or 'Geographic area description' in current_desc:
        source_fields = geo_cfg.get('description_source_fields', ['locality'])
        values = []
        for field in source_fields:
            if field in occurrence_df.columns:
                vals = [str(v).strip() for v in occurrence_df[field].dropna().unique() if str(v).strip()]
                values.extend(vals)
        if values:
            unique_vals = sorted(pd.unique(pd.Series(values)))
            location_str = ', '.join(unique_vals)
            eml_config['geographic_coverage'][geo_desc_field] = location_str
            reporter.add_success("Populated geographic description from configured source field(s).")
            # Update Title with location
            current_title = eml_config.get('title', '')
            if '[LOCATION]' in current_title:
                eml_config['title'] = current_title.replace('[LOCATION]', location_str)
                reporter.add_success(f"Updated EML title with location: {location_str}")
        else:
            reporter.add_warning("Could not populate geographic description: none of the configured fields found or non-empty.")

    # Temporal Coverage
    event_date_field = str(eml_config.get('temporal_coverage', {}).get('event_date_field', 'eventDate'))
    if event_date_field in occurrence_df.columns:
        dates = pd.to_datetime(occurrence_df[event_date_field], errors='coerce').dropna()
        if not dates.empty:
            eml_config.setdefault('temporal_coverage', {})
            temporal = eml_config['temporal_coverage']
            # Always override with calculated dates from data
            temporal['begin_date'] = dates.min().strftime('%Y-%m-%d')
            temporal['end_date'] = dates.max().strftime('%Y-%m-%d')
            reporter.add_success("Calculated temporal coverage dates from data.")
            # Update Title with timeframe
            current_title = eml_config.get('title', '')
            if '[TIMEFRAME]' in current_title:
                begin_year = dates.min().year
                end_year = dates.max().year
                timeframe = str(begin_year) if begin_year == end_year else f"{begin_year}-{end_year}"
                eml_config['title'] = current_title.replace('[TIMEFRAME]', timeframe)
                reporter.add_success(f"Updated EML title with timeframe: {timeframe}")
    return eml_config

def _auto_populate_taxonomic_coverage(eml_config: Dict, occurrence_df: pd.DataFrame, reporter) -> Dict:
    """Auto-populates taxonomic coverage from the occurrence DataFrame using multiple higher ranks (OBIS/GBIF)."""
    if occurrence_df.empty:
        return eml_config

    eml_config.setdefault('taxonomic_coverage', {})
    tax_coverage = eml_config['taxonomic_coverage']

    # Populate classifications with higher-level ranks only, across multiple ranks
    if not tax_coverage.get('classifications'):
        # Corrected default rank priority
        rank_priority = tax_coverage.get('rank_fields', ['kingdom', 'phylum', 'class', 'order', 'family', 'genus'])
        classifications = []
        rank_counts = {}
        per_rank_limit = 10  # Limit per rank to avoid one rank dominating
        total_limit = 10  # Set total limit to 10 as requested
        total_truncated = False

        for rank in rank_priority:
            if rank not in occurrence_df.columns:
                continue

            values = _unique_clean_values_from_series(occurrence_df[rank])

            truncated_this_rank = False
            if len(values) > per_rank_limit:
                values = values[:per_rank_limit]
                truncated_this_rank = True

            # Respect total limit
            remaining = max(0, total_limit - len(classifications))
            if remaining == 0:
                total_truncated = True
                break
            if len(values) > remaining:
                values = values[:remaining]
                total_truncated = True

            for v in values:
                classifications.append({'taxon_rank_name': rank, 'taxon_rank_value': str(v)})

            rank_counts[rank] = len(values)

            if truncated_this_rank:
                reporter.add_text(f"Note: Truncated {rank} list to {per_rank_limit} unique values for EML brevity.")

        if classifications:
            tax_coverage['classifications'] = classifications
            parts = [f"{r}={c}" for r, c in rank_counts.items() if c > 0]
            reporter.add_success(f"Populated taxonomic coverage with {len(classifications)} classifications across ranks: " + ", ".join(parts))
            if total_truncated:
                reporter.add_text(f"Note: Reached total cap of {total_limit} taxonomic classifications.")
        else:
            reporter.add_warning("Could not populate taxonomic coverage: no higher-rank columns (phylum/class/order/kingdom/domain) found.")

    return eml_config

def _auto_populate_from_project_data(eml_config: Dict, project_df: pd.DataFrame, reporter) -> Dict:
    """Auto-populates EML fields from project metadata (FAIRe projectMetadata sheet)."""
    if project_df.empty: return eml_config

    # Auto-populate project title if it's empty or template
    current_proj_title = eml_config.get('project', {}).get('title', '')
    if not current_proj_title or "Your Project Title" in current_proj_title:
        title_src_field = str(eml_config.get('project', {}).get('title_source_field', 'project_name'))
        # Try robust fetch across columnar and term_name layouts
        project_name = _get_project_value(project_df, title_src_field)
        if project_name:
            eml_config.setdefault('project', {})['title'] = project_name
            reporter.add_success(f"Populated EML project title with: {project_name}")
    return eml_config


def _safe_get(data: Dict, path: str, default: str = "") -> str:
    """Safely get nested dictionary values with dot notation"""
    keys = path.split('.')
    current = data
    
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    
    return str(current) if current is not None else default


def _create_person_element(parent: ET.Element, tag_name: str, person_data: Dict) -> None:
    """Create a person/organization element in EML format"""
    if not person_data:
        return
        
    person_elem = ET.SubElement(parent, tag_name)
    
    # Individual name
    given_name = _safe_get(person_data, 'individual_name.given_name')
    surname = _safe_get(person_data, 'individual_name.surname')
    
    if given_name or surname:
        individual_name = ET.SubElement(person_elem, 'individualName')
        if given_name:
            ET.SubElement(individual_name, 'givenName').text = given_name
        if surname:
            ET.SubElement(individual_name, 'surName').text = surname
    
    # Organization name
    org_name = _safe_get(person_data, 'organization_name')
    if org_name:
        ET.SubElement(person_elem, 'organizationName').text = org_name
    
    # Position name (alternative to individual name)
    position = _safe_get(person_data, 'position_name')
    if position:
        ET.SubElement(person_elem, 'positionName').text = position
    
    # Address
    address_data = person_data.get('address', {})
    if any(address_data.values()):
        address = ET.SubElement(person_elem, 'address')
        
        delivery_point = _safe_get(address_data, 'delivery_point')
        if delivery_point:
            ET.SubElement(address, 'deliveryPoint').text = delivery_point
            
        city = _safe_get(address_data, 'city')
        if city:
            ET.SubElement(address, 'city').text = city
            
        admin_area = _safe_get(address_data, 'administrative_area')
        if admin_area:
            ET.SubElement(address, 'administrativeArea').text = admin_area
            
        postal_code = _safe_get(address_data, 'postal_code')
        if postal_code:
            ET.SubElement(address, 'postalCode').text = postal_code
            
        country = _safe_get(address_data, 'country')
        if country:
            ET.SubElement(address, 'country').text = country
    
    # Contact info
    phone = _safe_get(person_data, 'phone')
    if phone:
        ET.SubElement(person_elem, 'phone').text = phone
        
    email = _safe_get(person_data, 'email')
    if email:
        ET.SubElement(person_elem, 'electronicMailAddress').text = email
        
    online_url = _safe_get(person_data, 'online_url')
    if online_url:
        ET.SubElement(person_elem, 'onlineUrl').text = online_url
        
    # User ID (e.g., ORCID)
    user_id = _safe_get(person_data, 'user_id')
    if user_id:
        user_id_elem = ET.SubElement(person_elem, 'userId')
        user_id_elem.text = user_id
        # Set directory attribute for ORCID
        if 'orcid.org' in user_id.lower():
            user_id_elem.set('directory', 'https://orcid.org')


def _create_coverage_element(parent: ET.Element, eml_config: Dict) -> None:
    """Create the coverage element with geographic, temporal, and taxonomic coverage"""
    coverage = ET.SubElement(parent, 'coverage')
    
    # Geographic coverage
    geo_config = eml_config.get('geographic_coverage', {})
    if geo_config:
        geo_coverage = ET.SubElement(coverage, 'geographicCoverage')
        
        description = _safe_get(geo_config, 'description')
        if description:
            ET.SubElement(geo_coverage, 'geographicDescription').text = description
        
        bounds_config = geo_config.get('bounding_coordinates', {})
        if bounds_config:
            bounding_coords = ET.SubElement(geo_coverage, 'boundingCoordinates')
            
            west = _safe_get(bounds_config, 'west')
            if west:
                ET.SubElement(bounding_coords, 'westBoundingCoordinate').text = str(west)
                
            east = _safe_get(bounds_config, 'east')
            if east:
                ET.SubElement(bounding_coords, 'eastBoundingCoordinate').text = str(east)
                
            north = _safe_get(bounds_config, 'north')
            if north:
                ET.SubElement(bounding_coords, 'northBoundingCoordinate').text = str(north)
                
            south = _safe_get(bounds_config, 'south')
            if south:
                ET.SubElement(bounding_coords, 'southBoundingCoordinate').text = str(south)
    
    # Temporal coverage
    temporal_config = eml_config.get('temporal_coverage', {})
    if temporal_config:
        temporal_coverage = ET.SubElement(coverage, 'temporalCoverage')
        
        begin_date = _safe_get(temporal_config, 'begin_date')
        end_date = _safe_get(temporal_config, 'end_date')
        
        if begin_date and end_date and begin_date != end_date:
            # Date range
            range_of_dates = ET.SubElement(temporal_coverage, 'rangeOfDates')
            begin_elem = ET.SubElement(range_of_dates, 'beginDate')
            ET.SubElement(begin_elem, 'calendarDate').text = begin_date
            end_elem = ET.SubElement(range_of_dates, 'endDate')
            ET.SubElement(end_elem, 'calendarDate').text = end_date
        elif begin_date:
            # Single date
            single_date = ET.SubElement(temporal_coverage, 'singleDateTime')
            ET.SubElement(single_date, 'calendarDate').text = begin_date
    
    # Taxonomic coverage
    taxonomic_config = eml_config.get('taxonomic_coverage', {})
    if taxonomic_config:
        taxonomic_coverage = ET.SubElement(coverage, 'taxonomicCoverage')
        
        general_desc = _safe_get(taxonomic_config, 'general_description')
        if general_desc:
            ET.SubElement(taxonomic_coverage, 'generalTaxonomicCoverage').text = general_desc
        
        # Taxonomic classifications
        classifications = taxonomic_config.get('classifications', [])
        for classification in classifications:
            if isinstance(classification, dict):
                tax_class = ET.SubElement(taxonomic_coverage, 'taxonomicClassification')
                
                rank_name = _safe_get(classification, 'taxon_rank_name')
                if rank_name:
                    ET.SubElement(tax_class, 'taxonRankName').text = rank_name
                    
                rank_value = _safe_get(classification, 'taxon_rank_value')
                if rank_value:
                    ET.SubElement(tax_class, 'taxonRankValue').text = rank_value
                    
                common_name = _safe_get(classification, 'common_name')
                if common_name:
                    ET.SubElement(tax_class, 'commonName').text = common_name


def _create_keywords_element(parent: ET.Element, keywords: List[Dict]) -> None:
    """Create keyword sets"""
    if not keywords:
        return
        
    for keyword_data in keywords:
        if isinstance(keyword_data, dict):
            keyword_set = ET.SubElement(parent, 'keywordSet')
            
            keyword = _safe_get(keyword_data, 'keyword')
            if keyword:
                ET.SubElement(keyword_set, 'keyword').text = keyword
                
            thesaurus = _safe_get(keyword_data, 'thesaurus')
            if thesaurus:
                ET.SubElement(keyword_set, 'keywordThesaurus').text = thesaurus
            else:
                # Required field, use placeholder if not provided
                ET.SubElement(keyword_set, 'keywordThesaurus').text = "N/A"


def _create_methods_element(parent: ET.Element, methods_config: Dict) -> None:
    """Create methods element with sampling and quality control information"""
    if not methods_config:
        return
        
    methods = ET.SubElement(parent, 'methods')
    
    # Method steps
    method_steps = methods_config.get('method_steps', [])
    for step_data in method_steps:
        if isinstance(step_data, dict):
            method_step = ET.SubElement(methods, 'methodStep')
            description = _safe_get(step_data, 'description')
            if description:
                desc_elem = ET.SubElement(method_step, 'description')
                ET.SubElement(desc_elem, 'para').text = description
    
    # Quality control
    quality_control = _safe_get(methods_config, 'quality_control')
    if quality_control:
        qc_elem = ET.SubElement(methods, 'qualityControl')
        qc_desc = ET.SubElement(qc_elem, 'description')
        ET.SubElement(qc_desc, 'para').text = quality_control
    
    # Sampling
    study_extent = _safe_get(methods_config, 'study_extent')
    sampling_desc = _safe_get(methods_config, 'sampling_description')
    
    if study_extent or sampling_desc:
        sampling = ET.SubElement(methods, 'sampling')
        
        if study_extent:
            study_extent_elem = ET.SubElement(sampling, 'studyExtent')
            ET.SubElement(study_extent_elem, 'description').text = study_extent
            
        if sampling_desc:
            sampling_desc_elem = ET.SubElement(sampling, 'samplingDescription')
            ET.SubElement(sampling_desc_elem, 'para').text = sampling_desc


def _create_project_element(parent: ET.Element, project_config: Dict) -> None:
    """Create project element with project information"""
    if not project_config:
        return
        
    project = ET.SubElement(parent, 'project')
    
    # Project title
    title = _safe_get(project_config, 'title')
    if title:
        ET.SubElement(project, 'title').text = title
    
    # Project personnel
    personnel_list = project_config.get('personnel', [])
    for person_data in personnel_list:
        if isinstance(person_data, dict):
            personnel = ET.SubElement(project, 'personnel')
            
            # Individual name
            given_name = _safe_get(person_data, 'individual_name.given_name')
            surname = _safe_get(person_data, 'individual_name.surname')
            
            if given_name or surname:
                individual_name = ET.SubElement(personnel, 'individualName')
                if given_name:
                    ET.SubElement(individual_name, 'givenName').text = given_name
                if surname:
                    ET.SubElement(individual_name, 'surName').text = surname
            
            # Role
            role = _safe_get(person_data, 'role')
            if role:
                ET.SubElement(personnel, 'role').text = role
    
    # Project description
    description = _safe_get(project_config, 'description')
    if description:
        desc_elem = ET.SubElement(project, 'description')
        ET.SubElement(desc_elem, 'para').text = description
    
    # Funding
    funding = _safe_get(project_config, 'funding')
    if funding:
        funding_elem = ET.SubElement(project, 'funding')
        ET.SubElement(funding_elem, 'para').text = funding


def _get_license_text(license_type: str) -> str:
    """Get the appropriate license text based on license type"""
    license_texts = {
        'CC0': """This work is licensed under a Creative Commons CCZero 1.0 License http://creativecommons.org/publicdomain/zero/1.0/. 
To the extent possible under law, the publisher has waived all copyright and related or neighboring rights to this dataset and dedicated it to the public domain worldwide. 
Users may copy, modify, distribute and use the dataset, including for commercial purposes, without restriction.""",
        
        'CC-BY': """This work is licensed under a Creative Commons Attribution 4.0 License http://creativecommons.org/licenses/by/4.0/legalcode. 
You are free to share, copy and redistribute the material in any medium or format, and adapt, remix, transform, and build upon the material 
for any purpose, even commercially. You must give appropriate credit, provide a link to the license, and indicate if changes were made. 
You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.""",
        
        'CC-BY-NC': """This work is licensed under a Creative Commons Attribution Non Commercial 4.0 License http://creativecommons.org/licenses/by-nc/4.0/legalcode. 
You are free to share, copy and redistribute the material in any medium or format, and adapt, remix, transform, and build upon the material. 
You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may not use the material for commercial purposes. 
You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use."""
    }
    
    return license_texts.get(license_type, license_texts['CC0'])


def _prettify_xml(elem: ET.Element) -> str:
    """Return a pretty-printed XML string for the Element"""
    rough_string = ET.tostring(elem, 'unicode')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def _validate_required_fields(eml_config: Dict, reporter) -> bool:
    """Validate that required EML fields are present"""
    required_fields = [
        ('title', 'Dataset title is required'),
        ('abstract', 'Dataset abstract is required'),
        ('contact.email', 'Contact email is required'),
    ]
    
    missing_fields = []
    
    for field_path, error_msg in required_fields:
        value = _safe_get(eml_config, field_path)
        if not value or value.strip() == "":
            missing_fields.append(error_msg)
    
    if missing_fields:
        reporter.add_error("Missing required EML fields:")
        for error in missing_fields:
            reporter.add_error(f"  - {error}")
        return False
    
    return True


def _validate_recommended_fields(eml_config: Dict, reporter) -> None:
    """Check for recommended EML fields and provide warnings with guidance"""
    
    # Check for recommended fields with helpful guidance
    recommendations = [
        {
            'field': 'creator.individual_name.given_name',
            'message': 'Creator given name is recommended for proper dataset citation',
            'guidance': 'Add creator -> individual_name -> given_name in config.yaml'
        },
        {
            'field': 'creator.individual_name.surname', 
            'message': 'Creator surname is recommended for proper dataset citation',
            'guidance': 'Add creator -> individual_name -> surname in config.yaml'
        },
        {
            'field': 'creator.organization_name',
            'message': 'Creator organization name is recommended',
            'guidance': 'Add creator -> organization_name in config.yaml'
        },
        {
            'field': 'geographic_coverage.description',
            'message': 'Geographic description improves dataset discoverability',
            'guidance': 'Add geographic_coverage -> description (e.g., "Gulf of Mexico", "North Atlantic")'
        },
        {
            'field': 'geographic_coverage.bounding_coordinates.west',
            'message': 'Geographic bounding coordinates help users understand spatial coverage',
            'guidance': 'Add geographic_coverage -> bounding_coordinates with west, east, north, south values'
        },
        {
            'field': 'temporal_coverage.begin_date',
            'message': 'Temporal coverage helps users understand when data was collected',
            'guidance': 'Add temporal_coverage -> begin_date and end_date (YYYY-MM-DD format)'
        },
        {
            'field': 'keywords',
            'message': 'Keywords improve dataset discoverability in OBIS/GBIF',
            'guidance': 'Add keywords list with relevant terms like "eDNA", "metabarcoding", "marine biodiversity"'
        },
        {
            'field': 'methods.sampling_description',
            'message': 'Sampling methods description helps users understand data collection procedures',
            'guidance': 'Add methods -> sampling_description with details about sample collection'
        },
        {
            'field': 'purpose',
            'message': 'Dataset purpose helps users understand research objectives',
            'guidance': 'Add purpose field describing why this dataset was created'
        }
    ]
    
    missing_recommendations = []
    
    for rec in recommendations:
        value = _safe_get(eml_config, rec['field'])
        
        # Special handling for different field types
        if rec['field'] == 'keywords':
            keywords = eml_config.get('keywords', [])
            if not keywords or len(keywords) == 0:
                missing_recommendations.append(rec)
        elif rec['field'] == 'geographic_coverage.bounding_coordinates.west':
            # Check if any bounding coordinate is present
            bounds = eml_config.get('geographic_coverage', {}).get('bounding_coordinates', {})
            if not any(bounds.values()):
                missing_recommendations.append(rec)
        else:
            if not value or value.strip() == "":
                missing_recommendations.append(rec)
    
    # Report recommendations
    if missing_recommendations:
        reporter.add_warning("EML Recommendations - Consider adding these fields for better metadata quality:")
        for rec in missing_recommendations:
            reporter.add_warning(f"• {rec['message']}")
            reporter.add_text(f"  Guidance: {rec['guidance']}")
    else:
        reporter.add_success("EML metadata appears comprehensive - all recommended fields are present!")
    
    # Check for template placeholders that need updating
    template_warnings = []
    title = _safe_get(eml_config, 'title')
    if '[LOCATION]' in title or '[TIMEFRAME]' in title:
        template_warnings.append("Title contains placeholder text [LOCATION] or [TIMEFRAME] - please update with specific information")
    
    abstract = _safe_get(eml_config, 'abstract')
    if 'Brief description of your dataset' in abstract:
        template_warnings.append("Abstract appears to contain template text - please customize with your specific dataset description")
    
    contact_email = _safe_get(eml_config, 'contact.email')
    if 'your.email@' in contact_email or '@institution.edu' in contact_email:
        template_warnings.append("Contact email appears to be template text - please update with actual email address")
    
    if template_warnings:
        reporter.add_warning("Template placeholders detected - please customize:")
        for warning in template_warnings:
            reporter.add_warning(f"• {warning}")
    
    # Provide usage guidance
    reporter.add_text("<h4>EML Best Practices:</h4>")
    reporter.add_list([
        "Use descriptive titles that include taxonomic, geographic, and temporal scope",
        "Provide detailed abstracts that help users assess fitness-for-use", 
        "Include complete contact information for data inquiries",
        "Add relevant keywords to improve discoverability",
        "Document sampling methods similar to a journal article methods section",
        "Use CC0 license when possible for maximum data reusability"
    ])


def create_eml_file(params: Dict, data: Dict[str, pd.DataFrame], reporter) -> str:
    """
    Create an EML XML file from configuration parameters.
    
    Args:
        params: Configuration parameters from config.yaml
        data: Dictionary containing project metadata DataFrames
        reporter: HTML reporter for logging
    
    Returns:
        Path to the generated EML file
    """
    reporter.add_section("Creating EML (Ecological Metadata Language) File")
    
    try:
        # Load EML configuration from separate file
        eml_config_path = params.get('eml_config_path', 'EML_config.yaml')
        
        if not os.path.exists(eml_config_path):
            reporter.add_error(f"EML configuration file not found: {eml_config_path}")
            raise FileNotFoundError(f"EML config file not found: {eml_config_path}")
        
        reporter.add_text(f"Loading EML configuration from: {eml_config_path}")
        
        try:
            with open(eml_config_path, 'r', encoding='utf-8') as f:
                eml_full_config = yaml.safe_load(f)
        except Exception as e:
            reporter.add_error(f"Failed to parse EML configuration file: {e}")
            raise ValueError(f"Invalid EML config file: {e}")
        
        # Extract the eml_metadata section
        eml_config = eml_full_config.get('eml_metadata', {})
        if not eml_config:
            reporter.add_error("No 'eml_metadata' section found in EML_config.yaml")
            raise ValueError("EML metadata configuration missing from EML_config.yaml")
        
        # Validate required fields
        if not _validate_required_fields(eml_config, reporter):
            raise ValueError("Required EML fields are missing")
        
        # Check recommended fields and provide guidance
        _validate_recommended_fields(eml_config, reporter)
        
        # Auto-populate EML fields from data, if available
        occurrence_df = data.get('occurrence', pd.DataFrame())
        if occurrence_df.empty:
            try:
                output_dir = params.get('output_dir', 'processed-v3/')
                api_choice = str(params.get('taxonomic_api_source', 'GBIF')).lower()
                occ_filename = f"occurrence_core_{api_choice}.csv"
                occ_path = os.path.join(output_dir, occ_filename)
                if os.path.exists(occ_path):
                    occurrence_df = pd.read_csv(occ_path)
                    reporter.add_text(f"Using occurrence data from file: {occ_path}")
                else:
                    reporter.add_warning(f"Occurrence file not found at expected path: {occ_path}. Skipping auto-population from data.")
            except Exception as e:
                reporter.add_warning(f"Could not load occurrence data for EML auto-population: {e}")
        project_df = data.get('projectMetadata', pd.DataFrame())

        # Dataset-level title from eml_metadata.title_source_field (if configured)
        dataset_title_src = eml_full_config.get('eml_metadata', {}).get('title_source_field')
        if dataset_title_src:
            current_title = str(eml_config.get('title', ''))
            needs_fill = (not current_title) or ('[LOCATION]' in current_title) or ('[TIMEFRAME]' in current_title)
            if needs_fill:
                val = _get_project_value(project_df, str(dataset_title_src))
                if val:
                    eml_config['title'] = val
                    reporter.add_success(f"Set EML dataset title from projectMetadata.{dataset_title_src}")

        eml_config = _auto_populate_coverage(eml_config, occurrence_df, reporter)
        eml_config = _auto_populate_taxonomic_coverage(eml_config, occurrence_df, reporter)
        eml_config = _auto_populate_from_project_data(eml_config, project_df, reporter)
        
        reporter.add_text("Building EML XML structure...")
        
        # Report what will be included in the EML
        sections_included = []
        if _safe_get(eml_config, 'creator.individual_name.given_name') or _safe_get(eml_config, 'creator.organization_name'):
            sections_included.append("Creator information")
        if _safe_get(eml_config, 'contact.email'):
            sections_included.append("Contact information")
        if eml_config.get('geographic_coverage'):
            sections_included.append("Geographic coverage")
        if eml_config.get('temporal_coverage'):
            sections_included.append("Temporal coverage") 
        if eml_config.get('taxonomic_coverage'):
            sections_included.append("Taxonomic coverage")
        if eml_config.get('keywords'):
            sections_included.append(f"Keywords ({len(eml_config.get('keywords', []))} terms)")
        if eml_config.get('methods'):
            sections_included.append("Methods and sampling")
        if eml_config.get('project'):
            sections_included.append("Project information")
        
        if sections_included:
            reporter.add_text("EML sections to be included:")
            reporter.add_list(sections_included)
        
        # Generate unique package ID
        package_id = f"{uuid.uuid4()}/eml-1.xml"
        
        # Create root EML element with proper namespaces
        root = ET.Element('eml:eml')
        root.set('packageId', package_id)
        root.set('system', 'https://doi.org')
        root.set('xmlns:eml', 'https://eml.ecoinformatics.org/eml-2.2.0')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        root.set('xmlns:stmml', 'http://www.xml-cml.org/schema/stmml-1.1')
        root.set('xsi:schemaLocation', 'https://eml.ecoinformatics.org/eml-2.2.0 http://rs.gbif.org/schema/eml-gbif-profile/1.1/eml.xsd')
        
        # Create dataset element
        dataset = ET.SubElement(root, 'dataset')
        
        # Basic metadata
        title = _safe_get(eml_config, 'title')
        if title:
            ET.SubElement(dataset, 'title').text = title
        
        # Creator (required)
        creator_config = eml_config.get('creator', {})
        if creator_config:
            _create_person_element(dataset, 'creator', creator_config)
        
        # Metadata provider (optional)
        metadata_provider_config = eml_config.get('metadata_provider', {})
        if metadata_provider_config and any(metadata_provider_config.values()):
            _create_person_element(dataset, 'metadataProvider', metadata_provider_config)
        
        # Associated parties
        associated_parties = eml_config.get('associated_parties', [])
        for party in associated_parties:
            if isinstance(party, dict):
                assoc_party = ET.SubElement(dataset, 'associatedParty')
                
                # Add person details directly to associatedParty element
                given_name = _safe_get(party, 'individual_name.given_name')
                surname = _safe_get(party, 'individual_name.surname')
                
                if given_name or surname:
                    individual_name = ET.SubElement(assoc_party, 'individualName')
                    if given_name:
                        ET.SubElement(individual_name, 'givenName').text = given_name
                    if surname:
                        ET.SubElement(individual_name, 'surName').text = surname
                
                org_name = _safe_get(party, 'organization_name')
                if org_name:
                    ET.SubElement(assoc_party, 'organizationName').text = org_name
                
                email = _safe_get(party, 'email')
                if email:
                    ET.SubElement(assoc_party, 'electronicMailAddress').text = email
                
                role = _safe_get(party, 'role')
                if role:
                    ET.SubElement(assoc_party, 'role').text = role
        
        # Contact (required)
        contact_config = eml_config.get('contact', {})
        if contact_config:
            _create_person_element(dataset, 'contact', contact_config)
        
        # Publication date
        pub_date = _safe_get(eml_config, 'publication_date')
        if pub_date:
            ET.SubElement(dataset, 'pubDate').text = pub_date
        
        # Language
        language = _safe_get(eml_config, 'language', 'en')
        ET.SubElement(dataset, 'language').text = language
        
        # Abstract (required)
        abstract = _safe_get(eml_config, 'abstract')
        if abstract:
            abstract_elem = ET.SubElement(dataset, 'abstract')
            ET.SubElement(abstract_elem, 'para').text = abstract.strip()
        
        # Keywords
        keywords = eml_config.get('keywords', [])
        _create_keywords_element(dataset, keywords)
        
        # Coverage (geographic, temporal, taxonomic)
        _create_coverage_element(dataset, eml_config)
        
        # Purpose
        purpose = _safe_get(eml_config, 'purpose')
        if purpose:
            purpose_elem = ET.SubElement(dataset, 'purpose')
            ET.SubElement(purpose_elem, 'para').text = purpose.strip()
        
        # Intellectual rights (license)
        license_type = _safe_get(eml_config, 'license', 'CC0')
        license_text = _get_license_text(license_type)
        intellectual_rights = ET.SubElement(dataset, 'intellectualRights')
        ET.SubElement(intellectual_rights, 'para').text = license_text
        
        # Methods
        methods_config = eml_config.get('methods', {})
        if methods_config:
            _create_methods_element(dataset, methods_config)
        
        # Project
        project_config = eml_config.get('project', {})
        if project_config:
            _create_project_element(dataset, project_config)
        
        # Additional metadata
        additional_info = _safe_get(eml_config, 'additional_info')
        if additional_info:
            ET.SubElement(dataset, 'additionalInfo').text = additional_info.strip()
        
        # Generate pretty XML
        xml_string = _prettify_xml(root)
        
        # Save to file
        output_dir = params.get('output_dir', 'processed-v3/')
        os.makedirs(output_dir, exist_ok=True)
        eml_path = os.path.join(output_dir, 'eml.xml')
        
        with open(eml_path, 'w', encoding='utf-8') as f:
            f.write(xml_string)
        
        # Get file size for reporting
        file_size_kb = os.path.getsize(eml_path) / 1024
        
        reporter.add_success(f"EML file successfully created: {eml_path}")
        reporter.add_text(f"Package ID: {package_id}")
        reporter.add_text(f"File size: {file_size_kb:.1f} KB")
        
        # Count XML elements for reporting
        try:
            tree = ET.parse(eml_path)
            root = tree.getroot()
            element_count = len(list(root.iter()))
            reporter.add_text(f"XML structure: {element_count} total elements")
        except:
            pass
        
        # Validate the XML structure
        try:
            ET.parse(eml_path)
            reporter.add_success("EML XML validation passed - file is well-formed and ready for OBIS/GBIF submission")
        except ET.ParseError as e:
            reporter.add_warning(f"EML XML validation warning: {e}")
        
        # Provide next steps guidance
        reporter.add_text("<h4>Next Steps for OBIS/GBIF Submission:</h4>")
        reporter.add_list([
            "Review the generated EML file to ensure all information is accurate",
            "Upload the EML file along with your Darwin Core data files to IPT",
            "Or submit directly to your OBIS node manager for review",
            "The EML file conforms to GBIF Metadata Profile v1.1 standards"
        ])
        
        return eml_path
        
    except Exception as e:
        error_msg = f"Failed to create EML file: {str(e)}"
        reporter.add_error(error_msg)
        raise Exception(error_msg)
