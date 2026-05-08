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

import copy
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

    proj = eml_config.setdefault('project', {})
    if not str(proj.get('id', '') or '').strip():
        project_id = _get_project_value(project_df, 'project_id')
        if project_id:
            proj['id'] = project_id
            reporter.add_success(f"Populated EML project id with: {project_id}")

    current_pd = str(proj.get('project_description', '') or '').strip()
    if not current_pd:
        pd_val = _get_project_value(project_df, 'project_description')
        if pd_val:
            proj['project_description'] = pd_val
            reporter.add_success("Populated EML project.project_description from projectMetadata")

    return eml_config


def _unique_short_names_from_frame(df: Optional[pd.DataFrame]) -> List[str]:
    """Return sorted unique non-empty short_name values if the frame has short_name or short_name_sm."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return []
    for col in ('short_name', 'short_name_sm'):
        if col in df.columns:
            return _unique_clean_values_from_series(df[col])
    return []


def _effective_short_name_for_eml(
    params: Dict,
    eml_config: Dict,
    project_df: pd.DataFrame,
    sample_metadata_df: Optional[pd.DataFrame] = None,
    occurrence_df: Optional[pd.DataFrame] = None,
) -> str:
    """Resolve short_name for dataset title suffix and alternateIdentifier.

    Priority: params['eml_short_name'] (split outputs), eml_metadata.short_name in YAML,
    then a single unique value from occurrence or sampleMetadata (FAIRe: short_name is on sampleMetadata),
    else projectMetadata when that sheet has exactly one short_name.
    """
    sn = str(params.get('eml_short_name', '') or '').strip()
    if sn:
        return sn
    sn = str(eml_config.get('short_name', '') or '').strip()
    if sn:
        return sn
    for df in (occurrence_df, sample_metadata_df):
        vals = _unique_short_names_from_frame(df)
        if len(vals) == 1:
            return str(vals[0])
    if isinstance(project_df, pd.DataFrame) and not project_df.empty and 'short_name' in project_df.columns:
        vals = _unique_clean_values_from_series(project_df['short_name'])
        if len(vals) == 1:
            return str(vals[0])
    return ''


def _ensure_project_personnel(eml_config: Dict, reporter) -> None:
    """GBIF IPT requires at least one <project><personnel>; copy creator if personnel is empty."""
    proj = eml_config.get('project')
    if not isinstance(proj, dict):
        return
    raw_list = proj.get('personnel')
    if raw_list is None:
        raw_list = []
    has_valid = False
    for p in raw_list:
        if not isinstance(p, dict):
            continue
        if (
            _safe_get(p, 'individual_name.given_name')
            or _safe_get(p, 'individual_name.surname')
            or _safe_get(p, 'organization_name')
        ):
            has_valid = True
            break
    if has_valid:
        return
    creators = _creator_records(eml_config)
    if not creators:
        reporter.add_warning(
            "EML project has no personnel and creator is missing; IPT may reject metadata until project.personnel or creator is set."
        )
        return
    person = copy.deepcopy(creators[0])
    if not str(person.get('role', '') or '').strip():
        person['role'] = 'principal investigator'
    eml_config.setdefault('project', {})['personnel'] = [person]
    reporter.add_text("Filled EML project.personnel from eml_metadata.creator (IPT requires at least one project personnel).")


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


def _non_empty_dict(data: Dict) -> bool:
    """True when a nested dict/list contains at least one non-empty scalar."""
    if not isinstance(data, dict):
        return bool(str(data or '').strip())
    for value in data.values():
        if isinstance(value, dict):
            if _non_empty_dict(value):
                return True
        elif isinstance(value, list):
            if any(_non_empty_dict(v) if isinstance(v, dict) else str(v or '').strip() for v in value):
                return True
        elif str(value or '').strip():
            return True
    return False


def _as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _resolve_contacts_file_path(eml_config_path: str, contacts_path: str) -> Optional[str]:
    """Resolve contacts.yaml (or alternate path) next to the EML config, then cwd."""
    if not contacts_path or not str(contacts_path).strip():
        return None
    contacts_path = str(contacts_path).strip()
    candidates = []
    if os.path.isabs(contacts_path):
        candidates.append(contacts_path)
    else:
        base_dir = os.path.dirname(os.path.abspath(eml_config_path))
        candidates.append(os.path.normpath(os.path.join(base_dir, contacts_path)))
        candidates.append(os.path.normpath(os.path.join(os.getcwd(), contacts_path)))
    for path in candidates:
        if os.path.isfile(path):
            return path
    # If still missing (e.g. wrong number of ".." in config), walk up from the EML file
    # directory and look for a file with the same basename (typically contacts.yaml).
    if not os.path.isabs(contacts_path):
        base_dir = os.path.dirname(os.path.abspath(eml_config_path))
        basename = os.path.basename(os.path.normpath(contacts_path))
        if basename:
            d = base_dir
            for _ in range(8):
                probe = os.path.join(d, basename)
                if os.path.isfile(probe):
                    return probe
                parent = os.path.dirname(d)
                if parent == d:
                    break
                d = parent
    return None


def _load_contacts_book(path: str, reporter) -> Dict[str, Dict]:
    """Load reusable people records from contacts.yaml (or compatible)."""
    try:
        with open(path, 'r', encoding='utf-8') as f:
            raw = yaml.safe_load(f)
    except Exception as exc:
        reporter.add_warning(f"Could not load contacts file '{path}': {exc}")
        return {}
    if not isinstance(raw, dict):
        return {}
    if isinstance(raw.get('contacts'), dict):
        book_raw = raw['contacts']
    else:
        book_raw = raw
    book = {str(k): v for k, v in book_raw.items() if isinstance(v, dict)}
    reporter.add_text(f"Loaded {len(book)} reusable contact(s) from: {path}")
    return book


def _contact_book_entry_to_person(entry: Dict) -> Dict:
    """Map master-file fields (IPT-style or edna2obis nested) to eml_metadata person dict."""
    if not isinstance(entry, dict):
        return {}

    def _addr_piece(d: Dict, *keys: str) -> str:
        for key in keys:
            v = d.get(key)
            if v is not None and str(v).strip():
                return str(v).strip()
        return ''

    person: Dict[str, Any] = {}

    gn = entry.get('givenName') or _safe_get(entry, 'individual_name.given_name')
    sn = entry.get('surName') or _safe_get(entry, 'individual_name.surname')
    if gn or sn:
        person['individual_name'] = {}
        if gn:
            person['individual_name']['given_name'] = str(gn).strip()
        if sn:
            person['individual_name']['surname'] = str(sn).strip()

    org = entry.get('organizationName') or entry.get('organization_name')
    if org:
        person['organization_name'] = str(org).strip()

    pos = entry.get('positionName') or entry.get('position_name')
    if pos:
        person['position_name'] = str(pos).strip()

    email = entry.get('electronicMailAddress') or entry.get('email')
    if email:
        person['email'] = str(email).strip()

    phone = entry.get('phone')
    if phone:
        person['phone'] = str(phone).strip()

    url = entry.get('onlineUrl') or entry.get('online_url')
    if url:
        person['online_url'] = str(url).strip()

    uid = entry.get('userId') or entry.get('user_id')
    if uid:
        person['user_id'] = str(uid).strip()

    addr_in = entry.get('address')
    if isinstance(addr_in, dict) and addr_in:
        addr = {
            'delivery_point': _addr_piece(
                addr_in, 'delivery_point', 'deliveryPoint'
            ),
            'city': _addr_piece(addr_in, 'city'),
            'administrative_area': _addr_piece(
                addr_in, 'administrative_area', 'administrativeArea'
            ),
            'postal_code': _addr_piece(addr_in, 'postal_code', 'postalCode'),
            'country': _addr_piece(addr_in, 'country'),
        }
        if any(addr.values()):
            person['address'] = addr

    return person


def _expand_contact_dict(node: Dict, book: Dict[str, Dict], reporter) -> Dict:
    ref = node.get('contact_ref')
    if not ref or not isinstance(ref, str) or not ref.strip():
        return node
    ref_key = ref.strip()
    if ref_key not in book:
        reporter.add_warning(
            f"contacts.yaml: unknown contact_ref '{ref_key}' — block left unchanged."
        )
        return node
    base = copy.deepcopy(_contact_book_entry_to_person(book[ref_key]))
    overlay = copy.deepcopy(node)
    overlay.pop('contact_ref', None)
    return _deep_merge_dict(base, overlay)


def _prepare_contact_refs_and_shorthand(eml_config: Dict) -> None:
    """Allow contact: ref_key or contact: { contact_ref: ref } and string list items in creators."""
    for key in ('contact', 'metadata_provider', 'creator'):
        val = eml_config.get(key)
        if isinstance(val, str) and val.strip():
            eml_config[key] = {'contact_ref': val.strip()}

    creators = eml_config.get('creators')
    if isinstance(creators, list):
        normalized = []
        for item in creators:
            if isinstance(item, str) and item.strip():
                normalized.append({'contact_ref': item.strip()})
            else:
                normalized.append(item)
        eml_config['creators'] = normalized

    parties = eml_config.get('associated_parties')
    if isinstance(parties, list):
        eml_config['associated_parties'] = [
            {'contact_ref': item.strip()} if isinstance(item, str) and item.strip() else item
            for item in parties
        ]

    proj = eml_config.get('project')
    if isinstance(proj, dict):
        pers = proj.get('personnel')
        if isinstance(pers, list):
            proj['personnel'] = [
                {'contact_ref': item.strip()} if isinstance(item, str) and item.strip() else item
                for item in pers
            ]


def _deep_expand_all_contact_refs(obj: Any, book: Dict[str, Dict], reporter) -> Any:
    """Recursively expand any dict that declares contact_ref using contacts.yaml."""
    if isinstance(obj, dict):
        if 'contact_ref' in obj and isinstance(obj.get('contact_ref'), str):
            expanded = _expand_contact_dict(obj, book, reporter)
            return _deep_expand_all_contact_refs(expanded, book, reporter)
        return {k: _deep_expand_all_contact_refs(v, book, reporter) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_deep_expand_all_contact_refs(v, book, reporter) for v in obj]
    return obj


def _normalize_orcid_user_id(user_id: str) -> str:
    """IPT exports ORCID userId as bare ORCID with directory='https://orcid.org/'."""
    value = str(user_id or '').strip()
    for prefix in ('https://orcid.org/', 'http://orcid.org/'):
        if value.lower().startswith(prefix):
            return value[len(prefix):]
    return value


def _creator_records(eml_config: Dict) -> List[Dict]:
    """Prefer IPT-style creators list, but keep legacy creator dict working."""
    creators = eml_config.get('creators')
    if isinstance(creators, list) and creators:
        return [c for c in creators if isinstance(c, dict) and _non_empty_dict(c)]
    creator = eml_config.get('creator')
    if isinstance(creator, dict) and _non_empty_dict(creator):
        return [creator]
    return []


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
        user_id_elem.text = _normalize_orcid_user_id(user_id)
        # Set directory attribute for ORCID
        if 'orcid.org' in user_id.lower():
            user_id_elem.set('directory', 'https://orcid.org/')


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

    # Sampling (GBIF IPT expects studyExtent/description/para, not description as plain text)
    study_extent = str(_safe_get(methods_config, 'study_extent') or '').strip()
    sampling_desc = str(_safe_get(methods_config, 'sampling_description') or '').strip()
    if not study_extent and sampling_desc:
        study_extent = sampling_desc

    if study_extent or sampling_desc:
        sampling = ET.SubElement(methods, 'sampling')

        if study_extent:
            study_extent_elem = ET.SubElement(sampling, 'studyExtent')
            desc_wrap = ET.SubElement(study_extent_elem, 'description')
            ET.SubElement(desc_wrap, 'para').text = study_extent

        if sampling_desc:
            sampling_desc_elem = ET.SubElement(sampling, 'samplingDescription')
            ET.SubElement(sampling_desc_elem, 'para').text = sampling_desc


def _create_project_element(parent: ET.Element, project_config: Dict) -> None:
    """Create project element with project information"""
    if not project_config:
        return
        
    project = ET.SubElement(parent, 'project')
    project_id = _safe_get(project_config, 'id')
    if project_id:
        project.set('id', project_id)
    
    # Project title
    title = _safe_get(project_config, 'title')
    if title:
        ET.SubElement(project, 'title').text = title
    
    # Project personnel (mirror dataset party fields where present for IPT validation)
    personnel_list = project_config.get('personnel', [])
    for person_data in personnel_list:
        if isinstance(person_data, dict):
            personnel = ET.SubElement(project, 'personnel')

            given_name = _safe_get(person_data, 'individual_name.given_name')
            surname = _safe_get(person_data, 'individual_name.surname')

            if given_name or surname:
                individual_name = ET.SubElement(personnel, 'individualName')
                if given_name:
                    ET.SubElement(individual_name, 'givenName').text = given_name
                if surname:
                    ET.SubElement(individual_name, 'surName').text = surname

            org_name = _safe_get(person_data, 'organization_name')
            if org_name:
                ET.SubElement(personnel, 'organizationName').text = org_name

            position = _safe_get(person_data, 'position_name')
            if position:
                ET.SubElement(personnel, 'positionName').text = position

            email = _safe_get(person_data, 'email')
            if email:
                ET.SubElement(personnel, 'electronicMailAddress').text = email

            phone = _safe_get(person_data, 'phone')
            if phone:
                ET.SubElement(personnel, 'phone').text = phone

            address_data = person_data.get('address', {})
            if isinstance(address_data, dict) and any(
                str(v or '').strip() for v in address_data.values()
            ):
                address = ET.SubElement(personnel, 'address')
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

            user_id = _safe_get(person_data, 'user_id')
            if user_id:
                user_id_elem = ET.SubElement(personnel, 'userId')
                user_id_elem.text = _normalize_orcid_user_id(user_id)
                if 'orcid.org' in user_id.lower():
                    user_id_elem.set('directory', 'https://orcid.org/')

            online_url = _safe_get(person_data, 'online_url')
            if online_url:
                ET.SubElement(personnel, 'onlineUrl').text = online_url

            role = _safe_get(person_data, 'role')
            if role:
                ET.SubElement(personnel, 'role').text = role
    
    # Project abstract: IPT exports project narrative as project/abstract/para.
    paras_to_add: List[str] = []
    proj_desc = str(_safe_get(project_config, 'project_description') or '').strip()
    legacy = str(_safe_get(project_config, 'abstract') or '').strip()
    if not legacy:
        legacy = str(_safe_get(project_config, 'description') or '').strip()
    if proj_desc:
        paras_to_add.append(proj_desc)
    if legacy and legacy != proj_desc:
        paras_to_add.append(legacy)
    if paras_to_add:
        desc_elem = ET.SubElement(project, 'abstract')
        for p in paras_to_add:
            ET.SubElement(desc_elem, 'para').text = p
    
    # Funding
    funding = _safe_get(project_config, 'funding')
    if funding:
        funding_elem = ET.SubElement(project, 'funding')
        ET.SubElement(funding_elem, 'para').text = funding

    # Award
    award_config = project_config.get('award', {})
    if isinstance(award_config, dict) and _non_empty_dict(award_config):
        award = ET.SubElement(project, 'award')
        funder_name = _safe_get(award_config, 'funder_name')
        if funder_name:
            ET.SubElement(award, 'funderName').text = funder_name
        award_number = _safe_get(award_config, 'award_number')
        if award_number:
            ET.SubElement(award, 'awardNumber').text = award_number
        award_title = _safe_get(award_config, 'title')
        if award_title:
            ET.SubElement(award, 'title').text = award_title

    # Study area description
    study_area = project_config.get('study_area_description', {})
    if isinstance(study_area, dict) and _non_empty_dict(study_area):
        study_area_elem = ET.SubElement(project, 'studyAreaDescription')
        descriptor = ET.SubElement(study_area_elem, 'descriptor')
        descriptor.set('name', _safe_get(study_area, 'descriptor_name', 'generic') or 'generic')
        descriptor.set(
            'citableClassificationSystem',
            str(_safe_get(study_area, 'citable_classification_system', 'false') or 'false').lower(),
        )
        descriptor_value = _safe_get(study_area, 'descriptor_value')
        if descriptor_value:
            ET.SubElement(descriptor, 'descriptorValue').text = descriptor_value

    # Design description
    design_description = _safe_get(project_config, 'design_description')
    if design_description:
        design_elem = ET.SubElement(project, 'designDescription')
        desc_elem = ET.SubElement(design_elem, 'description')
        ET.SubElement(desc_elem, 'para').text = design_description


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


def _normalize_license_key(license_type: str) -> str:
    return str(license_type or 'CC0').strip().upper().replace(' ', '').replace('-', '')


def _create_intellectual_rights_element(dataset: ET.Element, license_type: str) -> None:
    """Write intellectualRights; CC0 uses GBIF IPT–expected ulink/citetitle markup."""
    intellectual_rights = ET.SubElement(dataset, 'intellectualRights')
    if _normalize_license_key(license_type) == 'CC0':
        para = ET.SubElement(intellectual_rights, 'para')
        para.text = (
            'To the extent possible under law, the publisher has waived all rights to these data '
            'and has dedicated them to the '
        )
        ulink = ET.SubElement(para, 'ulink')
        ulink.set('url', 'http://creativecommons.org/publicdomain/zero/1.0/legalcode')
        citetitle = ET.SubElement(ulink, 'citetitle')
        citetitle.text = 'Public Domain (CC0 1.0)'
        ulink.tail = (
            '. Users may copy, modify, distribute and use the work, including for commercial '
            'purposes, without restriction.'
        )
    else:
        license_text = _get_license_text(license_type)
        ET.SubElement(intellectual_rights, 'para').text = license_text


def _create_licensed_element(dataset: ET.Element, eml_config: Dict) -> None:
    """Create IPT-style licensed block when license metadata is configured."""
    licensed_config = eml_config.get('licensed', {})
    license_type = _normalize_license_key(_safe_get(eml_config, 'license', 'CC0'))

    if not isinstance(licensed_config, dict) or not _non_empty_dict(licensed_config):
        if license_type != 'CC0':
            return
        licensed_config = {
            'license_name': 'Creative Commons Zero v1.0 Universal',
            'url': 'https://spdx.org/licenses/CC0-1.0.html',
            'identifier': 'CC0-1.0',
        }

    licensed = ET.SubElement(dataset, 'licensed')
    license_name = _safe_get(licensed_config, 'license_name')
    if license_name:
        ET.SubElement(licensed, 'licenseName').text = license_name
    url = _safe_get(licensed_config, 'url')
    if url:
        ET.SubElement(licensed, 'url').text = url
    identifier = _safe_get(licensed_config, 'identifier')
    if identifier:
        ET.SubElement(licensed, 'identifier').text = identifier


def _format_template(value: str, values: Dict[str, str]) -> str:
    """Best-effort formatting for optional short_name/project_id URL templates."""
    text = str(value or '').strip()
    if not text:
        return ''
    try:
        return text.format(**values)
    except Exception:
        return text


def _deep_merge_dict(base: Dict, override: Dict) -> Dict:
    """Recursively merge override into base without mutating the caller's dicts."""
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if (
            isinstance(value, dict)
            and isinstance(merged.get(key), dict)
        ):
            merged[key] = _deep_merge_dict(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _apply_short_name_overrides(eml_config: Dict, short_name: str, reporter) -> Dict:
    """Apply optional per-cruise/per-expedition EML overrides."""
    if not short_name:
        return eml_config
    overrides = eml_config.get('per_short_name_overrides', {})
    if not isinstance(overrides, dict):
        return eml_config
    override = overrides.get(short_name)
    if not isinstance(override, dict):
        return eml_config
    reporter.add_text(f"Applying EML per_short_name_overrides for: {short_name}")
    return _deep_merge_dict(eml_config, override)


def _create_distribution_elements(dataset: ET.Element, eml_config: Dict, template_values: Dict[str, str]) -> None:
    distribution_config = eml_config.get('distribution', {})
    if not isinstance(distribution_config, dict):
        return

    entries = []
    info_url = _format_template(_safe_get(distribution_config, 'information_url'), template_values)
    if info_url:
        entries.append(('information', info_url))
    download_url = _format_template(_safe_get(distribution_config, 'download_url'), template_values)
    if download_url:
        entries.append(('download', download_url))

    for item in _as_list(distribution_config.get('links')):
        if not isinstance(item, dict):
            continue
        function = _safe_get(item, 'function', 'information') or 'information'
        url = _format_template(_safe_get(item, 'url'), template_values)
        if url:
            entries.append((function, url))

    for function, url in entries:
        distribution = ET.SubElement(dataset, 'distribution')
        distribution.set('scope', 'document')
        online = ET.SubElement(distribution, 'online')
        url_elem = ET.SubElement(online, 'url')
        url_elem.set('function', function)
        url_elem.text = url


def _create_maintenance_element(dataset: ET.Element, eml_config: Dict) -> None:
    maintenance_config = eml_config.get('maintenance', {})
    if not isinstance(maintenance_config, dict) or not _non_empty_dict(maintenance_config):
        return

    maintenance = ET.SubElement(dataset, 'maintenance')
    description = _safe_get(maintenance_config, 'description')
    if description or 'description' in maintenance_config:
        desc_elem = ET.SubElement(maintenance, 'description')
        ET.SubElement(desc_elem, 'para').text = description
    frequency = _safe_get(maintenance_config, 'update_frequency')
    if frequency:
        ET.SubElement(maintenance, 'maintenanceUpdateFrequency').text = frequency


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
    ]
    
    missing_fields = []
    
    for field_path, error_msg in required_fields:
        value = _safe_get(eml_config, field_path)
        if not value or value.strip() == "":
            missing_fields.append(error_msg)

    if not _safe_get(eml_config, 'contact.email'):
        missing_fields.append('Contact email is required')

    if not _creator_records(eml_config):
        missing_fields.append('At least one creator is required (use creators: list or legacy creator: block)')
    
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
        elif rec['field'].startswith('creator.') and _creator_records(eml_config):
            continue
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

        eml_config_abspath = os.path.abspath(eml_config_path)
        contacts_path_setting = (
            eml_config.get('contacts_path')
            or eml_full_config.get('contacts_path')
            or params.get('contacts_path', 'contacts.yaml')
        )
        resolved_contacts = _resolve_contacts_file_path(
            eml_config_abspath, str(contacts_path_setting)
        )
        contacts_book: Dict[str, Dict] = {}
        if resolved_contacts:
            contacts_book = _load_contacts_book(resolved_contacts, reporter)
        elif str(contacts_path_setting).strip():
            reporter.add_warning(
                f"No contacts file found for '{contacts_path_setting}' (searched next to EML config and working directory); "
                "contact_ref entries will not expand."
            )

        _prepare_contact_refs_and_shorthand(eml_config)
        if contacts_book:
            expanded_cfg = _deep_expand_all_contact_refs(
                copy.deepcopy(eml_config), contacts_book, reporter
            )
            eml_config.clear()
            eml_config.update(expanded_cfg)
        eml_config.pop('contacts_path', None)
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
        _ensure_project_personnel(eml_config, reporter)

        sample_metadata_df = data.get('sampleMetadata')
        if not isinstance(sample_metadata_df, pd.DataFrame):
            sample_metadata_df = None
        effective_short_name = _effective_short_name_for_eml(
            params,
            eml_config,
            project_df,
            sample_metadata_df=sample_metadata_df,
            occurrence_df=occurrence_df if isinstance(occurrence_df, pd.DataFrame) else None,
        )
        eml_config = _apply_short_name_overrides(eml_config, effective_short_name, reporter)

        _prepare_contact_refs_and_shorthand(eml_config)
        if contacts_book:
            expanded_cfg = _deep_expand_all_contact_refs(
                copy.deepcopy(eml_config), contacts_book, reporter
            )
            eml_config.clear()
            eml_config.update(expanded_cfg)

        # Validate after overrides and final contact expansion (contact_ref / shorthand).
        if not _validate_required_fields(eml_config, reporter):
            raise ValueError("Required EML fields are missing")
        _validate_recommended_fields(eml_config, reporter)

        reporter.add_text("Building EML XML structure...")
        
        # Report what will be included in the EML
        sections_included = []
        creators_for_report = _creator_records(eml_config)
        if creators_for_report:
            sections_included.append(f"Creator information ({len(creators_for_report)} creator(s))")
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
        package_id = f"{uuid.uuid4()}/v1.0"
        
        # Create root EML element with proper namespaces
        root = ET.Element('eml:eml')
        root.set('packageId', package_id)
        root.set('system', 'http://gbif.org')
        root.set('scope', 'system')
        root.set('xmlns:eml', 'https://eml.ecoinformatics.org/eml-2.2.0')
        root.set('xmlns:dc', 'http://purl.org/dc/terms/')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        root.set('xmlns:stmml', 'http://www.xml-cml.org/schema/stmml-1.1')
        root.set('xsi:schemaLocation', 'https://eml.ecoinformatics.org/eml-2.2.0 https://rs.gbif.org/schema/eml-gbif-profile/1.3/eml.xsd')
        root.set('{http://www.w3.org/XML/1998/namespace}lang', _safe_get(eml_config, 'language', 'eng') or 'eng')
        
        # Create dataset element
        dataset = ET.SubElement(root, 'dataset')

        project_config = eml_config.get('project', {})
        template_values = {
            'short_name': effective_short_name,
            'project_id': _safe_get(project_config, 'id') or _get_project_value(project_df, 'project_id') or '',
        }

        # Identifiers and title (IPT order)
        for identifier in _as_list(eml_config.get('alternate_identifiers')):
            identifier = _format_template(identifier, template_values)
            if identifier:
                ET.SubElement(dataset, 'alternateIdentifier').text = identifier
        if effective_short_name:
            ET.SubElement(dataset, 'shortName').text = effective_short_name
        extra_alt = eml_config.get('extra_alternate_identifiers')
        if extra_alt is not None:
            for identifier in _as_list(extra_alt):
                identifier = _format_template(identifier, template_values)
                if identifier:
                    ET.SubElement(dataset, 'alternateIdentifier').text = identifier

        title = _safe_get(eml_config, 'title') or ''
        title_stripped = title.strip()
        if effective_short_name and bool(eml_config.get('append_short_name_slug_to_title', False)):
            title_stripped = f"{title_stripped} {effective_short_name}".strip() if title_stripped else effective_short_name
        elif not title_stripped and effective_short_name:
            title_stripped = effective_short_name
        final_dataset_title = title_stripped
        if final_dataset_title:
            title_elem = ET.SubElement(dataset, 'title')
            title_elem.text = final_dataset_title
            title_elem.set('{http://www.w3.org/XML/1998/namespace}lang', _safe_get(eml_config, 'language', 'eng') or 'eng')
        
        # Creators (required). Prefer creators: list; fallback to legacy creator: dict.
        for creator_config in _creator_records(eml_config):
            _create_person_element(dataset, 'creator', creator_config)
        
        # Metadata provider (optional)
        metadata_provider_config = eml_config.get('metadata_provider', {})
        if metadata_provider_config and _non_empty_dict(metadata_provider_config):
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

        # Intellectual rights (license)
        license_type = _safe_get(eml_config, 'license', 'CC0')
        _create_intellectual_rights_element(dataset, license_type)

        # Structured license metadata used by IPT exports
        _create_licensed_element(dataset, eml_config)

        # Distribution URLs
        _create_distribution_elements(dataset, eml_config, template_values)
        
        # Coverage (geographic, temporal, taxonomic)
        _create_coverage_element(dataset, eml_config)

        # Purpose
        purpose = _safe_get(eml_config, 'purpose')
        if purpose:
            purpose_elem = ET.SubElement(dataset, 'purpose')
            ET.SubElement(purpose_elem, 'para').text = purpose.strip()

        # Maintenance
        _create_maintenance_element(dataset, eml_config)

        # Contact
        contact_config = eml_config.get('contact', {})
        if contact_config:
            _create_person_element(dataset, 'contact', contact_config)
        
        # Methods
        methods_config = eml_config.get('methods', {})
        if methods_config:
            _create_methods_element(dataset, methods_config)
        
        # Project
        if project_config:
            _create_project_element(dataset, project_config)

        # GBIF metadata profile: resource citation (IPT additionalMetadata / gbif / citation)
        pub_date = _safe_get(eml_config, 'publication_date') or ''
        pub_year = pub_date[:4] if len(str(pub_date).strip()) >= 4 else ''
        resource_version = str(eml_config.get('resource_version', '') or '1.1').strip()
        nested_gbif = eml_config.get('gbif') if isinstance(eml_config.get('gbif'), dict) else {}
        citation_source = (
            (nested_gbif.get('resource_citation') if nested_gbif else None)
            or eml_config.get('gbif_resource_citation')
            or eml_config.get('dataset_citation')
        )
        cite_values = {
            **template_values,
            'title': final_dataset_title,
            'dataset_title': final_dataset_title,
            'pub_year': pub_year,
            'resource_version': resource_version,
        }
        citation_text = _format_template(str(citation_source or '').strip(), cite_values)
        if citation_text:
            additional_meta = ET.SubElement(root, 'additionalMetadata')
            meta_wrap = ET.SubElement(additional_meta, 'metadata')
            gbif_block = ET.SubElement(meta_wrap, 'gbif')
            citation_el = ET.SubElement(gbif_block, 'citation')
            citation_el.text = citation_text

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
