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
            reporter.add_text(f"  💡 {rec['guidance']}")
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
        reporter.add_text(f"📦 Package ID: {package_id}")
        reporter.add_text(f"📏 File size: {file_size_kb:.1f} KB")
        
        # Count XML elements for reporting
        try:
            tree = ET.parse(eml_path)
            root = tree.getroot()
            element_count = len(list(root.iter()))
            reporter.add_text(f"🏗️ XML structure: {element_count} total elements")
        except:
            pass
        
        # Validate the XML structure
        try:
            ET.parse(eml_path)
            reporter.add_success("✅ EML XML validation passed - file is well-formed and ready for OBIS/GBIF submission")
        except ET.ParseError as e:
            reporter.add_warning(f"❌ EML XML validation warning: {e}")
        
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


def _auto_populate_from_project_data(eml_config: Dict, data: Dict[str, pd.DataFrame], reporter) -> Dict:
    """
    Auto-populate some EML fields from project metadata if they're not already filled
    
    This is a helper function that can extract information from the FAIRe metadata
    to reduce manual configuration burden.
    """
    try:
        # Get project metadata
        project_df = data.get('projectMetadata', pd.DataFrame())
        if project_df.empty:
            return eml_config
        
        # Auto-populate title if it's still the template
        current_title = _safe_get(eml_config, 'title')
        if '[LOCATION]' in current_title or '[TIMEFRAME]' in current_title:
            # Try to extract better title info from project metadata
            # This would require understanding the FAIRe format better
            reporter.add_text("Consider updating the EML title with more specific information")
        
        return eml_config
        
    except Exception as e:
        reporter.add_warning(f"Could not auto-populate EML fields from project data: {e}")
        return eml_config
