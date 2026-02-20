import csv
import os
from xml.dom import minidom
from xml.etree.ElementTree import Element, SubElement, tostring


_DWC_TEXT_NS = "http://rs.tdwg.org/dwc/text/"
_DWC_TEXT_SCHEMA = "http://rs.tdwg.org/dwc/text/tdwg_dwc_text.xsd"
_DWC_TERMS_NS = "http://rs.tdwg.org/dwc/terms/"
_GBIF_TERMS_NS = "http://rs.gbif.org/terms/1.0/"


def _read_csv_header(csv_path):
    with open(csv_path, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if row:
                return [c.strip() for c in row]
    return []


def _term_uri(term_name, namespace_uri):
    return f"{namespace_uri}{term_name}"


def _write_pretty_xml(elem, out_path):
    xml_str = minidom.parseString(tostring(elem, encoding="utf-8")).toprettyxml(
        indent="    ", encoding="utf-8"
    )
    with open(out_path, "wb") as f:
        f.write(xml_str)


def create_meta_xml(
    output_dir,
    core_filename,
    extension_filenames=None,
    metadata_filename=None,
    reporter=None,
):
    """
    Create a Darwin Core Archive meta.xml for the files in output_dir.
    """
    extension_filenames = extension_filenames or []

    core_path = os.path.join(output_dir, core_filename)
    if not os.path.exists(core_path):
        raise FileNotFoundError(f"Core file not found: {core_path}")

    core_header = _read_csv_header(core_path)
    if not core_header:
        raise ValueError(f"Core file has no header row: {core_path}")

    if "occurrenceID" not in core_header:
        raise ValueError(f"Core file missing required 'occurrenceID' column: {core_path}")

    occ_id_idx = core_header.index("occurrenceID")

    archive_attrib = {
        "xmlns": _DWC_TEXT_NS,
        "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
        "xsi:schemaLocation": f"{_DWC_TEXT_NS} {_DWC_TEXT_SCHEMA}",
    }
    if metadata_filename:
        metadata_path = os.path.join(output_dir, metadata_filename)
        if os.path.exists(metadata_path):
            archive_attrib["metadata"] = metadata_filename

    archive = Element("archive", archive_attrib)

    core = SubElement(
        archive,
        "core",
        {
            "encoding": "UTF-8",
            "fieldsTerminatedBy": ",",
            "linesTerminatedBy": "\\n",
            "fieldsEnclosedBy": '"',
            "ignoreHeaderLines": "1",
            "rowType": _term_uri("Occurrence", _DWC_TERMS_NS),
        },
    )
    files = SubElement(core, "files")
    SubElement(files, "location").text = core_filename
    SubElement(core, "id", {"index": str(occ_id_idx)})

    for idx, col in enumerate(core_header):
        SubElement(core, "field", {"index": str(idx), "term": _term_uri(col, _DWC_TERMS_NS)})

    for ext_filename in extension_filenames:
        ext_path = os.path.join(output_dir, ext_filename)
        if not os.path.exists(ext_path):
            continue

        ext_header = _read_csv_header(ext_path)
        if not ext_header:
            continue

        if "occurrenceID" not in ext_header:
            continue

        coreid_idx = ext_header.index("occurrenceID")

        if ext_filename.lower().startswith("emof"):
            row_type = _term_uri("MeasurementOrFact", _DWC_TERMS_NS)
            default_ns = _DWC_TERMS_NS
        elif ext_filename.lower().startswith("dna_derived_extension"):
            row_type = _term_uri("DNADerivedData", _GBIF_TERMS_NS)
            default_ns = _GBIF_TERMS_NS
        else:
            row_type = _term_uri("MeasurementOrFact", _DWC_TERMS_NS)
            default_ns = _DWC_TERMS_NS

        ext = SubElement(
            archive,
            "extension",
            {
                "encoding": "UTF-8",
                "fieldsTerminatedBy": ",",
                "linesTerminatedBy": "\\n",
                "fieldsEnclosedBy": '"',
                "ignoreHeaderLines": "1",
                "rowType": row_type,
            },
        )
        ext_files = SubElement(ext, "files")
        SubElement(ext_files, "location").text = ext_filename
        SubElement(ext, "coreid", {"index": str(coreid_idx)})

        dna_core_terms = {"occurrenceID", "eventID", "parentEventID", "materialSampleID"}
        for idx, col in enumerate(ext_header):
            ns = default_ns
            if default_ns == _GBIF_TERMS_NS and col in dna_core_terms:
                ns = _DWC_TERMS_NS
            SubElement(ext, "field", {"index": str(idx), "term": _term_uri(col, ns)})

    out_path = os.path.join(output_dir, "meta.xml")
    _write_pretty_xml(archive, out_path)

    if reporter:
        reporter.add_text(f"Created Darwin Core Archive metafile: {out_path}")

    return out_path

