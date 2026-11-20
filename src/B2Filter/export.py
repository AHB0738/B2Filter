"""Excel report writer for B²Filter results."""
from __future__ import annotations
from pathlib import Path
from copy import copy
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

HEADER_FILL = PatternFill("solid", fgColor="119C9C")
HEADER_FONT = Font(color="FFFFFF", bold=True)
CATEGORY_ORDER = ["HV-H", "HV-MH", "HET-bH", "HET-bMH", "Kn-bML", "Kn-nb", "CNV"]

def _style_sheet(ws):
    for col in ws.columns:
        letter = col[0].column_letter
        max_len = max((len(str(c.value)) for c in col if c.value is not None), default=0)
        ws.column_dimensions[letter].width = min(max(max_len + 2, 12), 55)
        for cell in col:
            cell.alignment = Alignment(horizontal="center", vertical="center")
    if ws.max_row:
        for cell in ws[1]:
            cell.fill = HEADER_FILL
            cell.font = HEADER_FONT
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

def add_resource_links(path):
    """Add convenient external database links when the required columns exist."""
    from openpyxl import load_workbook
    wb = load_workbook(path)
    blue = Font(color="0000FF", underline="single")
    for ws in wb.worksheets:
        if ws.title == "Conclusion" or ws.max_row < 2:
            continue
        headers = {c.value: c.column for c in ws[1]}
        symbol_col = headers.get("SYMBOL")
        position_col, ref_col, alt_col = headers.get("POSITION"), headers.get("REF"), headers.get("ALT")
        mane_col, hgvsc_col = headers.get("MANE_SELECT"), headers.get("HGVSc")
        def append_link_column(title, url_builder):
            col = ws.max_column + 1
            ws.cell(1, col, title).fill = HEADER_FILL
            ws.cell(1, col).font = HEADER_FONT
            for r in range(2, ws.max_row + 1):
                url = url_builder(r)
                if url:
                    cell = ws.cell(r, col, ws.cell(r, symbol_col).value if symbol_col else title)
                    cell.hyperlink = url
                    cell.font = blue
        if symbol_col:
            append_link_column("GTEx", lambda r: f"https://gtexportal.org/home/gene/{ws.cell(r, symbol_col).value}" if ws.cell(r, symbol_col).value else None)
            append_link_column("OMIM", lambda r: f"https://www.omim.org/search/?search={ws.cell(r, symbol_col).value}" if ws.cell(r, symbol_col).value else None)
            append_link_column("PubMed", lambda r: f"https://pubmed.ncbi.nlm.nih.gov/?term={ws.cell(r, symbol_col).value}" if ws.cell(r, symbol_col).value else None)
        if position_col and ref_col and alt_col:
            append_link_column("gnomAD", lambda r: _gnomad_url(ws, r, position_col, ref_col, alt_col))
            append_link_column("UCSC", lambda r: _ucsc_url(ws, r, position_col, alt_col))
        if symbol_col and mane_col and hgvsc_col:
            append_link_column("ClinVar", lambda r: _clinvar_url(ws, r, symbol_col, mane_col, hgvsc_col))
        _style_sheet(ws)
    wb.save(path)

def _gnomad_url(ws, r, pcol, refcol, altcol):
    pos, ref, alt = ws.cell(r,pcol).value, ws.cell(r,refcol).value, ws.cell(r,altcol).value
    if not pos or not ref or not alt: return None
    cleaned = str(pos).strip().replace("chr", "").replace(":", "-")
    return f"https://gnomad.broadinstitute.org/variant/{cleaned}-{ref}-{alt}"

def _ucsc_url(ws, r, pcol, altcol):
    pos, alt = ws.cell(r,pcol).value, ws.cell(r,altcol).value
    if not pos or not alt or ":" not in str(pos): return None
    chrom, start = str(pos).split(":", 1)
    if str(alt).startswith("<DEL:") and str(alt).endswith(">"):
        region = f"{chrom}:{start}-{str(alt)[5:-1]}"
    else:
        region = f"{chrom}:{start}"
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={region}"

def _clinvar_url(ws, r, symbol_col, mane_col, hgvsc_col):
    gene, mane, hgvsc = (ws.cell(r,c).value for c in (symbol_col,mane_col,hgvsc_col))
    if not gene or not mane or not hgvsc: return None
    from urllib.parse import quote
    term = f"{mane}({gene}):{hgvsc}"
    return "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + quote(term)

def write_excel(output_path, filtered, ranked, original):
    output_path = Path(output_path)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        if not ranked.empty:
            ranked.to_excel(writer, sheet_name="B²Pred", index=False)
        for category in CATEGORY_ORDER:
            df = filtered.get(category, pd.DataFrame())
            df.to_excel(writer, sheet_name=category, index=False)
        coverage_cols = [c for c in original.columns if c.lower().startswith("coverage")]
        summary = [{"Couverture": "Conclusion"}, {"Couverture": "B²Filter / B²Pred"}]
        if coverage_cols:
            for col in coverage_cols:
                values = original[col].dropna()
                summary[0][col] = values.iloc[0] if not values.empty else ""
        pd.DataFrame(summary).to_excel(writer, sheet_name="Conclusion", index=False)
        for ws in writer.book.worksheets:
            _style_sheet(ws)
    add_resource_links(output_path)
    return output_path
