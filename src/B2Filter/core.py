"""Core B²Filter/B²Pred rules. No GUI or file-system interaction in this module."""
from __future__ import annotations
import re
import pandas as pd

REQUIRED_COLUMNS = {
    "REF", "ALT", "POSITION", "SYMBOL", "GENOTYPE", "IMPACT",
    "BIALLELIC", "CANONICAL", "Consequence", "gnomADe_AF",
    "GTEX_testis_RATIO", "SIFT", "PolyPhen",
    "KNOWN_CANDIDATE_GENE",
}

CATEGORY_ORDER = ["HV-H", "HV-MH", "HET-bH", "HET-bMH", "Kn-bML", "Kn-nb", "CNV"]

def extract_number_from_parentheses(value):
    if pd.isna(value):
        return None
    match = re.search(r"\(([0-9.,]+)\)", str(value))
    if match:
        try:
            return float(match.group(1).replace(",", "."))
        except ValueError:
            return None
    try:
        return float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return None

def score_from_frequency(af, thresholds=(0.0001, 0.001, 0.01, 0.1), scores=(5, 4, 3, 2, 1)):
    if pd.isna(af):
        return 5
    try:
        af = float(af)
    except (TypeError, ValueError):
        return 5
    for threshold, score in zip(thresholds, scores[:-1]):
        if af < threshold:
            return score
    return scores[-1]

def score_sift(value):
    value = extract_number_from_parentheses(value)
    if value is None:
        return 0
    if value < .2: return 5
    if value < .4: return 4
    if value < .6: return 3
    if value < .8: return 2
    return 1

def score_polyphen_am_spliceai(value, inverse=True):
    value = extract_number_from_parentheses(value) if isinstance(value, str) else value
    if value is None or pd.isna(value):
        return 0
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 0
    if value < .2: return 1
    if value < .4: return 2
    if value < .6: return 3
    if value < .8: return 4
    return 5

def get_genotype_category(row):
    ref = row.get("REF")
    impact = row.get("IMPACT")
    genotype = row.get("GENOTYPE")
    canonical = row.get("CANONICAL")
    biallelic = row.get("BIALLELIC")
    known = pd.notna(row.get("KNOWN_CANDIDATE_GENE"))
    if ref != "N" and impact == "HIGH" and genotype == "HV" and canonical == "YES":
        return "HV-H"
    if ref != "N" and impact == "MODHIGH" and genotype == "HV" and canonical == "YES":
        return "HV-MH"
    if ref != "N" and biallelic == "HIGH" and genotype == "HET" and canonical == "YES":
        return "HET-bH"
    if ref != "N" and biallelic == "MODHIGH" and genotype == "HET" and canonical == "YES":
        return "HET-bMH"
    if ref != "N" and known and genotype in ("HV", "HET") and biallelic in ("LOW", "MODERATE") and canonical == "YES":
        return "Kn-bML"
    if ref != "N" and known and biallelic == "NO" and genotype == "HET" and impact in ("HIGH", "MODHIGH") and canonical == "YES":
        return "Kn-nb"
    if ref == "N":
        return "CNV"
    return None

def genotype_score(row):
    category = row.get("Category")
    if category == "Kn-nb":
        return 0
    genotype = str(row.get("GENOTYPE", "")).upper()
    impact = str(row.get("IMPACT", "")).upper()
    return {"HV": 1, "HET": 0}.get(genotype, 0) + {
        "HIGH": 4, "MODHIGH": 3, "MODERATE": 2, "LOW": 1, "NO": 0
    }.get(impact, 0)

def gtex_score(value):
    if pd.isna(value):
        return 1
    try:
        value = float(value)
    except (TypeError, ValueError):
        return 1
    if value < 1: return 1
    if value < 2: return 2
    if value < 5: return 3
    if value < 10: return 4
    return 5

def consequence_score(row):
    category = get_genotype_category(row)
    consequence = str(row.get("Consequence", "")).strip()
    if category in ("HV-H", "HV-MH", "HET-bH"):
        return 12
    lof = {
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
        "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
        "transcript_amplification", "feature_elongation", "feature_truncation"
    }
    if category == "Kn-bML":
        terms = set(consequence.split("&"))
        if terms & lof:
            return 12
        if any(term.startswith("inframe_") for term in terms):
            return 8
    if "missense_variant" in consequence:
        return (
            score_sift(row.get("SIFT"))
            + score_polyphen_am_spliceai(row.get("PolyPhen"))
            + score_polyphen_am_spliceai(row.get("am_pathogenicity"))
            + score_polyphen_am_spliceai(row.get("SpliceAI_DS"))
        )
    if consequence.startswith("splice_") or consequence == "synonymous_variant":
        return score_polyphen_am_spliceai(row.get("SpliceAI_DS"))
    return 4

def biallelic_score(value):
    if pd.isna(value):
        return 0
    return {"HIGH": 5, "MODHIGH": 4, "MODERATE": 3, "LOW": 2, "NO": 0}.get(str(value).upper(), 0)

def compute_b2patho_thresholds(max_score=None):
    """Return the original B²Patho mapping (fixed cut-offs: 9, 18, 27, 36)."""
    def patho_score(score, category=None):
        if score <= 9: return 1
        if score <= 18: return 2
        if score <= 27: return 3
        if score <= 36: return 4
        return 5
    return patho_score
