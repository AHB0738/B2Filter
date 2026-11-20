"""High-level input, filtering, scoring and report-generation API."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from .core import (
    REQUIRED_COLUMNS, CATEGORY_ORDER, get_genotype_category, genotype_score,
    gtex_score, score_from_frequency, consequence_score, biallelic_score,
    compute_b2patho_thresholds
)
from .annotations import compute_pubscore_strg
from .export import write_excel

def _read_input(input_file):
    path = Path(input_file)
    # VCF-like tabular exports and Strelka/GATK TSV files are tab-delimited.
    sep = "\t" if path.suffix.lower() in {".vcf", ".tsv"} or ".strelka" in path.name.lower() or ".gatk" in path.name.lower() else ","
    return pd.read_csv(path, sep=sep, comment="#" if path.suffix.lower() == ".vcf" else None)

def _validate_columns(df):
    missing = sorted(REQUIRED_COLUMNS - set(df.columns))
    if missing:
        raise ValueError("Input is missing required columns: " + ", ".join(missing))

def _filter_categories(df):
    return {
        "HV-H": df[(df["REF"] != "N") & (df["IMPACT"] == "HIGH") & (df["GENOTYPE"] == "HV") & (df["CANONICAL"] == "YES")].copy(),
        "HV-MH": df[(df["REF"] != "N") & (df["IMPACT"] == "MODHIGH") & (df["GENOTYPE"] == "HV") & (df["CANONICAL"] == "YES")].copy(),
        "HET-bH": df[(df["REF"] != "N") & (df["BIALLELIC"] == "HIGH") & (df["GENOTYPE"] == "HET") & (df["CANONICAL"] == "YES")].copy(),
        "HET-bMH": df[(df["REF"] != "N") & (df["BIALLELIC"] == "MODHIGH") & (df["GENOTYPE"] == "HET") & (df["CANONICAL"] == "YES")].copy(),
        "Kn-bML": df[(df["REF"] != "N") & df["KNOWN_CANDIDATE_GENE"].notna() & df["GENOTYPE"].isin(["HV", "HET"]) & df["BIALLELIC"].isin(["LOW", "MODERATE"]) & (df["CANONICAL"] == "YES")].copy(),
        "Kn-nb": df[(df["REF"] != "N") & df["KNOWN_CANDIDATE_GENE"].notna() & (df["BIALLELIC"] == "NO") & (df["GENOTYPE"] == "HET") & df["IMPACT"].isin(["HIGH", "MODHIGH"]) & (df["CANONICAL"] == "YES")].copy(),
        "CNV": df[df["REF"] == "N"].copy(),
    }

def run_dataframe(data, output_dir=None, output_name="B2Filter_results.xlsx",
                  pubmed=True, request_delay=0.2, top_n=5, write_report=True):
    """Run B²Filter on a pandas DataFrame.

    Returns a dictionary containing filtered category DataFrames, the ranked
    B²Pred DataFrame, and (if requested) the generated Excel report path.
    """
    df = data.copy()
    _validate_columns(df)
    categories = _filter_categories(df)
    selected = pd.concat(categories.values(), ignore_index=True)
    symbols = selected["SYMBOL"].dropna().astype(str).unique().tolist()
    pub_scores = compute_pubscore_strg(symbols, request_delay=request_delay) if pubmed else {s: 0 for s in symbols}
    df["PubScoreStrg"] = df["SYMBOL"].map(pub_scores).fillna(0).astype(int)

    ranked = df[df["BIALLELIC"].fillna("").astype(str).str.upper() != "NO"].copy()
    if not ranked.empty:
        ranked["Category"] = ranked.apply(get_genotype_category, axis=1)
        ranked = ranked[
            (ranked["Category"] != "CNV") |
            ((ranked["Category"] == "CNV") & (ranked["PubScoreStrg"] >= 4))
        ].copy()
        if not ranked.empty:
            ranked["GtexScore"] = ranked["GTEX_testis_RATIO"].apply(gtex_score)
            ranked["GenotypeScore"] = ranked.apply(genotype_score, axis=1)
            ranked["GnomadScore"] = ranked["gnomADe_AF"].apply(score_from_frequency)
            ranked["ConsequenceScore"] = ranked.apply(consequence_score, axis=1)
            ranked["BiallelicScore"] = ranked["BIALLELIC"].apply(biallelic_score)
            cnv_mask = ranked["Category"] == "CNV"
            def _cnv_score(row):
                if row["Category"] != "CNV":
                    return 0
                genotype = str(row.get("GENOTYPE", ""))
                impact = str(row.get("IMPACT", ""))
                biallelic = str(row.get("BIALLELIC", ""))
                pub = row.get("PubScoreStrg", 0)
                gtex = gtex_score(row.get("GTEX_testis_RATIO"))
                if genotype == "HV" and impact == "HIGH" and biallelic == "HIGH":
                    return 35 + pub + gtex
                if genotype == "HET" and impact == "HIGH" and biallelic == "HIGH":
                    return 32 + pub + gtex
                if genotype == "HET" and impact == "HIGH" and biallelic == "MODHIGH":
                    return 30 + pub + gtex
                return pub + gtex
            ranked["CNVScore"] = ranked.apply(_cnv_score, axis=1)
            ranked["TotalScore"] = (
                ranked["GtexScore"] + ranked["GenotypeScore"] +
                ranked["BiallelicScore"] + ranked["GnomadScore"] +
                ranked["ConsequenceScore"] + ranked["PubScoreStrg"]
            )
            ranked.loc[cnv_mask, "TotalScore"] = ranked.loc[cnv_mask, "CNVScore"]
            ranked["Gnomad_Pub_Sum"] = ranked["GnomadScore"] + ranked["PubScoreStrg"]
            ranked = ranked.sort_values(["TotalScore", "Gnomad_Pub_Sum"], ascending=[False, False])
            ranked = ranked.drop(columns=["Gnomad_Pub_Sum"])
            ranked["Rank"] = range(1, len(ranked) + 1)
            ranked["CNV_Top"] = cnv_mask.loc[ranked.index] & (ranked["Rank"] <= 10)
            patho = compute_b2patho_thresholds(ranked["TotalScore"].max())
            ranked["B²Patho"] = ranked.apply(lambda row: patho(row["TotalScore"], row["Category"]), axis=1)
            ranked = ranked.head(top_n).copy()
            # Preserve the original HGVSc/HGVSp display convention without mutating input.
            for col in ("HGVSc", "HGVSp"):
                if col in ranked.columns:
                    ranked[col] = ranked[col].astype("string").str.split(":", n=1).str[-1]

    report = None
    if write_report:
        if output_dir is None:
            raise ValueError("output_dir is required when write_report=True")
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        report = write_excel(output_path / output_name, categories, ranked, df)
    return {"filtered": categories, "ranked": ranked, "report": report}

def run(input_file, output_dir=".", output_name=None, pubmed=True,
        request_delay=0.2, top_n=5):
    """Run B²Filter from a CSV/TSV/VCF-like file and write an Excel report."""
    input_path = Path(input_file)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")
    data = _read_input(input_path)
    name = output_name or f"{input_path.stem}_B2Filter.xlsx"
    return run_dataframe(
        data, output_dir=output_dir, output_name=name, pubmed=pubmed,
        request_delay=request_delay, top_n=top_n, write_report=True
    )
