
"""High-level input, filtering, scoring and report-generation API."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .core import (
    REQUIRED_COLUMNS,
    get_genotype_category,
    genotype_score,
    gtex_score,
    score_from_frequency,
    consequence_score,
    biallelic_score,
    compute_b2patho_thresholds,
)
from .annotations import compute_pubscore_strg
from .export import write_excel


def _read_input(input_file):
    """Read CSV, TSV, Strelka or GATK tabular input."""
    path = Path(input_file)

    sep = (
        "\t"
        if (
            path.suffix.lower() in {".vcf", ".tsv"}
            or ".strelka" in path.name.lower()
            or ".gatk" in path.name.lower()
        )
        else ","
    )

    return pd.read_csv(
        path,
        sep=sep,
        comment="#" if path.suffix.lower() == ".vcf" else None,
        low_memory=False,
    )


def _normalize_annotations(df):
    """
    Normalize gene symbols and optional annotations.

    Supports:
    - Cleaning SYMBOL values from Strelka/VCF exports.
    - Native SpliceAI_DS.
    - Ensembl VEP SpliceAI delta scores.
    - Missing AlphaMissense annotations.

    Missing annotations remain NaN.
    Missing gene symbols remain pd.NA.
    """
    df = df.copy()

    # Step 1: Clean gene symbols.
    if "SYMBOL" in df.columns:
        df["SYMBOL"] = (
            df["SYMBOL"]
            .astype("string")
            .str.strip()
            .str.lstrip("'")
            .str.strip()
            .replace("", pd.NA)
        )

    # Step 2: Normalize SpliceAI annotations.
    spliceai_columns = [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
    ]

    # Preserve native SpliceAI_DS when already available.
    if "SpliceAI_DS" not in df.columns:

        available = [
            col for col in spliceai_columns
            if col in df.columns
        ]

        if available:
            scores = df[available].apply(
                pd.to_numeric,
                errors="coerce",
            )

            df["SpliceAI_DS"] = scores.max(
                axis=1,
                skipna=True,
            )
        else:
            df["SpliceAI_DS"] = float("nan")

    else:
        df["SpliceAI_DS"] = pd.to_numeric(
            df["SpliceAI_DS"],
            errors="coerce",
        )

    # Step 3: AlphaMissense is optional.
    if "am_pathogenicity" not in df.columns:
        df["am_pathogenicity"] = float("nan")
    else:
        df["am_pathogenicity"] = pd.to_numeric(
            df["am_pathogenicity"],
            errors="coerce",
        )

    return df


def _validate_columns(df):
    """Validate the required input columns."""
    missing = sorted(
        REQUIRED_COLUMNS - set(df.columns)
    )

    if missing:
        raise ValueError(
            "Input is missing required columns: "
            + ", ".join(missing)
        )


def _filter_categories(df):
    """Apply the seven original B²Filter classification rules."""
    return {
        "HV-H": df[
            (df["REF"] != "N")
            & (df["IMPACT"] == "HIGH")
            & (df["GENOTYPE"] == "HV")
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "HV-MH": df[
            (df["REF"] != "N")
            & (df["IMPACT"] == "MODHIGH")
            & (df["GENOTYPE"] == "HV")
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "HET-bH": df[
            (df["REF"] != "N")
            & (df["BIALLELIC"] == "HIGH")
            & (df["GENOTYPE"] == "HET")
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "HET-bMH": df[
            (df["REF"] != "N")
            & (df["BIALLELIC"] == "MODHIGH")
            & (df["GENOTYPE"] == "HET")
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "Kn-bML": df[
            (df["REF"] != "N")
            & df["KNOWN_CANDIDATE_GENE"].notna()
            & df["GENOTYPE"].isin(["HV", "HET"])
            & df["BIALLELIC"].isin(["LOW", "MODERATE"])
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "Kn-nb": df[
            (df["REF"] != "N")
            & df["KNOWN_CANDIDATE_GENE"].notna()
            & (df["BIALLELIC"] == "NO")
            & (df["GENOTYPE"] == "HET")
            & df["IMPACT"].isin(["HIGH", "MODHIGH"])
            & (df["CANONICAL"] == "YES")
        ].copy(),

        "CNV": df[
            df["REF"] == "N"
        ].copy(),
    }


def run_dataframe(
    data,
    output_dir=None,
    output_name="B2Filter_results.xlsx",
    pubmed=True,
    request_delay=0.2,
    top_n=5,
    write_report=True,
):
    """
    Run B²Filter on a pandas DataFrame.

    B²Pred ranks exclusively the variants selected by the seven
    B²Filter categories.

    Returns:
        - filtered: category DataFrames
        - ranked: B²Pred ranked DataFrame
        - report: generated Excel report path
    """
    df = data.copy()

    # Step 0: Normalize symbols and optional annotations.
    df = _normalize_annotations(df)

    # Validate mandatory columns.
    _validate_columns(df)

    # Step 1: Apply the seven original categories.
    categories = _filter_categories(df)

    # Step 2: Combine selected variants.
    selected = pd.concat(
        categories.values(),
        ignore_index=True,
    )

    # Step 3: Retrieve PubMed scores for selected genes.
    symbols = (
        selected["SYMBOL"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    pub_scores = (
        compute_pubscore_strg(
            symbols,
            request_delay=request_delay,
        )
        if pubmed
        else {symbol: 0 for symbol in symbols}
    )

    # Add PubMed scores to original input.
    df["PubScoreStrg"] = (
        df["SYMBOL"]
        .map(pub_scores)
        .fillna(0)
        .astype(int)
    )

    # Add PubMed scores to selected variants.
    selected["PubScoreStrg"] = (
        selected["SYMBOL"]
        .map(pub_scores)
        .fillna(0)
        .astype(int)
    )

    # Step 4: B²Pred processes selected variants only.
    ranked = selected.copy()

    if not ranked.empty:

        ranked["Category"] = ranked.apply(
            get_genotype_category,
            axis=1,
        )

        # Preserve original CNV PubMed eligibility rule.
        ranked = ranked[
            (ranked["Category"] != "CNV")
            | (
                (ranked["Category"] == "CNV")
                & (ranked["PubScoreStrg"] >= 4)
            )
        ].copy()

        if not ranked.empty:

            # Step 5: Calculate individual scores.
            ranked["GtexScore"] = (
                ranked["GTEX_testis_RATIO"].apply(
                    gtex_score
                )
            )

            ranked["GenotypeScore"] = ranked.apply(
                genotype_score,
                axis=1,
            )

            ranked["GnomadScore"] = (
                ranked["gnomADe_AF"].apply(
                    score_from_frequency
                )
            )

            ranked["ConsequenceScore"] = (
                ranked.apply(
                    consequence_score,
                    axis=1,
                )
            )

            ranked["BiallelicScore"] = (
                ranked["BIALLELIC"].apply(
                    biallelic_score
                )
            )

            # Step 6: Original CNV-specific score.
            cnv_mask = ranked["Category"] == "CNV"

            def _cnv_score(row):
                if row["Category"] != "CNV":
                    return 0

                genotype = str(row.get("GENOTYPE", ""))
                impact = str(row.get("IMPACT", ""))
                biallelic = str(row.get("BIALLELIC", ""))

                pub = row.get("PubScoreStrg", 0)
                gtex = gtex_score(
                    row.get("GTEX_testis_RATIO")
                )

                if (
                    genotype == "HV"
                    and impact == "HIGH"
                    and biallelic == "HIGH"
                ):
                    return 35 + pub + gtex

                if (
                    genotype == "HET"
                    and impact == "HIGH"
                    and biallelic == "HIGH"
                ):
                    return 32 + pub + gtex

                if (
                    genotype == "HET"
                    and impact == "HIGH"
                    and biallelic == "MODHIGH"
                ):
                    return 30 + pub + gtex

                return pub + gtex

            ranked["CNVScore"] = ranked.apply(
                _cnv_score,
                axis=1,
            )

            # Step 7: Calculate total scores.
            ranked["TotalScore"] = (
                ranked["GtexScore"]
                + ranked["GenotypeScore"]
                + ranked["BiallelicScore"]
                + ranked["GnomadScore"]
                + ranked["ConsequenceScore"]
                + ranked["PubScoreStrg"]
            )

            ranked.loc[cnv_mask, "TotalScore"] = (
                ranked.loc[cnv_mask, "CNVScore"]
            )

            # Step 8: Sort variants.
            ranked["Gnomad_Pub_Sum"] = (
                ranked["GnomadScore"]
                + ranked["PubScoreStrg"]
            )

            ranked = ranked.sort_values(
                ["TotalScore", "Gnomad_Pub_Sum"],
                ascending=[False, False],
            )

            ranked = ranked.drop(
                columns=["Gnomad_Pub_Sum"]
            )

            # Step 9: Assign ranks and CNV top indicators.
            ranked["Rank"] = range(
                1,
                len(ranked) + 1,
            )

            ranked["CNV_Top"] = (
                (ranked["Category"] == "CNV")
                & (ranked["Rank"] <= 10)
            )

            # Step 10: B²Patho classification.
            patho = compute_b2patho_thresholds(
                ranked["TotalScore"].max()
            )

            ranked["B²Patho"] = ranked.apply(
                lambda row: patho(
                    row["TotalScore"],
                    row["Category"],
                ),
                axis=1,
            )

            # Step 11: Keep requested top variants.
            ranked = ranked.head(top_n).copy()

            # Preserve original HGVSc/HGVSp convention.
            for col in ("HGVSc", "HGVSp"):
                if col in ranked.columns:
                    ranked[col] = (
                        ranked[col]
                        .astype("string")
                        .str.split(":", n=1)
                        .str[-1]
                    )

    # Step 12: Generate Excel report.
    report = None

    if write_report:

        if output_dir is None:
            raise ValueError(
                "output_dir is required when write_report=True"
            )

        output_path = Path(output_dir)

        output_path.mkdir(
            parents=True,
            exist_ok=True,
        )

        report = write_excel(
            output_path / output_name,
            categories,
            ranked,
            df,
        )

    return {
        "filtered": categories,
        "ranked": ranked,
        "report": report,
    }


def run(
    input_file,
    output_dir=".",
    output_name=None,
    pubmed=True,
    request_delay=0.2,
    top_n=5,
):
    """Run B²Filter from CSV/TSV/VCF-like input."""
    input_path = Path(input_file)

    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input file does not exist: {input_path}"
        )

    data = _read_input(input_path)

    name = output_name or (
        f"{input_path.stem}_B2Filter.xlsx"
    )

    return run_dataframe(
        data,
        output_dir=output_dir,
        output_name=name,
        pubmed=pubmed,
        request_delay=request_delay,
        top_n=top_n,
        write_report=True,
    )
