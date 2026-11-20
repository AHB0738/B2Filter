import pandas as pd

from B2Filter.pipeline import run_dataframe



def test_kn_nb_variant_is_included_in_b2pred():
    data = pd.DataFrame([{
        "REF": "A",
        "ALT": "G",
        "POSITION": "chr1:100000",
        "SYMBOL": "GENE_TEST",
        "GENOTYPE": "HET",
        "IMPACT": "HIGH",
        "BIALLELIC": "NO",
        "CANONICAL": "YES",
        "Consequence": "missense_variant",
        "gnomADe_AF": 0.00001,
        "GTEX_testis_RATIO": 2.0,
        "SIFT": "deleterious(0.02)",
        "PolyPhen": "probably_damaging(0.99)",
        "am_pathogenicity": 0.9,
        "SpliceAI_DS": 0.8,
        "KNOWN_CANDIDATE_GENE": "NOA:3",
    }])

    result = run_dataframe(
        data,
        pubmed=False,
        write_report=False,
        top_n=5,
    )

    assert len(result["filtered"]["Kn-nb"]) == 1
    assert "GENE_TEST" in result["ranked"]["SYMBOL"].values
