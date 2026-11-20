import pandas as pd
from B2Filter.core import (
    get_genotype_category, score_sift, gtex_score,
    score_from_frequency, biallelic_score, compute_b2patho_thresholds
)

def test_category_hv_high():
    row = pd.Series({
        "REF": "A", "IMPACT": "HIGH", "GENOTYPE": "HV",
        "CANONICAL": "YES", "BIALLELIC": "NO",
        "KNOWN_CANDIDATE_GENE": None
    })
    assert get_genotype_category(row) == "HV-H"

def test_score_functions():
    assert score_sift("deleterious(0.02)") == 5
    assert gtex_score(6) == 4
    assert score_from_frequency(0.00001) == 5
    assert biallelic_score("HIGH") == 5
    assert compute_b2patho_thresholds()(37) == 5
