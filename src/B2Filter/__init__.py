"""B²Filter: programmatic variant filtering and prioritization."""

from .pipeline import run, run_dataframe
from .core import (
    get_genotype_category,
    genotype_score,
    consequence_score,
    gtex_score,
    biallelic_score,
    score_from_frequency,
    score_sift,
    score_polyphen_am_spliceai,
    compute_b2patho_thresholds,
)

__version__ = "1.0.0"
__all__ = [
    "run", "run_dataframe", "get_genotype_category", "genotype_score",
    "consequence_score", "gtex_score", "biallelic_score",
    "score_from_frequency", "score_sift",
    "score_polyphen_am_spliceai", "compute_b2patho_thresholds",
]
