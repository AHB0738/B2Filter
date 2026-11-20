import pandas as pd
import pytest

from B2Filter.core import get_genotype_category


@pytest.mark.parametrize(
    "row, expected",
    [
        ({"REF": "A", "IMPACT": "HIGH", "GENOTYPE": "HV",
          "CANONICAL": "YES", "BIALLELIC": "NO",
          "KNOWN_CANDIDATE_GENE": None}, "HV-H"),

        ({"REF": "A", "IMPACT": "MODHIGH", "GENOTYPE": "HV",
          "CANONICAL": "YES", "BIALLELIC": "NO",
          "KNOWN_CANDIDATE_GENE": None}, "HV-MH"),

        ({"REF": "A", "IMPACT": "HIGH", "GENOTYPE": "HET",
          "CANONICAL": "YES", "BIALLELIC": "HIGH",
          "KNOWN_CANDIDATE_GENE": None}, "HET-bH"),

        ({"REF": "A", "IMPACT": "MODHIGH", "GENOTYPE": "HET",
          "CANONICAL": "YES", "BIALLELIC": "MODHIGH",
          "KNOWN_CANDIDATE_GENE": None}, "HET-bMH"),

        ({"REF": "A", "IMPACT": "LOW", "GENOTYPE": "HET",
          "CANONICAL": "YES", "BIALLELIC": "LOW",
          "KNOWN_CANDIDATE_GENE": "NOA:3"}, "Kn-bML"),

        ({"REF": "A", "IMPACT": "HIGH", "GENOTYPE": "HET",
          "CANONICAL": "YES", "BIALLELIC": "NO",
          "KNOWN_CANDIDATE_GENE": "NOA:3"}, "Kn-nb"),

        ({"REF": "N", "IMPACT": "HIGH", "GENOTYPE": "HV",
          "CANONICAL": "YES", "BIALLELIC": "HIGH",
          "KNOWN_CANDIDATE_GENE": None}, "CNV"),
    ],
)
def test_all_seven_categories(row, expected):
    assert get_genotype_category(pd.Series(row)) == expected


def test_unclassified_variant():
    row = pd.Series({
        "REF": "A",
        "IMPACT": "LOW",
        "GENOTYPE": "HET",
        "CANONICAL": "NO",
        "BIALLELIC": "NO",
        "KNOWN_CANDIDATE_GENE": None,
    })
    assert get_genotype_category(row) is None
