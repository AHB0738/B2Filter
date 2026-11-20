# B²Filter

**B²Filter** is a Python package for research-oriented filtering and prioritization of genetic variants, with a B²Pred ranking module. It is designed for programmatic use and Linux command-line workflows; it has no graphical interface and does not import Tkinter.

> **Research-use notice:** B²Filter is a research software tool, not a medical device or a stand-alone clinical diagnostic system. Candidate prioritization requires expert review, validation, and interpretation in the context of phenotype, inheritance, sequencing quality, and current evidence.

## Installation

Python 3.9 or newer is required.

```bash
python -m pip install b2filter
```

The package name on PyPI is lowercase (`b2filter`), while the Python import is case-sensitive:

```python
import B2Filter
```

## Python API

```python
import B2Filter

result = B2Filter.run(
    input_file="sample_variants.csv",
    output_dir="results",
    pubmed=True,
    top_n=5,
)

print(result["report"])
print(result["ranked"].head())
```

You can also pass a pandas DataFrame:

```python
import pandas as pd
import B2Filter

variants = pd.read_csv("sample_variants.csv")
result = B2Filter.run_dataframe(
    variants,
    output_dir="results",
    pubmed=False,  # no external HGNC/PubMed requests
)
```

## Linux command line

```bash
b2filter sample_variants.csv --output-dir results
```

Disable external PubMed/HGNC requests:

```bash
b2filter sample_variants.csv --output-dir results --no-pubmed
```

## Input

The input is expected to be a tabular annotation file containing the columns used by the original B²Filter workflow. The required fields are validated before processing. Common CSV, TSV and tab-delimited VCF-like exports are supported; raw VCF files with INFO-field annotations are not automatically converted into the expected annotation table.

Required columns include `REF`, `ALT`, `POSITION`, `SYMBOL`, `GENOTYPE`, `IMPACT`, `BIALLELIC`, `CANONICAL`, `Consequence`, `gnomADe_AF`, `GTEX_testis_RATIO`, `SIFT`, `PolyPhen`, `am_pathogenicity`, `SpliceAI_DS`, and `KNOWN_CANDIDATE_GENE`.

## Output

The Excel report contains:
- `B²Pred`: ranked candidate variants (Top 5 by default)
- `HV-H`, `HV-MH`, `HET-bH`, `HET-bMH`, `Kn-bML`, `Kn-nb`, `CNV`: category-specific tables
- `Conclusion`: summary sheet

Where the necessary columns are available, database hyperlinks are added for gnomAD, GTEx, ClinVar, OMIM, PubMed and UCSC.

## External services

When `pubmed=True`, B²Filter queries HGNC and NCBI PubMed over the internet to calculate `PubScoreStrg`. Network failures are handled by assigning a score of zero for the affected gene; this means results may differ when external services are unavailable or their databases change. For reproducible analyses, archive input files, package version, parameters, date, and external annotation results. You can set `pubmed=False` to disable these requests.

## Compatibility

The package metadata declares Python 3.9+. Compatibility should be verified in the target environment using the included tests. It is not technically possible to guarantee compatibility with every Python/Linux distribution or every version. The package avoids Tkinter and platform-specific GUI dependencies.

## Development and tests

```bash
git clone https://github.com/YOUR_GITHUB_USERNAME/B2Filter.git
cd B2Filter
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
pytest
```

## Citation

If you use B²Filter in a scientific publication, cite the versioned release and its archived DOI (for example, Zenodo) once created. Do not cite an untagged development branch as the stable release.

## License

MIT. Review authorship, institutional ownership, and permissions before publishing the repository under this license.
