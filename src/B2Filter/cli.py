"""Command-line interface."""
import argparse
import json
from .pipeline import run

def main():
    parser = argparse.ArgumentParser(prog="b2filter", description="Run B²Filter/B²Pred variant prioritization.")
    parser.add_argument("input", help="Input CSV, TSV or VCF-like tabular file")
    parser.add_argument("-o", "--output-dir", default="results", help="Output directory (default: results)")
    parser.add_argument("--no-pubmed", action="store_true", help="Disable HGNC/PubMed annotation queries")
    parser.add_argument("--top", type=int, default=5, help="Number of ranked variants in B²Pred (default: 5)")
    parser.add_argument("--delay", type=float, default=0.2, help="Delay between PubMed queries in seconds")
    args = parser.parse_args()
    result = run(args.input, args.output_dir, pubmed=not args.no_pubmed, top_n=args.top, request_delay=args.delay)
    print(f"Report written to: {result['report']}")
    print(f"Variants in B²Pred: {len(result['ranked'])}")
    for category, frame in result["filtered"].items():
        print(f"{category}: {len(frame)}")
if __name__ == "__main__":
    main()
