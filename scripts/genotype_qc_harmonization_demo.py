#!/usr/bin/env python3
"""
Run a synthetic genotype QC and harmonization readiness workflow.

Created by: Mahdi Akbarzadeh
"""

from __future__ import annotations

import argparse
import csv
import html
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

VALID_ALLELES = {"A", "C", "G", "T"}
MISSING_GENOTYPES = {"--", "", "NA", "N/A", "."}
AMBIGUOUS_ALLELE_SETS = {frozenset({"A", "T"}), frozenset({"C", "G"})}
AUTOSOMES = {str(i) for i in range(1, 23)}


def read_tsv(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def is_missing_genotype(genotype: str) -> bool:
    return genotype.strip().upper() in MISSING_GENOTYPES


def is_valid_genotype(genotype: str) -> bool:
    genotype = genotype.strip().upper()
    if is_missing_genotype(genotype):
        return True
    if len(genotype) != 2:
        return False
    return all(allele in VALID_ALLELES for allele in genotype)


def genotype_allele_set(genotype: str) -> frozenset:
    genotype = genotype.strip().upper()
    if is_missing_genotype(genotype):
        return frozenset()
    return frozenset(genotype)


def is_ambiguous_variant(effect_allele: str, other_allele: str) -> bool:
    return frozenset({effect_allele.upper(), other_allele.upper()}) in AMBIGUOUS_ALLELE_SETS


def classify_harmonization(row: dict, reference_by_rsid: Dict[str, dict], genome_build: str) -> Tuple[str, str]:
    rsid = row["rsid"]
    genotype = row["genotype"].strip().upper()

    if rsid not in reference_by_rsid:
        return "WARN", "Variant is absent from the synthetic reference table."

    ref = reference_by_rsid[rsid]
    if ref.get("genome_build", "") != genome_build:
        return "FAIL", "Reference genome build does not match the declared input genome build."

    if is_missing_genotype(genotype):
        return "SKIP", "Missing genotype; harmonization skipped."

    if not is_valid_genotype(genotype):
        return "FAIL", "Invalid genotype encoding."

    ref_alleles = frozenset({ref["effect_allele"].upper(), ref["other_allele"].upper()})
    observed_alleles = genotype_allele_set(genotype)

    if is_ambiguous_variant(ref["effect_allele"], ref["other_allele"]):
        return "WARN", "Ambiguous A/T or C/G variant; real workflows require allele-frequency-aware handling."

    if observed_alleles.issubset(ref_alleles):
        return "PASS", "Observed genotype alleles are compatible with the reference alleles."

    return "FAIL", "Observed genotype alleles are not compatible with the reference alleles."


def build_variant_status(genotype_rows: List[dict], reference_rows: List[dict], genome_build: str) -> List[dict]:
    reference_by_rsid = {row["rsid"]: row for row in reference_rows}
    rsid_counts = Counter(row["rsid"] for row in genotype_rows)
    output_rows: List[dict] = []

    for row in genotype_rows:
        genotype = row["genotype"].strip().upper()
        rsid = row["rsid"]
        chromosome = row["chromosome"].strip()

        duplicate_status = "WARN" if rsid_counts[rsid] > 1 else "PASS"
        chromosome_status = "PASS" if chromosome in AUTOSOMES else "WARN"

        try:
            position = int(row["position"])
            position_status = "PASS" if position > 0 else "FAIL"
        except ValueError:
            position_status = "FAIL"

        genotype_status = "PASS" if is_valid_genotype(genotype) else "FAIL"
        missing_status = "SKIP" if is_missing_genotype(genotype) else "PASS"
        harmonization_status, harmonization_reason = classify_harmonization(
            row=row,
            reference_by_rsid=reference_by_rsid,
            genome_build=genome_build,
        )

        output_rows.append(
            {
                "rsid": rsid,
                "chromosome": chromosome,
                "position": row["position"],
                "genotype": genotype,
                "duplicate_status": duplicate_status,
                "chromosome_status": chromosome_status,
                "position_status": position_status,
                "genotype_status": genotype_status,
                "missing_status": missing_status,
                "harmonization_status": harmonization_status,
                "harmonization_reason": harmonization_reason,
            }
        )

    return output_rows


def build_qc_summary(variant_rows: List[dict], genome_build: str) -> List[dict]:
    total = len(variant_rows)
    missing = sum(row["missing_status"] == "SKIP" for row in variant_rows)
    invalid = sum(row["genotype_status"] == "FAIL" for row in variant_rows)
    duplicates = sum(row["duplicate_status"] == "WARN" for row in variant_rows)
    non_autosomal = sum(row["chromosome_status"] == "WARN" for row in variant_rows)
    invalid_positions = sum(row["position_status"] == "FAIL" for row in variant_rows)
    harmonization_pass = sum(row["harmonization_status"] == "PASS" for row in variant_rows)

    call_rate = 1.0 - (missing / total if total else 0.0)
    harmonization_pass_rate = harmonization_pass / total if total else 0.0

    return [
        {"metric": "created_by", "value": "Mahdi Akbarzadeh", "status": "PASS", "interpretation": "Creator attribution recorded."},
        {"metric": "declared_genome_build", "value": genome_build, "status": "PASS" if genome_build else "FAIL", "interpretation": "Genome build declaration used for reference compatibility checks."},
        {"metric": "variant_count", "value": str(total), "status": "PASS", "interpretation": "Total number of synthetic variant rows processed."},
        {"metric": "call_rate", "value": f"{call_rate:.4f}", "status": "PASS" if call_rate >= 0.95 else "WARN", "interpretation": "Synthetic sample-level call rate estimate."},
        {"metric": "invalid_genotype_count", "value": str(invalid), "status": "PASS" if invalid == 0 else "FAIL", "interpretation": "Invalid genotype encodings should be reviewed."},
        {"metric": "duplicate_marker_rows", "value": str(duplicates), "status": "PASS" if duplicates == 0 else "WARN", "interpretation": "Duplicate marker rows were flagged."},
        {"metric": "non_autosomal_rows", "value": str(non_autosomal), "status": "PASS" if non_autosomal == 0 else "WARN", "interpretation": "Non-autosomal markers are flagged in this autosomal-focused demonstration."},
        {"metric": "invalid_position_rows", "value": str(invalid_positions), "status": "PASS" if invalid_positions == 0 else "FAIL", "interpretation": "Positions must be positive integers."},
        {"metric": "harmonization_pass_rate", "value": f"{harmonization_pass_rate:.4f}", "status": "PASS" if harmonization_pass_rate >= 0.70 else "WARN", "interpretation": "Proportion of rows compatible with the synthetic reference alleles."},
    ]


def build_harmonization_status(variant_rows: List[dict]) -> List[dict]:
    counts = Counter(row["harmonization_status"] for row in variant_rows)
    total = len(variant_rows)
    meanings = {
        "PASS": "Compatible with the synthetic reference table.",
        "WARN": "Reportable only with caution or additional review.",
        "FAIL": "Not harmonization-ready.",
        "SKIP": "Skipped because genotype was missing or not applicable.",
    }
    return [
        {
            "status_label": status,
            "count": str(counts.get(status, 0)),
            "proportion": f"{(counts.get(status, 0) / total if total else 0.0):.4f}",
            "interpretation": meanings[status],
        }
        for status in ["PASS", "WARN", "FAIL", "SKIP"]
    ]


def rows_to_html_table(rows: List[dict]) -> str:
    if not rows:
        return "<p>No rows available.</p>"
    headers = list(rows[0].keys())
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        cells = "".join(f"<td>{html.escape(str(row.get(header, '')))}</td>" for header in headers)
        body_rows.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def render_html_report(out_path: Path, qc_rows: List[dict], harmonization_rows: List[dict]) -> None:
    qc_table = rows_to_html_table(qc_rows)
    harmonization_table = rows_to_html_table(harmonization_rows)
    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Genotype QC and Harmonization Readiness Report</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 40px; color: #222; line-height: 1.6; }}
    header {{ border-bottom: 2px solid #222; margin-bottom: 24px; }}
    h1 {{ margin-bottom: 4px; }}
    h2 {{ border-bottom: 1px solid #ddd; padding-bottom: 4px; margin-top: 28px; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0 24px 0; }}
    th, td {{ border: 1px solid #ddd; padding: 8px; vertical-align: top; }}
    th {{ background: #f2f2f2; text-align: left; }}
    .warning {{ background: #fff3cd; padding: 12px; border: 1px solid #ffeeba; }}
    .small {{ color: #555; font-size: 0.95em; }}
  </style>
</head>
<body>
<header>
  <h1>Genotype QC and Harmonization Readiness Report</h1>
  <p><strong>Created by:</strong> Mahdi Akbarzadeh</p>
  <p class="small">Synthetic demonstration only. Not intended for clinical or consumer-facing genetic interpretation.</p>
</header>
<section>
  <h2>QC Summary</h2>
  {qc_table}
</section>
<section>
  <h2>Harmonization Status Summary</h2>
  {harmonization_table}
</section>
<section class="warning">
  <h2>Limitations</h2>
  <p>This report uses synthetic data only. Real genotype harmonization requires verified genome build, platform-specific marker metadata, allele-frequency-aware strand resolution, and strict audit logging.</p>
</section>
<footer>
  <p class="small">Created by: Mahdi Akbarzadeh</p>
</footer>
</body>
</html>
"""
    out_path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a synthetic genotype QC and harmonization readiness workflow.")
    parser.add_argument("--genotype", required=True, help="Path to synthetic genotype TSV.")
    parser.add_argument("--reference", required=True, help="Path to synthetic reference panel TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--genome-build", required=True, help="Declared genome build, such as GRCh37.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genotype_rows = read_tsv(Path(args.genotype))
    reference_rows = read_tsv(Path(args.reference))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    variant_rows = build_variant_status(genotype_rows, reference_rows, args.genome_build)
    qc_rows = build_qc_summary(variant_rows, args.genome_build)
    harmonization_rows = build_harmonization_status(variant_rows)
    write_tsv(outdir / "variant_status_table.tsv", variant_rows, ["rsid", "chromosome", "position", "genotype", "duplicate_status", "chromosome_status", "position_status", "genotype_status", "missing_status", "harmonization_status", "harmonization_reason"])
    write_tsv(outdir / "qc_summary.tsv", qc_rows, ["metric", "value", "status", "interpretation"])
    write_tsv(outdir / "harmonization_status.tsv", harmonization_rows, ["status_label", "count", "proportion", "interpretation"])
    render_html_report(outdir / "harmonization_readiness_report.html", qc_rows, harmonization_rows)
    print("Created by: Mahdi Akbarzadeh")
    print(f"Output directory: {outdir}")
    print("Created outputs:")
    print("- qc_summary.tsv")
    print("- variant_status_table.tsv")
    print("- harmonization_status.tsv")
    print("- harmonization_readiness_report.html")


if __name__ == "__main__":
    main()
