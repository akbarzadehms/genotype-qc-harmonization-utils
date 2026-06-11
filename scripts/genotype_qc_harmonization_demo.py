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


def classify_harmonization(
    row: dict,
    reference_by_rsid: Dict[str, dict],
    genome_build: str,
) -> Tuple[str, str]:
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


def build_variant_status(
    genotype_rows: List[dict],
    reference_rows: List[dict],
    genome_build: str,
) -> List[dict]:
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
        {
            "metric": "created_by",
            "value": "Mahdi Akbarzadeh",
            "status": "PASS",
            "interpretation": "Creator attribution recorded.",
        },
        {
            "metric": "declared_genome_build",
            "value": genome_build,
            "status": "PASS" if genome_build else "FAIL",
            "interpretation": "Genome build declaration used for reference compatibility checks.",
        },
        {
            "metric": "variant_count",
            "value": str(total),
            "status": "PASS",
            "interpretation": "Total number of synthetic variant rows processed.",
        },
        {
            "metric": "call_rate",
            "value": f"{call_rate:.4f}",
            "status": "PASS" if call_rate >= 0.95 else "WARN",
            "interpretation": "Synthetic sample-level call rate estimate.",
        },
        {
            "metric": "invalid_genotype_count",
            "value": str(invalid),
            "status": "PASS" if invalid == 0 else "FAIL",
            "interpretation": "Invalid genotype encodings should be reviewed.",
        },
        {
            "metric": "duplicate_marker_rows",
            "value": str(duplicates),
            "status": "PASS" if duplicates == 0 else "WARN",
            "interpretation": "Duplicate marker rows were flagged.",
        },
        {
            "metric": "non_autosomal_rows",
            "value": str(non_autosomal),
            "status": "PASS" if non_autosomal == 0 else "WARN",
            "interpretation": "Non-autosomal markers are flagged in this autosomal-focused demonstration.",
        },
        {
            "metric": "invalid_position_rows",
            "value": str(invalid_positions),
            "status": "PASS" if invalid_positions == 0 else "FAIL",
            "interpretation": "Positions must be positive integers.",
        },
        {
            "metric": "harmonization_pass_rate",
            "value": f"{harmonization_pass_rate:.4f}",
            "status": "PASS" if harmonization_pass_rate >= 0.70 else "WARN",
            "interpretation": "Proportion of rows compatible with the synthetic reference alleles.",
        },
    ]


def build_harmonization_status(variant_rows: List[dict]) -> List[dict]:
    counts = Counter(row["harmonization_status"] for row in variant_rows)
    total = len(variant_rows)
    return [
        {
            "status_label": status,
            "count": str(counts.get(status, 0)),
            "proportion": f"{(counts.get(status, 0) / total if total else 0.0):.4f}",
            "interpretation": {
                "PASS": "Compatible with the synthetic reference table.",
                "WARN": "Reportable only with caution or additional review.",
                "FAIL": "Not harmonization-ready.",
                "SKIP": "Skipped because genotype was missing or not applicable.",
            }[status],
        }
        for status in ["PASS", "WARN", "FAIL", "SKIP"]
    ]


def status_badge(status: str) -> str:
    status_clean = html.escape(status.upper())
    css_class = {
        "PASS": "pass",
        "WARN": "warn",
        "FAIL": "fail",
        "SKIP": "skip",
        "BLOCK": "fail",
    }.get(status_clean, "skip")
    return f'<span class="badge {css_class}">{status_clean}</span>'


def get_metric(qc_rows: List[dict], metric: str, default: str = "NA") -> str:
    for row in qc_rows:
        if row.get("metric") == metric:
            return row.get("value", default)
    return default


def get_metric_status(qc_rows: List[dict], metric: str, default: str = "SKIP") -> str:
    for row in qc_rows:
        if row.get("metric") == metric:
            return row.get("status", default)
    return default


def rows_to_html_table(rows: List[dict], status_columns: set[str] | None = None) -> str:
    if not rows:
        return "<p>No rows available.</p>"

    status_columns = status_columns or set()
    headers = list(rows[0].keys())
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []

    for row in rows:
        cells = []
        for header in headers:
            value = str(row.get(header, ""))
            if header in status_columns:
                cells.append(f"<td>{status_badge(value)}</td>")
            else:
                cells.append(f"<td>{html.escape(value)}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")

    return f"<table><thead><tr>{header_html}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def render_status_bars(harmonization_rows: List[dict]) -> str:
    rows_html = []
    for row in harmonization_rows:
        status = row["status_label"]
        count = html.escape(row["count"])
        proportion = float(row["proportion"])
        width = max(0.0, min(proportion * 100.0, 100.0))
        css_class = {
            "PASS": "bar-fill-pass",
            "WARN": "bar-fill-warn",
            "FAIL": "bar-fill-fail",
            "SKIP": "bar-fill-skip",
        }.get(status, "bar-fill-skip")

        rows_html.append(
            f"""
            <div class="status-row">
              {status_badge(status)}
              <div class="bar"><div class="{css_class}" style="width: {width:.2f}%"></div></div>
              <strong>{count}</strong>
            </div>
            """
        )
    return "\n".join(rows_html)


def render_html_report(
    out_path: Path,
    qc_rows: List[dict],
    harmonization_rows: List[dict],
    genome_build: str,
) -> None:
    variant_count = get_metric(qc_rows, "variant_count")
    call_rate = float(get_metric(qc_rows, "call_rate", "0"))
    harmonization_pass_rate = float(get_metric(qc_rows, "harmonization_pass_rate", "0"))
    call_rate_status = get_metric_status(qc_rows, "call_rate")
    harmonization_status = get_metric_status(qc_rows, "harmonization_pass_rate")
    invalid_count = int(get_metric(qc_rows, "invalid_genotype_count", "0"))
    invalid_position_count = int(get_metric(qc_rows, "invalid_position_rows", "0"))
    overall_label = "Review" if invalid_count or invalid_position_count else "Ready"

    qc_table = rows_to_html_table(qc_rows, status_columns={"status"})
    harmonization_table = rows_to_html_table(harmonization_rows)
    status_bars = render_status_bars(harmonization_rows)

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Genotype QC and Harmonization Readiness Report</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <style>
    :root {{
      --bg: #f6f8fb;
      --panel: #ffffff;
      --ink: #14213d;
      --muted: #5c677d;
      --line: #d9e2ec;
      --navy: #0b2545;
      --green: #166534;
      --amber: #92400e;
      --red: #991b1b;
      --gray: #374151;
      --green-bg: #dcfce7;
      --amber-bg: #fef3c7;
      --red-bg: #fee2e2;
      --gray-bg: #f3f4f6;
      --shadow: 0 12px 30px rgba(15, 23, 42, 0.08);
    }}

    * {{ box-sizing: border-box; }}

    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
      color: var(--ink);
      background: var(--bg);
      line-height: 1.55;
    }}

    .hero {{
      background: linear-gradient(135deg, #081f3a 0%, #123c69 55%, #1e6091 100%);
      color: white;
      padding: 44px 56px;
    }}

    .hero-inner {{
      max-width: 1180px;
      margin: 0 auto;
    }}

    .eyebrow {{
      letter-spacing: 0.08em;
      text-transform: uppercase;
      font-size: 0.78rem;
      color: #c7d2fe;
      font-weight: 700;
      margin-bottom: 12px;
    }}

    h1 {{
      margin: 0;
      font-size: 2.35rem;
      line-height: 1.15;
      font-weight: 800;
    }}

    .subtitle {{
      max-width: 920px;
      margin: 16px 0 0;
      color: #e0f2fe;
      font-size: 1.04rem;
    }}

    .meta {{
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
      margin-top: 24px;
    }}

    .meta-item {{
      background: rgba(255,255,255,0.12);
      border: 1px solid rgba(255,255,255,0.18);
      border-radius: 999px;
      padding: 7px 12px;
      font-size: 0.9rem;
    }}

    .container {{
      max-width: 1180px;
      margin: -24px auto 48px;
      padding: 0 20px;
    }}

    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 16px;
      margin-bottom: 22px;
    }}

    .card, .section {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 18px;
      box-shadow: var(--shadow);
    }}

    .card {{
      padding: 20px;
    }}

    .metric-label {{
      color: var(--muted);
      font-size: 0.86rem;
      font-weight: 700;
      letter-spacing: 0.02em;
      text-transform: uppercase;
      margin-bottom: 8px;
    }}

    .metric-value {{
      font-size: 1.8rem;
      font-weight: 800;
      margin-bottom: 6px;
    }}

    .metric-note {{
      color: var(--muted);
      font-size: 0.92rem;
    }}

    .layout {{
      display: grid;
      grid-template-columns: 1.2fr 0.8fr;
      gap: 22px;
      align-items: start;
    }}

    .section {{
      padding: 24px;
      margin-bottom: 22px;
    }}

    h2 {{
      margin: 0 0 14px;
      font-size: 1.35rem;
      color: var(--navy);
    }}

    h3 {{
      margin: 18px 0 10px;
      color: var(--navy);
      font-size: 1.05rem;
    }}

    .badge {{
      display: inline-flex;
      align-items: center;
      border-radius: 999px;
      padding: 5px 10px;
      font-size: 0.82rem;
      font-weight: 800;
      white-space: nowrap;
    }}

    .pass {{ color: var(--green); background: var(--green-bg); }}
    .warn {{ color: var(--amber); background: var(--amber-bg); }}
    .fail {{ color: var(--red); background: var(--red-bg); }}
    .skip {{ color: var(--gray); background: var(--gray-bg); }}
    .technical {{ color: #1e3a8a; background: #dbeafe; }}

    .status-stack {{
      display: grid;
      gap: 10px;
    }}

    .status-row {{
      display: grid;
      grid-template-columns: 90px 1fr 80px;
      align-items: center;
      gap: 10px;
    }}

    .bar {{
      height: 11px;
      background: #e5e7eb;
      border-radius: 999px;
      overflow: hidden;
    }}

    .bar-fill-pass {{ height: 100%; background: #22c55e; }}
    .bar-fill-warn {{ height: 100%; background: #f59e0b; }}
    .bar-fill-fail {{ height: 100%; background: #ef4444; }}
    .bar-fill-skip {{ height: 100%; background: #6b7280; }}

    table {{
      border-collapse: collapse;
      width: 100%;
      overflow: hidden;
      border-radius: 12px;
      border: 1px solid var(--line);
      font-size: 0.94rem;
    }}

    th, td {{
      padding: 11px 12px;
      border-bottom: 1px solid var(--line);
      vertical-align: top;
      text-align: left;
    }}

    th {{
      background: #f1f5f9;
      color: #0f172a;
      font-size: 0.82rem;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }}

    tr:last-child td {{ border-bottom: none; }}

    .callout {{
      border-left: 5px solid #f59e0b;
      background: #fffbeb;
      padding: 14px 16px;
      border-radius: 12px;
      color: #78350f;
      margin-top: 12px;
    }}

    .note {{
      border-left: 5px solid #2563eb;
      background: #eff6ff;
      padding: 14px 16px;
      border-radius: 12px;
      color: #1e3a8a;
      margin-top: 12px;
    }}

    .footer {{
      color: var(--muted);
      text-align: center;
      font-size: 0.9rem;
      padding: 14px 0 0;
    }}

    .mono {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, "Liberation Mono", monospace;
      font-size: 0.9em;
    }}

    @media (max-width: 900px) {{
      .summary-grid, .layout {{ grid-template-columns: 1fr; }}
      .hero {{ padding: 34px 24px; }}
      h1 {{ font-size: 1.85rem; }}
    }}
  </style>
</head>
<body>
  <header class="hero">
    <div class="hero-inner">
      <div class="eyebrow">Synthetic genomic reporting demonstration</div>
      <h1>Genotype QC and Harmonization Readiness Report</h1>
      <p class="subtitle">
        A technical demonstration report for genotype input validation, quality control,
        allele/reference compatibility checks, and downstream report-readiness assessment.
      </p>
      <div class="meta">
        <div class="meta-item"><strong>Created by:</strong> Mahdi Akbarzadeh</div>
        <div class="meta-item"><strong>Input type:</strong> Synthetic genotype table</div>
        <div class="meta-item"><strong>Declared build:</strong> {html.escape(genome_build)}</div>
        <div class="meta-item"><strong>Use:</strong> Technical demonstration only</div>
      </div>
    </div>
  </header>

  <main class="container">
    <section class="summary-grid">
      <div class="card">
        <div class="metric-label">Variant rows</div>
        <div class="metric-value">{html.escape(variant_count)}</div>
        <div class="metric-note">Synthetic records processed</div>
      </div>
      <div class="card">
        <div class="metric-label">Call rate</div>
        <div class="metric-value">{call_rate * 100:.2f}%</div>
        <div class="metric-note">{status_badge(call_rate_status)}</div>
      </div>
      <div class="card">
        <div class="metric-label">Harmonization pass rate</div>
        <div class="metric-value">{harmonization_pass_rate * 100:.2f}%</div>
        <div class="metric-note">{status_badge(harmonization_status)}</div>
      </div>
      <div class="card">
        <div class="metric-label">Overall readiness</div>
        <div class="metric-value">{overall_label}</div>
        <div class="metric-note"><span class="badge technical">TECHNICAL ONLY</span></div>
      </div>
    </section>

    <section class="layout">
      <div>
        <section class="section">
          <h2>Executive Technical Summary</h2>
          <p>
            This synthetic example demonstrates how a genomic reporting system can audit genotype-level
            data before downstream PRS, pharmacogenomics, ancestry, or clinical-facing modules are rendered.
          </p>
          <div class="callout">
            <strong>Interpretation:</strong> The demonstration input contains intentionally inserted QC issues.
            This is expected and is useful for showing report-readiness logic.
          </div>
        </section>

        <section class="section">
          <h2>QC Summary</h2>
          {qc_table}
        </section>
      </div>

      <aside>
        <section class="section">
          <h2>Harmonization Status</h2>
          <div class="status-stack">
            {status_bars}
          </div>

          <h3>Readiness interpretation</h3>
          <p>
            The input is suitable for technical demonstration, but not for production or clinical reporting.
            Failed and warning-level records should be resolved before downstream analysis.
          </p>
        </section>

        <section class="section">
          <h2>Report Controls</h2>
          <table>
            <tbody>
              <tr><th>Clinical use</th><td><span class="badge fail">Not allowed</span></td></tr>
              <tr><th>Customer-facing use</th><td><span class="badge fail">Not allowed</span></td></tr>
              <tr><th>Technical demo</th><td><span class="badge pass">Allowed</span></td></tr>
              <tr><th>Data type</th><td>Synthetic only</td></tr>
            </tbody>
          </table>
        </section>
      </aside>
    </section>

    <section class="section">
      <h2>Harmonization Status Summary</h2>
      {harmonization_table}
      <div class="note">
        <strong>Design note:</strong> This report is generated by the Python workflow from validated TSV outputs,
        rather than being manually edited as a static HTML artifact.
      </div>
    </section>

    <section class="section">
      <h2>Limitations</h2>
      <ul>
        <li>This report uses synthetic demonstration data only.</li>
        <li>It does not provide medical, clinical, pharmacogenomic, ancestry, or consumer genetic interpretation.</li>
        <li>Real genotype harmonization requires verified genome build, platform-specific marker metadata, allele-frequency-aware strand resolution, and audit logging.</li>
        <li>Downstream PRS, PGx, and ancestry modules should not be rendered unless QC and harmonization readiness are explicitly documented.</li>
      </ul>
    </section>

    <div class="footer">
      Created by: Mahdi Akbarzadeh · Genotype QC and Harmonization Utilities · Synthetic demonstration report
    </div>
  </main>
</body>
</html>
"""
    out_path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a synthetic genotype QC and harmonization readiness workflow."
    )
    parser.add_argument("--genotype", required=True, help="Path to synthetic genotype TSV.")
    parser.add_argument("--reference", required=True, help="Path to synthetic reference panel TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--genome-build", required=True, help="Declared genome build, such as GRCh37.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    genotype_path = Path(args.genotype)
    reference_path = Path(args.reference)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genotype_rows = read_tsv(genotype_path)
    reference_rows = read_tsv(reference_path)

    variant_rows = build_variant_status(
        genotype_rows=genotype_rows,
        reference_rows=reference_rows,
        genome_build=args.genome_build,
    )
    qc_rows = build_qc_summary(variant_rows=variant_rows, genome_build=args.genome_build)
    harmonization_rows = build_harmonization_status(variant_rows=variant_rows)

    write_tsv(
        outdir / "variant_status_table.tsv",
        variant_rows,
        [
            "rsid",
            "chromosome",
            "position",
            "genotype",
            "duplicate_status",
            "chromosome_status",
            "position_status",
            "genotype_status",
            "missing_status",
            "harmonization_status",
            "harmonization_reason",
        ],
    )
    write_tsv(
        outdir / "qc_summary.tsv",
        qc_rows,
        ["metric", "value", "status", "interpretation"],
    )
    write_tsv(
        outdir / "harmonization_status.tsv",
        harmonization_rows,
        ["status_label", "count", "proportion", "interpretation"],
    )
    render_html_report(
        out_path=outdir / "harmonization_readiness_report.html",
        qc_rows=qc_rows,
        harmonization_rows=harmonization_rows,
        genome_build=args.genome_build,
    )

    print("Created by: Mahdi Akbarzadeh")
    print(f"Output directory: {outdir}")
    print("Created outputs:")
    print("- qc_summary.tsv")
    print("- variant_status_table.tsv")
    print("- harmonization_status.tsv")
    print("- harmonization_readiness_report.html")


if __name__ == "__main__":
    main()
