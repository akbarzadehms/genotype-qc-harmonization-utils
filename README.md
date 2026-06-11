# Genotype QC and Harmonization Utilities

**A lightweight synthetic-data demonstration of genotype input validation, quality control, and harmonization readiness reporting.**

Created by: **Mahdi Akbarzadeh**

---

## Why this project matters

Genomic reports are only as reliable as the genotype data and harmonization layer beneath them.

Before a genotype-derived report is generated, a workflow should check:

- Whether the input schema is valid
- Whether genotype encoding is interpretable
- Whether sample-level and variant-level call rates are acceptable
- Whether duplicate markers exist
- Whether chromosome and position fields are plausible
- Whether genome build is declared
- Whether alleles can be aligned to a reference table
- Whether ambiguous strand variants need special handling
- Whether the input is ready for downstream PRS, PGx, or ancestry modules

This repository provides a minimal, public, synthetic-data demonstration of these checks.

---

## Core workflow

```text
Synthetic genotype input
→ input schema validation
→ genotype QC
→ duplicate marker audit
→ chromosome/position audit
→ allele/reference matching
→ ambiguous strand flagging
→ harmonization readiness status
→ TSV outputs
→ HTML technical report
```

---

## Repository structure

```text
genotype-qc-harmonization-utils/
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── pyproject.toml
├── .gitignore
├── scripts/
│   └── genotype_qc_harmonization_demo.py
├── examples/
│   ├── input/
│   │   ├── synthetic_genotype.tsv
│   │   └── synthetic_reference_panel.tsv
│   └── output/
├── docs/
│   └── qc_harmonization_logic.md
└── tests/
```

---

## Quick start

Run the demonstration workflow:

```bash
python scripts/genotype_qc_harmonization_demo.py \
  --genotype examples/input/synthetic_genotype.tsv \
  --reference examples/input/synthetic_reference_panel.tsv \
  --outdir examples/output \
  --genome-build GRCh37
```

Expected outputs:

```text
examples/output/qc_summary.tsv
examples/output/variant_status_table.tsv
examples/output/harmonization_status.tsv
examples/output/harmonization_readiness_report.html
```

---

## Input format

The synthetic genotype file uses this minimal schema:

| Column | Description |
|---|---|
| `rsid` | Variant identifier |
| `chromosome` | Chromosome label |
| `position` | Base-pair position |
| `genotype` | Diploid genotype encoded as two alleles, such as `AG`, `CC`, or `--` |
| `sample_id` | Synthetic sample identifier |

The synthetic reference table uses this minimal schema:

| Column | Description |
|---|---|
| `rsid` | Variant identifier |
| `chromosome` | Reference chromosome |
| `position` | Reference position |
| `effect_allele` | Demonstration effect allele |
| `other_allele` | Demonstration non-effect allele |
| `genome_build` | Declared reference genome build |

---

## Status labels

| Status | Meaning |
|---|---|
| `PASS` | Check passed and is reportable |
| `WARN` | Check passed with a limitation or caution |
| `FAIL` | Check failed and requires review |
| `SKIP` | Check was not applicable |

---

## Scientific limitations

This repository is a synthetic technical demonstration. It does not provide clinical interpretation, medical advice, consumer genetic interpretation, or production-grade harmonization.

Real-world harmonization requires platform-specific validation, genome-build verification, allele-frequency-aware strand resolution, reference-panel matching, and strict audit logging.

---

## Created by

Mahdi Akbarzadeh
