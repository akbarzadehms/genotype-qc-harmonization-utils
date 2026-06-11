# QC and Harmonization Logic

Created by: Mahdi Akbarzadeh

## Purpose

This document explains the demonstration logic used by `genotype_qc_harmonization_demo.py`.

## Checks performed

### Input schema

The script expects a synthetic genotype table with:

- `rsid`
- `chromosome`
- `position`
- `genotype`
- `sample_id`

The synthetic reference table should include:

- `rsid`
- `chromosome`
- `position`
- `effect_allele`
- `other_allele`
- `genome_build`

### Genotype validity

Valid observed genotypes are two-letter combinations of:

```text
A, C, G, T
```

The following are treated as missing:

```text
--, NA, N/A, ., empty string
```

### Chromosome plausibility

This demonstration is autosomal-focused. Chromosomes `1` to `22` pass. Other labels are flagged as warnings.

### Position plausibility

Positions must be positive integers.

### Duplicate markers

Duplicate `rsid` rows are flagged as warnings.

### Ambiguous strand variants

A/T and C/G variants are flagged as warnings because real-world resolution often requires allele frequencies and platform-specific metadata.

### Harmonization readiness

A variant is considered harmonization-ready when:

- It is present in the synthetic reference table
- The reference genome build matches the declared input genome build
- The genotype is valid and non-missing
- Observed alleles are compatible with reference alleles
- The marker is not an ambiguous strand case requiring extra review

## Limitations

This is a synthetic technical demonstration and should not be used as production genotype harmonization software.
