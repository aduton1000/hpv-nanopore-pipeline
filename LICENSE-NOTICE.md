# Licensing notice — dual license by file type

This repository is **dual-licensed**, split by file type:

- **Code → MIT** — see [`LICENSE-MIT.txt`](LICENSE-MIT.txt).
  Covers all scripts and source: `*.py`, `*.sh`, `Snakefile`, `*.nf`, `*.lua`,
  `*.yml`/`*.yaml`, `nextflow.config`, and similar executable/source files (the bulk of this
  repository).
- **Data & documents → CC-BY-4.0** — see [`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt).
  Covers the manual under `docs/` (`.md`/`.pdf`/`.docx`) and its figures,
  `docs/installed_versions.txt`, `PIPELINE_STAGES.md`, and the synthetic
  `pipeline/testdata/`.

**CC-BY-4.0 attribution.** Bright Adu, ORCID
[0000-0002-3622-2137](https://orcid.org/0000-0002-3622-2137), Genomics and Bioinformatics
Core, Noguchi Memorial Institute for Medical Research, University of Ghana, Legon, Ghana.
Please **cite via [`CITATION.cff`](CITATION.cff)**.

## Genuinely ambiguous files (assignment stated)
- `envs/*.yml` conda env specs and `nextflow.config` — build / source config; licensed as
  **code (MIT)**.
- `docs/installed_versions.txt` — a captured `conda list` (data / provenance record);
  licensed **CC-BY-4.0**.
- `pipeline/testdata/*.fastq.gz` — synthetic test data; licensed **CC-BY-4.0**.
- `PIPELINE_STAGES.md`, `README.md` and this `LICENSE-NOTICE.md` are documentation
  (**CC-BY-4.0**); `CITATION.cff` and `.gitignore` are repository metadata.

> Note: GitHub may show **"View license"** rather than a single SPDX badge for a
> dual-licensed repository — this is expected, not an error.
