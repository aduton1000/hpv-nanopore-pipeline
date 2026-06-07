# HPV Nanopore Analysis Pipeline (shared cluster, CPU-only)

End-to-end reconstruction of **HPV whole genomes** from Oxford Nanopore **tiling-amplicon**
reads for **12 high-risk types** (HPV 16/18/31/33/35/39/45/52/56/58/59/66). Basecalling and
demultiplexing happen on the instrument (GridION / PromethION P2 Solo); this pipeline starts
from demultiplexed `fastq_pass/barcodeNN/*.fastq.gz` and produces per-type consensus genomes,
variant calls, coverage QC, and publication-quality figures — **entirely on CPU** (no GPU, no
basecaller).

> **Consumes the design from [`hpv-primer-design`](https://github.com/aduton1000/hpv-primer-design).**
> This pipeline ingests that project's **scheme v5** bundle (per-type `primer.bed` +
> `reference.fasta`, `scheme_summary.tsv`, the 12 canonical-accession references). The design
> repo is the source of the scheme; this repo is the analysis environment + validated pipeline.

## Status: validated (software)
- **End-to-end on the real caller path** — ARTIC **1.10.1** + Clair3 **2.0.1**
  (`engine=artic+clair3`): **100 % genome, 0 % N, origin-seam N = 0** across all sample×type
  combinations, and a HPV16/HPV58 **co-infection** correctly split.
- **All-users deployment verified** — `verify_env.sh` run under a genuine **non-admin**
  account: **53 passed, 0 failed** (named-env activation; deployed pipeline run end-to-end in
  the user's own workdir; no writes to the shared pool/conda; silent profile auto-load; Clair3
  2.x models present/readable).

> This is **software validation**. **Biological validation on control samples (characterised
> isolates / clinical positives) is the next, separate, wet-lab step** and is not claimed here.

## How it works (multi-environment, CPU-only)
Six isolated `hpv-*` conda envs (their tool versions would conflict if merged). The driver
`run_pipeline.py` runs in **`hpv-py`**; stages needing another env are called per-stage via
`conda run -n <env>` (e.g. consensus calls `hpv-artic` for `artic minion` + Clair3). The
origin-spanning junction amplicon (amp3) is handled by a **two-run** scheme (5 non-wrapping
amplicons on the canonical reference + amp3 on a rotated reference, then rotate-back-and-merge).
Full detail in [`PIPELINE_STAGES.md`](PIPELINE_STAGES.md) and [`docs/`](docs/).

## Install → verify → run

**Install (admin, once per cluster):** uploads the design bundle + this repo, creates the six
named `hpv-*` envs in the shared conda, fetches Clair3 models, deploys the pipeline code to
`/hpc/opt/hpv-bioinfo/pipeline`, and installs the auto-loading profile.
```bash
sudo ./setup_hpv_nanopore_env.sh --check  --bundle /path/to/hpv-analysis-bundle   # dry-run
sudo ./setup_hpv_nanopore_env.sh          --bundle /path/to/hpv-analysis-bundle   # install
```

**Verify (any user, any node):**
```bash
./verify_env.sh        # proves env activation, deployed run, no shared-pool writes; PASS/FAIL
```

**Run (any user, zero setup).** The login profile auto-loads `SCHEME_DIR`, `REF_DIR`,
`MODELS_DIR`, `HPV_PIPELINE` and puts `hpv-run` on `PATH` — no `source`, no clone, nothing in
`$HOME`. Outputs land in the user's own `--outdir`, never in `/hpc/opt`:
```bash
cd ~/my_run
hpv-run \
  --fastq-pass  fastq_pass \
  --samplesheet sample_sheet.csv \
  --scheme-dir  "$SCHEME_DIR" \
  --ref         "$REF_DIR/hpv_all_types.fasta" \
  --models-dir  "$MODELS_DIR" \
  --outdir      results \
  --use-artic --require-artic
```
(SLURM/non-login batch jobs must `source /etc/profile.d/hpv-bioinfo.sh` first — see the manual.)

## Shared-deployment model
- **GitHub** (this repo) = source / version control.
- **`/hpc/opt/hpv-bioinfo/pipeline`** (`$HPV_PIPELINE`) = the deployed copy every user runs via
  `hpv-run`. The profile auto-loads at login (vars + PATH only, silent).
- **Updating** = an admin pulls new code and re-runs the deploy step of
  `setup_hpv_nanopore_env.sh`. Ordinary users never touch this repo or `$HOME` to run it.

## Repository layout
```
setup_hpv_nanopore_env.sh     shared, all-users, CPU-only installer (idempotent, sudo-aware)
verify_env.sh                 non-admin verification (PASS/FAIL; the 53/0 run above)
envs/*.yml                    the six hpv-* conda env specs (version floors)
etc/profile.d/hpv-bioinfo.sh  auto-loading login profile (vars + PATH only, silent)
pipeline/                     hpvpipe.py, run_pipeline.py, simulate_test_data.py,
                              Snakefile, main.nf, nextflow.config, testdata/ (SYNTHETIC)
qc_lib.py, qc_report.py       publication QC figure module (placed at $INSTALL_ROOT on deploy)
PIPELINE_STAGES.md            full stage-by-stage map
docs/                         FINAL manual v1.1 (.md/.pdf/.docx) + figures/
docs/installed_versions.txt   captured `conda list` of all six hpv-* envs (provenance of the
                              validated run; complements the envs/*.yml floors)
```

## Tool versions (validated run)
artic 1.10.1 · clair3 2.0.1 · minimap2 2.31-r1302 · samtools/bcftools 1.23.1 · medaka 2.2.2 ·
nanoplot 1.46.2 · flye 2.9.6 · multiqc 1.35 · seqkit 2.13.0 · mosdepth 0.3.14 · bedtools
2.31.1 · chopper 0.13.0 · kraken2 2.17.1 · snakemake 9.22.0 · nextflow 26.04.3 · python 3.13.
Full per-env capture in `docs/installed_versions.txt`.

## What is NOT in this repo
Large / non-redistributable assets are fetched or built by the installer, never stored here:
Clair3 models, the minimap2 `.mmi` index, the conda envs / Miniforge, and the reference data
(consumed from the design bundle). **No real sequencing data or results are ever committed** —
only the small **synthetic** `pipeline/testdata/` and the `sample_sheet_TEMPLATE.csv` pattern
are included (see `.gitignore`).

## Full manual
See [`docs/HPV_Nanopore_Pipeline_Manual.pdf`](docs/HPV_Nanopore_Pipeline_Manual.pdf) (also
`.md`/`.docx`) — architecture, the six envs, every stage, origin-junction handling, Clair3
model selection, running real data, interpreting outputs, the validation record, and
troubleshooting.

## License
This repository is **dual-licensed by file type**:
- **Code → MIT** (`LICENSE-MIT.txt`): all scripts/source — `*.py`, `*.sh`, `Snakefile`,
  `*.nf`, `*.lua`, `*.yml`/`*.yaml`, `nextflow.config` and similar (the bulk of this repo).
- **Data & documents → CC-BY-4.0** (`LICENSE-CC-BY-4.0.txt`): the manual under `docs/`
  (`.md`/`.pdf`/`.docx`) and its figures, `docs/installed_versions.txt`, `PIPELINE_STAGES.md`,
  and the synthetic `pipeline/testdata/`.

**CC-BY attribution:** Bright Adu, ORCID 0000-0002-3622-2137, Genomics and Bioinformatics
Core, NMIMR, University of Ghana, Legon, Ghana. **Cite via `CITATION.cff`** (details in
`LICENSE-NOTICE.md`).

*Ambiguous files (assignment stated):* `envs/*.yml` and `nextflow.config` → **MIT**
(build config); `docs/installed_versions.txt` (captured `conda list`) → **CC-BY-4.0**
(provenance data); `pipeline/testdata/*.fastq.gz` (synthetic) → **CC-BY-4.0**;
`PIPELINE_STAGES.md`/`README.md` → **CC-BY-4.0** (docs). GitHub may show "View license"
rather than one SPDX badge for a dual-licensed repo — expected, not an error.
