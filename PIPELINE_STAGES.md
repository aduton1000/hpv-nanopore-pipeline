# PIPELINE_STAGES.md — HPV Nanopore tiling-amplicon end-to-end stage map

Entry point: demultiplexed, basecalled **`fastq_pass/barcodeNN/*.fastq.gz`** from the
instrument (GridION / PromethION P2 Solo — basecalling/demux already done there).
CPU-only throughout. All **core stages are implemented** in `pipeline/hpvpipe.py` and
driven by three equivalent runners: `pipeline/run_pipeline.py`, `pipeline/Snakefile`,
`pipeline/main.nf`. Validated end-to-end on synthetic reads (incl. origin-spanning amp3
and a co-infection): see "End-to-end test" below.

## Run the pipeline (real data)
The pipeline code is **deployed on the shared pool** at `$HPV_PIPELINE`
(`/hpc/opt/hpv-bioinfo/pipeline`) — every user runs it from there; no home-dir copy. Run
from your own working dir: every output lands there; nothing writes to `/hpc/opt`.
*(Source vs deployed: GitHub = source/version control; `$HPV_PIPELINE` = the deployed copy
users run; updating = an admin re-runs `setup_hpv_nanopore_env.sh`, refreshing it in place.)*

The profile **auto-loads at login** (vars + PATH only, silent) — interactive users need no
manual `source`. **SLURM / non-login batch jobs** don't read `/etc/profile.d`, so a job
script must `source /etc/profile.d/hpv-bioinfo.sh` (or `…/etc/hpv-bioinfo.env`) first.
```bash
# interactive login: SCHEME_DIR/REF_DIR/MODELS_DIR/HPV_PIPELINE/HPV_ARTIC_RUN already set, hpv-run on PATH
# batch job script only: source /etc/profile.d/hpv-bioinfo.sh
cd ~/run_2026-06-06                        # your own workdir

# (A) hpv-run — recommended (runs run_pipeline.py inside hpv-py via conda run; no activation)
hpv-run --fastq-pass /path/to/fastq_pass --samplesheet sample_sheet.csv \
    --scheme-dir "$SCHEME_DIR" --ref "$REF_DIR/hpv_all_types.fasta" --models-dir "$MODELS_DIR" \
    --outdir results --threads 16 --use-artic        # drop --use-artic for the portable caller
#   explicit equivalent:  conda activate hpv-py && python "$HPV_PIPELINE/run_pipeline.py" …

# (B) Snakemake (DAG-parallel across samples/types)
conda activate hpv-py
snakemake -s "$HPV_PIPELINE/Snakefile" -j 16 --config scheme_dir=$SCHEME_DIR \
    ref=$REF_DIR/hpv_all_types.fasta models_dir=$MODELS_DIR \
    fastq_pass=/path/to/fastq_pass samplesheet=sample_sheet.csv use_artic=true

# (C) Nextflow
nextflow run "$HPV_PIPELINE/main.nf" --fastq_pass /path/to/fastq_pass --samplesheet sample_sheet.csv \
    --scheme_dir $SCHEME_DIR --ref $REF_DIR/hpv_all_types.fasta --models $MODELS_DIR --use_artic true
```
`--use-artic` runs `artic minion` + Clair3 2.x for the canonical amplicons (model
auto-selected per the FASTQ `basecall_model_version_id`); the rotated-reference junction
fix (amp3) is layered on in all cases. Without it, a portable minimap2+samtools consensus
engine produces the identical merged output (no Clair3) — useful for quick checks.

### sample_sheet.csv
```
barcode,sample_id,expected_types,notes
barcode01,SAMPLE_001,,                       # expected_types optional (;-separated)
barcode02,SAMPLE_002,HPV16;HPV58,co-infection
```
`barcode` = MinKNOW `fastq_pass/barcodeNN/` folder; `sample_id` = your label.

## Stage map (all core stages IMPLEMENTED)

| # | Stage | Status | Env | What it does |
|---|-------|--------|-----|--------------|
| 1 | **Ingest** barcoded FASTQ | **IMPL** | hpv-qc | Concatenate `fastq_pass/barcodeNN/*.fastq.gz` per barcode; map barcode→sample_id via `sample_sheet.csv`. |
| 2 | **Read QC + figures** | **IMPL** | hpv-qc, hpv-py | per-read lengths/stats + `qc_report.py`: per-barcode reads/yield, length & quality dists, length-window marking 2–4 kb AND the **HPV52 1.1 kb** amplicon. (NanoPlot/NanoComp optional add-on.) |
| 3 | **Length/quality filter** | **IMPL** | hpv-qc | `chopper --minlength 900 --maxlength 4800 -q 9` (portable awk fallback). **900 bp lower bound keeps the HPV52 gap-fill amplicon.** |
| 4 | **HPV type assignment (binning)** | **IMPL** | hpv-py | `hpvpipe.py bin`: competitive `minimap2 -ax map-ont` vs `hpv_all_types.fasta`, best-hit (AS + identity ≥0.80) per read, tally per type, **split FASTQ per type**, types-present table (≥`min-reads`, co-infection aware). |
| 5 | **Per-type consensus (+junction fix)** | **IMPL** | hpv-artic / hpv-py | `hpvpipe.py consensus`: `artic minion` + Clair3 2.x (`--use-artic`) or portable engine, **+ rotated-reference amp3 recovery (Task 1)** → merged canonical consensus FASTA, VCF, trimmed BAM, per-amplicon coverage, `junction_proof.tsv`. |
| 6 | **(optional) De novo path** | optional | hpv-asm | `flye` + `medaka` polish for divergent/novel/co-infection cases. Not in the core path. |
| 7 | **Per-amplicon + per-pool coverage QC** | **IMPL** | hpv-py | per-amplicon depth (pysam over `amplicon.bed`, mosdepth-compatible) → `qc_report.py` **per-amplicon coverage heatmap** (dropout), depth tracks, pool-balance. |
| 8 | **Consensus QC + summary + figures** | **IMPL** | hpv-py | %genome/%N per (sample × type) completeness heatmap; type-assignment + run-summary tables; publication figures. (MultiQC combined report when `multiqc` present.) |
| 9 | **Lineage assignment** (per consensus) | follow-up | hpv-py | A/B/C/D sublineage (esp. **HPV58 sublineage C**) via Nextclade custom dataset/alignment — the one remaining optional add-on. |

## Model selection (Clair3 2.x / ARTIC) — GridION AND PromethION P2 Solo
The lab sequences on **both a GridION and a PromethION P2 Solo** (both **R10.4.1**
chemistry, but the Dorado version / basecall preset can differ between runs). The shared
`$MODELS_DIR` therefore holds the **full r10.4.1 Clair3 2.x (PyTorch) set — `fast`, `hac`
and `sup`** (`r1041_e82_400bps_{fast,hac,sup}_*`). Clair3 ≥2.0.1 uses the **new PyTorch
model format**; legacy 1.x (TensorFlow) models will fail at runtime, so the installer
verifies the 2.x format and FAILS LOUDLY if any of fast/hac/sup is missing/legacy.

- **Auto-select (preferred):** GridION *and* P2-Solo FASTQ headers carry
  `basecall_model_version_id`, so ARTIC picks the matching model automatically — just pass
  the directory:
  ```bash
  artic minion --model-dir "$MODELS_DIR" --scheme-directory "$SCHEME_DIR" \
        --scheme-name HPV58 --scheme-version V5 --threads 8 <reads> SAMPLE_HPV58
  ```
- **List what's available:**
  ```bash
  ls -1 "$MODELS_DIR"     # e.g. r1041_e82_400bps_fast_v500  _hac_v500  _sup_v500
  ```
- **Explicit `--model` (when a header lacks `basecall_model_version_id`, or to force one):**
  ```bash
  # GridION run basecalled with the SUP preset:
  artic minion --model r1041_e82_400bps_sup_v500  --model-dir "$MODELS_DIR" ...
  # PromethION P2 Solo run basecalled with the HAC preset:
  artic minion --model r1041_e82_400bps_hac_v500  --model-dir "$MODELS_DIR" ...
  # (use the fast model only for quick/low-accuracy checks)
  ```
  Match the `--model` name to the run's basecaller preset; use the exact directory names
  printed by `ls "$MODELS_DIR"` (version suffix may differ as ONT updates Rerio).

**Hard rule:** models are **read** from `$MODELS_DIR`; pass `--model-dir "$MODELS_DIR"` so
nothing is written under `$CONDA_PREFIX/bin/models` (read-only for non-admins). No runtime
writes to the shared install (`/hpc/opt/...`).

## ARTIC 1.10.1 scheme loading (offline) — how `--use-artic` resolves the scheme
ARTIC 1.10.1's `get_scheme()` is **remote-first** (fetches `quick-lab/primerschemes/index.json`)
and a compute node with no internet hits *"Failed to fetch manifest"*. Two things make it
work offline:
1. **The pipeline calls `artic minion --bed <primer.bed> --ref <reference.fasta>`** — ARTIC
   uses these directly (`if args.bed and args.ref: …`) and **never touches the manifest**.
   This is the path `hpvpipe.py consensus` uses; it is the guaranteed-offline route.
2. The bundle also ships the **ARTIC 1.10.1 layout** so the `--scheme-name/--scheme-length/--scheme-version`
   form resolves offline too: `primer-schemes/HPV{t}/2000/v5.0.0/{primer.bed,reference.fasta}`
   + a local `primer-schemes/index.json` (get_scheme falls back to it, then skips download
   because the files already exist) + legacy `HPV{t}/V5/HPV{t}.scheme.bed`.
   Example: `artic minion --scheme-directory $SCHEME_DIR --scheme-name HPV16 --scheme-length 2000 --scheme-version v5.0.0 …`
- **chrom fix:** reference/BED chrom is `HPV{t}_{accession}` (underscore). A pipe (`HPV16|K02718`)
  is **rejected** by primalbedtools/ARTIC (`chrom ^[A-Za-z0-9_.]+$`); this was a second cause
  of `--use-artic` failing.
- **Upgrading an already-installed scheme tree in place** (no full re-ingest):
  ```bash
  python /hpc/opt/hpv-bioinfo/.../convert_scheme_to_artic.py --scheme-dir $SCHEME_DIR        # --check for dry-run
  ```

## Validating the real ARTIC + Clair3 2.x path (cluster)
The portable engine is for quick checks; **validate the actual Clair3 path with strict mode**
(no silent fallback — a broken artic path then fails loudly instead of looking green):
```bash
source /etc/profile.d/hpv-bioinfo.sh             # sets HPV_PIPELINE + HPV_ARTIC_RUN
conda activate hpv-py
python "$HPV_PIPELINE/simulate_test_data.py" --scheme-dir "$SCHEME_DIR" --out fastq_pass   # 400bps headers
hpv-run --fastq-pass fastq_pass --samplesheet sample_sheet.csv \
    --scheme-dir "$SCHEME_DIR" --ref "$REF_DIR/hpv_all_types.fasta" --models-dir "$MODELS_DIR" \
    --outdir results --threads 16 --use-artic --require-artic
```
Each (sample × type) runs artic **twice** (main = amp3-excluded canonical; junction =
rotated amp3-only) — see the origin-junction section — so neither `align_trim` call sees a
negative length. `--require-artic` (strict) **aborts non-zero** if *either* artic run fails;
it does **not** fall back. The run prints a per-(sample × type) table — **engine** (must be
`artic+clair3`), **%genome**, **%N**, **seam N** (0 = junction recovered through the real
Clair3 path), and the **model** auto-selected from the 400bps header; `junction_proof.tsv`
records `junction_engine=artic+clair3`. Production runs that prefer resilience can add
`--allow-fallback` (artic failure → portable engine with a WARN; the engine column says so).

**Expected on the synthetic data:** `SAMPLE_COINF` → HPV16;HPV58 (co-infection); every type
`engine=artic+clair3`, `junction_engine=artic+clair3`, `%genome≈100`, `seam N=0`, no
`negative dimensions` error.

> **Off-cluster status:** the negative-length **geometry fix is verified here** (both beds
> all-positive) and the **two-run + merge orchestration is validated with a faithful artic
> stub** (engine/junction_engine=artic+clair3, seam N=0). Real **Clair3 variant calling**
> can only run where `hpv-artic` (clair3 2.0.1, linux) is installed — **run the command above
> on node1** to confirm the real Clair3 consensus. If it still errored it would now abort
> loudly (strict), not silently fall back.

## Origin-junction amplicon (amp3) — IMPLEMENTED (two-frame ARTIC integration)
Amp3 spans the circular origin, so in the linear primer.bed its right coord < left coord →
**negative amplicon length**. ARTIC 1.10.1 `align_trim` does `np.zeros(amp.length)` and
crashes ("negative dimensions are not allowed"); a single rotated frame only *moves* the
wrap to whichever amplicon spans the new origin, so it does not help. The fix runs artic
**twice, on frames where every amplicon align_trim sees is positive-length**, and merges:

1. **Main run** — `artic minion --bed <canonical bed with amp3 REMOVED> --ref <canonical>`:
   align_trim sees only the 5 non-wrapping amplicons (all positive) → no crash; Clair3 calls
   them. The origin seam is left as N.
2. **Junction run** — `artic minion --bed <rotated amp3-ONLY bed> --ref <rotated reference>`:
   on the rotated frame amp3 is a single **contiguous, positive-length** amplicon (offset
   ≈ L/2, cut kept outside amp3), so align_trim handles it normally and **Clair3 calls amp3**.
3. **Merge** — rotate the amp3 consensus back (`pos_canon = (pos_rot + offset) mod L`) and
   fill only the seam positions still N (no double-count, length unchanged).

So **both runs go through real artic/Clair3** and neither ever sees a negative length. (When
artic is unavailable and `--allow-fallback` is set, the junction run uses a portable rotated
pileup instead; `junction_engine` in the QC records which was used.) The amp3-excluded bed
is written per run as `<sample>.<type>.no_amp3.primer.bed`; the rotated scheme is
`<sample>.<type>.rot/`. Every run emits `<sample>.<type>.junction_proof.tsv`
(`junction_engine`, `canonical_only_seam_N`, `merged_seam_N` =0 when recovered,
`seam_identity_to_reference_pct`, `continuous_across_origin`, `length_unchanged`) and the
per-amplicon coverage shows amp3 covered. Example offsets: HPV16 3953, HPV58 3912, HPV52 3971.

**Geometry verified** (off-cluster): the amp3-excluded canonical bed and the rotated
amp3-only bed are all positive-length for every type → `align_trim` cannot hit the
negative-length crash. The two-run + merge orchestration was validated end-to-end with a
faithful artic stub (engine=artic+clair3, junction_engine=artic+clair3, seam N=0, co-infection
split) — see "Validating the real ARTIC + Clair3 path" for the on-cluster command that
exercises real Clair3.

## End-to-end test (synthetic; proves the junction fix + co-infection)
```bash
source /etc/profile.d/hpv-bioinfo.sh
conda activate hpv-py
python "$HPV_PIPELINE/simulate_test_data.py" --scheme-dir "$SCHEME_DIR" --out fastq_pass   # 3 barcodes incl co-infection + origin-spanning amp3 reads
hpv-run --fastq-pass fastq_pass --samplesheet sample_sheet.csv \
       --scheme-dir "$SCHEME_DIR" --ref "$REF_DIR/hpv_all_types.fasta" --outdir results
```
Verified result: SAMPLE_COINF → **HPV16;HPV58** (co-infection detected); every type 100%
genome; junction proof e.g. HPV58 — canonical-only left **102 N** at the origin seam,
rotated ref (offset 3912) filled 102 bp → **merged seam N=0, 100% identity to reference,
continuous across origin**; 10 publication figure groups generated. Reproduced identically
by the Snakemake runner.

## CPU / threading on the 72-core node
- The pipeline is embarrassingly parallel over (sample × type). Prefer **many samples in
  parallel, modest threads each** over maxing threads on one. With 72 cores, e.g.
  `snakemake -j 72` (or `run_pipeline.py --threads 8` × ~8 concurrent samples).
- Clair3 (`--use-artic`) is the heavy step; ~8 threads/sample is a good balance. minimap2
  binning/mapping and mosdepth/pysam coverage are light (2–4 threads).
- The rotated-reference junction fix roughly **doubles amp3's mapping cost only** (a second
  small map of the binned reads) — negligible vs Clair3.
- Memory: mapping/QC light; Clair3 and (optional) Flye are the peaks — keep concurrency
  within RAM. All scratch/work/caches go to the user's CWD/HOME (`NXF_HOME`,
  `APPTAINER_CACHEDIR`, `CONDA_PKGS_DIRS`), never the shared pool.

## Outputs (all in the user's working dir)
```
results/reads/      per-sample concatenated FASTQ
results/qc/         read_stats.tsv, lengths/      (+ results/figures publication QC)
results/filtered/   length/quality-filtered FASTQ
results/typing/     <sample>.types.tsv, type_counts.tsv, per-type FASTQ, <sample>.present.txt
results/consensus/  <sample>.<type>.{consensus.fasta, vcf.gz, canonical.bam, consensus_qc.tsv,
                    junction_proof.tsv}; consensus_qc.tsv (aggregate)
results/coverage/   <sample>.<type>.{amplicon.regions.bed.gz, per-base.bed.gz}
results/figures/    SVG+PDF+300dpi PNG (per-amplicon heatmap, depth tracks, type bars, …)
results/run_summary.tsv   types present per sample
```
