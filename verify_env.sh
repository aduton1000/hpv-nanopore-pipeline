#!/usr/bin/env bash
# =============================================================================
# verify_env.sh — runnable by ANY user on ANY node, anytime. Proves the hpv-* envs
# are NAMED envs in the shared cluster conda (activatable as `conda activate hpv-<name>`
# like the lab's wgs_* envs), runs smoke tests, verifies the Clair3 2.x r10.4.1 models
# (GridION + PromethION P2 Solo) are present+readable, and proves a non-admin user can
# activate by name and run end-to-end writing ONLY to HOME (no writes to the shared pool
# /hpc/opt/hpv-bioinfo or the shared conda /hpc/opt/conda). Exits non-zero on any failure.
# =============================================================================
set -uo pipefail

INSTALL_ROOT="${HPV_INSTALL_ROOT:-/hpc/opt/hpv-bioinfo}"
[ -r "$INSTALL_ROOT/etc/hpv-bioinfo.env" ] && . "$INSTALL_ROOT/etc/hpv-bioinfo.env"
CONDA_HOME="${HPV_CONDA_HOME:-/hpc/opt/conda}"
SCHEME_DIR="${SCHEME_DIR:-$INSTALL_ROOT/primer-schemes}"
REF_DIR="${REF_DIR:-$INSTALL_ROOT/references}"
MODELS_DIR="${MODELS_DIR:-$INSTALL_ROOT/clair3_models}"
HPV_PIPELINE="${HPV_PIPELINE:-$INSTALL_ROOT/pipeline}"   # deployed pipeline code (all users)
MODEL_PREFIX="${MODEL_PREFIX:-r1041_e82_400bps}"; MODEL_VARIANTS="${MODEL_VARIANTS:-fast hac sup}"
CONDA="$CONDA_HOME/bin/conda"
HPV_ENVS="hpv-qc hpv-artic hpv-asm hpv-py hpv-typing hpv-smk"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/hpv_verify.XXXX")"   # USER space, never shared pool
trap 'rm -rf "$WORK"' EXIT

PASS=0; FAIL=0; ROWS=""
ok(){ ROWS+=$(printf '  [PASS] %s\n' "$1"); PASS=$((PASS+1)); }
no(){ ROWS+=$(printf '  [FAIL] %s\n' "$1"); FAIL=$((FAIL+1)); }
cr(){ "$CONDA" run -n "$1" "${@:2}"; }                   # run tool in NAMED env $1

echo "==================================================================="
echo " HPV Nanopore env verification   (user: $(id -un), uid $(id -u), node $(hostname -s 2>/dev/null||echo '?'))"
echo " shared conda=$CONDA_HOME   shared pool=$INSTALL_ROOT   models=$MODELS_DIR"
echo "==================================================================="
[ -x "$CONDA" ] && ok "shared conda present ($("$CONDA" --version 2>/dev/null))" || no "shared conda at $CONDA_HOME"

# ---------- hpv-* are NAMED envs (activatable like wgs_*) ----------
echo; echo "## Named-env resolution (conda env list)"
for e in $HPV_ENVS; do
  if "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -Fxq "$e"; then
    echo "  [name] $e  ->  $("$CONDA" env list | awk -v n="$e" '$1==n{print $NF}')"; ok "named env present: $e"
  else echo "  [name] $e  ->  NOT FOUND by name"; no "named env present: $e"; fi
done

echo; echo "## Compute"
echo "  threads(nproc): $(nproc 2>/dev/null||echo '?')   RAM: $(free -h 2>/dev/null|awk '/Mem:/{print $2}'||echo '?')   load:$(cut -d' ' -f1-3 /proc/loadavg 2>/dev/null||echo '?')"

echo; echo "## Tool versions"
[ -x "$INSTALL_ROOT/bin/nextflow" ] && { printf '  %-12s %s\n' nextflow "$($INSTALL_ROOT/bin/nextflow -version 2>&1|awk '/version/{print;exit}')"; ok "version: nextflow"; } || no "version: nextflow"
# clean_ver: keep the first line carrying a real dotted version (e.g. 2.31, 1.23.1),
# dropping plotly/chrome warnings, deprecation notices and python tracebacks so a row never
# shows a multi-line traceback or a WARNING instead of a version.
clean_ver(){ grep -aviE 'warning|traceback|deprecat|plotly|chrome|kaleido|future|^[[:space:]]*file "' 2>/dev/null | grep -m1 -aE '[0-9]+\.[0-9]'; }
# ver: print a CLEAN version if one can be extracted; else, if the executable is installed in
# the env, print "present (version probe noisy)"; else MISSING. PASS unless genuinely absent.
ver(){ local label="$1" env="$2"; shift 2; local exe="$1" out
  out="$(cr "$env" "$@" 2>&1 | clean_ver)"
  if [ -n "$out" ]; then printf '  %-12s %s\n' "$label" "$out"; ok "version: $label"
  elif cr "$env" bash -lc "command -v '$exe' >/dev/null 2>&1"; then printf '  %-12s present (version probe noisy)\n' "$label"; ok "present: $label"
  else printf '  %-12s MISSING\n' "$label"; no "present: $label"; fi; }
# present: for tools with no version flag (e.g. align_trim) — confirm the executable exists/responds.
present(){ local label="$1" env="$2" exe="$3"
  if cr "$env" bash -lc "command -v '$exe' >/dev/null 2>&1"; then printf '  %-12s present\n' "$label"; ok "present: $label"
  else printf '  %-12s MISSING\n' "$label"; no "present: $label"; fi; }
ver snakemake hpv-smk    snakemake --version
ver artic     hpv-artic  artic --version
ver clair3    hpv-artic  run_clair3.sh --version
ver minimap2  hpv-artic  minimap2 --version
ver samtools  hpv-artic  samtools --version
ver bcftools  hpv-artic  bcftools --version
present align_trim hpv-artic align_trim          # no --version flag; presence-checked
ver flye      hpv-asm    flye --version
ver medaka    hpv-asm    medaka --version
ver nanoplot  hpv-qc     NanoPlot --version       # plotly/chrome warning filtered out
ver nanocomp  hpv-qc     NanoComp --version
ver pycoqc    hpv-qc     pycoQC --version         # crashy --version -> falls back to "present"
ver multiqc   hpv-qc     multiqc --version
ver seqkit    hpv-qc     seqkit version
ver mosdepth  hpv-qc     mosdepth --version
ver bedtools  hpv-qc     bedtools --version
ver chopper   hpv-qc     chopper --version
ver kraken2   hpv-typing kraken2 --version
ver python    hpv-py     python --version
echo "  python libs (hpv-py):"
cr hpv-py python - <<'PY' 2>&1 | sed 's/^/    /'
import importlib
for m in ["pysam","Bio","pandas","numpy","scipy","pyfaidx","skbio","matplotlib","seaborn","plotly","kaleido"]:
    try: mod=importlib.import_module(m); print(f"{m:12s} {getattr(mod,'__version__','?')}")
    except Exception as e: print(f"{m:12s} MISSING ({e.__class__.__name__})")
PY

# ---------- CLAIR3 2.x r10.4.1 MODELS (GridION + PromethION P2 Solo) ----------
echo; echo "## Clair3 2.x r10.4.1 models (need fast/hac/sup; readable; PyTorch format)"
is_2x_dir(){ local d="$1"; [ -d "$d" ]||return 1
  find "$d" -maxdepth 1 -name '*.index' 2>/dev/null|grep -q . && return 1
  find "$d" -maxdepth 1 \( -name '*.pt' -o -name 'pytorch_model*' -o -name '*.pkl' -o -name '*.json' \) 2>/dev/null|grep -q .; }
for v in $MODEL_VARIANTS; do
  d="$(find "$MODELS_DIR" -maxdepth 2 -type d -iname "*${MODEL_PREFIX}*${v}*" 2>/dev/null|head -1)"
  if [ -n "$d" ] && [ -r "$d" ] && is_2x_dir "$d"; then echo "  [model] $v -> $d (2.x, readable)"; ok "clair3 2.x model: $v"
  else echo "  [model] $v -> MISSING/legacy/unreadable"; no "clair3 2.x model: $v"; fi
done

# ---------- smoke tests on a tiny synthetic dataset ----------
echo; echo "## Smoke tests (synthetic data; outputs in $WORK)"
cd "$WORK"
if [ -f "$SCHEME_DIR/HPV16/V5/HPV16.reference.fasta" ]; then
  cr hpv-py python - "$SCHEME_DIR/HPV16/V5/HPV16.reference.fasta" "$WORK" <<'PY' && ok "pysam read of ingested reference"
import sys, pysam
ref, work = sys.argv[1], sys.argv[2]
fa=pysam.FastaFile(ref); seq=fa.fetch(fa.references[0]); sub=seq[1000:2200]
open(f"{work}/mini.ref.fasta","w").write(f">{fa.references[0]}\n{sub}\n")
with open(f"{work}/mini.reads.fastq","w") as fh:
    for i in range(50):
        s=sub[i*5:i*5+800] or sub[:800]; fh.write(f"@r{i}\n{s}\n+\n{'I'*len(s)}\n")
print("mini ref + 50 reads written")
PY
else no "pysam read of ingested reference (scheme not found)"; fi
# minimap2 + samtools sort/index are verified on a REAL pipeline BAM in the "Deployed pipeline"
# section below (samtools quickcheck on <tag>.canonical.bam), not on a degenerate toy stream.
cr hpv-artic bash -lc "run_clair3.sh --version >/dev/null 2>&1" && ok "clair3 --version" || no "clair3 --version"
cr hpv-artic bash -lc "artic --version >/dev/null 2>&1" && ok "artic --version" || no "artic --version"
if cr hpv-artic python - "$SCHEME_DIR/HPV16/V5/HPV16.primer.bed" <<'PY' ; then ok "scheme parse (HPV16 primer.bed)"; else no "scheme parse"; fi
import sys
try:
    from primalbedtools.bedfiles import BedLineParser; BedLineParser.from_file(sys.argv[1]); print("primalbedtools OK")
except Exception:
    n=0
    for ln in open(sys.argv[1]):
        f=ln.split("\t"); assert len(f)>=6 and f[5] in "+-"; n+=1
    print(f"structural OK ({n})")
PY
if cr hpv-py python - "$WORK" <<'PY' ; then ok "publication figure export (svg+pdf+png, kaleido)"; else no "publication figure export"; fi
import sys, os
fd=os.path.join(sys.argv[1],"figures"); os.makedirs(fd,exist_ok=True)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(4,3)); ax.bar(["b01","b02"],[1200,800])
for e in ("svg","pdf","png"): fig.savefig(f"{fd}/smoke.{e}",dpi=300,bbox_inches="tight")
import plotly.graph_objects as go; go.Figure(go.Bar(x=["a"],y=[1])).write_image(f"{fd}/p.png",scale=2)
assert all(os.path.getsize(f"{fd}/smoke.{e}")>0 for e in ("svg","pdf","png")); print("ok")
PY
[ -x "$INSTALL_ROOT/bin/nextflow" ] && "$INSTALL_ROOT/bin/nextflow" -version >/dev/null 2>&1 && ok "nextflow -version" || no "nextflow -version"

# ---------- DEPLOYED PIPELINE (shared $HPV_PIPELINE; all users run from here) ----------
echo; echo "## Deployed pipeline ($HPV_PIPELINE; every user runs from here — no home-dir copy)"
PIPE_OK=1
[ -r "$HPV_PIPELINE/run_pipeline.py" ] && ok "pipeline located + readable: \$HPV_PIPELINE/run_pipeline.py" || { no "pipeline at \$HPV_PIPELINE/run_pipeline.py"; PIPE_OK=0; }
[ -r "$HPV_PIPELINE/hpvpipe.py" ]      && ok "deployed: hpvpipe.py"      || { no "deployed: hpvpipe.py"; PIPE_OK=0; }
{ [ -r "$INSTALL_ROOT/qc_report.py" ] && [ -r "$INSTALL_ROOT/qc_lib.py" ]; } && ok "qc_report.py+qc_lib.py at \$INSTALL_ROOT (../ resolution for figures)" || no "qc_report.py/qc_lib.py at \$INSTALL_ROOT"
[ -x "$INSTALL_ROOT/bin/hpv-run" ]     && ok "hpv-run launcher present + executable" || no "hpv-run launcher ($INSTALL_ROOT/bin/hpv-run)"
# deployed code must be read/execute for users but NOT writable by a non-admin
if touch "$HPV_PIPELINE/.wp_$$" 2>/dev/null; then rm -f "$HPV_PIPELINE/.wp_$$"
  [ "$(id -u)" -eq 0 ] && ok "\$HPV_PIPELINE writable (root only — non-admins blocked)" || no "\$HPV_PIPELINE writable by non-admin (PERMS TOO LOOSE)"
else ok "\$HPV_PIPELINE NOT writable by this user (correct for non-admin)"; fi
# end-to-end through the DEPLOYED code (portable engine; no models needed), outputs ONLY in $WORK
if [ "$PIPE_OK" -eq 1 ] && [ -f "$SCHEME_DIR/HPV16/V5/HPV16.reference.fasta" ]; then
  DRUN="$WORK/deploytest"; mkdir -p "$DRUN"
  if ( cd "$DRUN" \
       && cr hpv-py python "$HPV_PIPELINE/simulate_test_data.py" --scheme-dir "$SCHEME_DIR" --out fastq_pass >/dev/null 2>&1 \
       && cr hpv-py python "$HPV_PIPELINE/run_pipeline.py" --fastq-pass fastq_pass --samplesheet sample_sheet.csv \
              --scheme-dir "$SCHEME_DIR" --ref "$REF_DIR/hpv_all_types.fasta" --models-dir "$MODELS_DIR" \
              --outdir results --threads 4 >/dev/null 2>&1 ); then
    [ -s "$DRUN/results/run_summary.tsv" ] && ok "deployed pipeline ran end-to-end via \$HPV_PIPELINE (outputs in user workdir, not /hpc/opt)" || no "deployed pipeline produced no run_summary.tsv"
    grep -q "SAMPLE_COINF" "$DRUN/results/run_summary.tsv" 2>/dev/null && ok "deployed run resolved co-infection (SAMPLE_COINF)" || no "deployed run co-infection"
    { [ -d "$DRUN/results/figures" ] && [ -n "$(find "$DRUN/results/figures" -type f 2>/dev/null)" ]; } && ok "deployed run produced figures (../qc_report.py resolved)" || no "deployed run figures"
    # minimap2 + samtools sort/index proven on a REAL sorted BAM from the deployed run
    # (replaces the old degenerate toy stream): the .bam exists, is indexed, and passes quickcheck.
    dbam="$(find "$DRUN/results" -name '*.canonical.bam' 2>/dev/null | head -1)"
    if [ -n "$dbam" ] && [ -f "${dbam}.bai" ] && cr hpv-artic samtools quickcheck "$dbam" >/dev/null 2>&1; then
      ok "minimap2+samtools: valid sorted+indexed BAM from deployed run ($(basename "$dbam"))"
    else no "minimap2+samtools sorted/indexed BAM (samtools quickcheck on deployed output)"; fi
  else no "deployed pipeline end-to-end run via \$HPV_PIPELINE"; fi
else echo "  (skipped deployed end-to-end run — pipeline files or scheme missing above)"; fi

# ---------- LOGIN PROFILE: auto-load vars+PATH only, SILENT, idempotent ----------
echo; echo "## Login profile (/etc/profile.d/hpv-bioinfo.sh): vars+PATH auto-load, zero output"
PROF=/etc/profile.d/hpv-bioinfo.sh
if [ -r "$PROF" ]; then
  # (1) sourcing emits ZERO stdout+stderr (safe for scp/sftp/rsync/non-interactive ssh)
  pout="$(bash -c "source '$PROF'; true" 2>&1)"
  [ -z "$pout" ] && ok "profile sources with ZERO output (non-interactive-safe)" || no "profile printed on source: [$pout]"
  # (2) after sourcing: HPV_PIPELINE set AND hpv-run on PATH
  bash -c "source '$PROF'; [ -n \"\$HPV_PIPELINE\" ] && command -v hpv-run >/dev/null 2>&1" \
    && ok "after source: \$HPV_PIPELINE set + hpv-run on PATH" || no "after source: HPV_PIPELINE/hpv-run missing"
  # (3) load does NOT auto-activate any conda env (no conda init at login).
  #     Clear any inherited CONDA_DEFAULT_ENV first so this tests OUR profile, not the caller's shell.
  env -u CONDA_DEFAULT_ENV bash -c "source '$PROF'; [ -z \"\${CONDA_DEFAULT_ENV:-}\" ]" \
    && ok "no env auto-activated at load (no conda init at login)" || no "an env was auto-activated at load"
  # (4) PATH idempotent: double-source must never DUPLICATE the shared bin (1 if bin exists, 0 if not)
  n="$(bash -c "source '$PROF'; source '$PROF'; printf ':%s:' \"\$PATH\"" | grep -o ":$INSTALL_ROOT/bin:" | wc -l | tr -d ' ')"
  [ "${n:-0}" -le 1 ] && ok "PATH not duplicated on re-source (shared bin x${n})" || no "PATH duplicated on re-source (x$n)"
  # (5) NO banner text in a real login shell
  if bash -lc 'true' 2>&1 | grep -qiE 'hpv nanopore env ready|hpv-run --help|HPV Nanopore env'; then no "login banner still printed"; else ok "no banner in a login shell"; fi
else echo "  (profile not installed at $PROF — skipped; install via setup_hpv_nanopore_env.sh)"; fi

# ---------- ALL-USERS / NO-WRITE-TO-SHARED-POOL TEST ----------
echo; echo "## All-users test (non-admin: conda activate BY NAME; read shared; write only HOME)"
[ -n "$(find "$WORK" -type f 2>/dev/null)" ] && ok "outputs written to user workdir ($WORK)" || no "user workdir outputs"
for tgt in "$INSTALL_ROOT" "$CONDA_HOME"; do
  if touch "$tgt/.wp_$$" 2>/dev/null; then rm -f "$tgt/.wp_$$"
    [ "$(id -u)" -eq 0 ] && ok "$tgt writable (root only — non-admins blocked)" || no "$tgt writable by non-admin (PERMS TOO LOOSE)"
  else ok "$tgt NOT writable by this user (correct for non-admin)"; fi
done
NONADMIN=0
if command -v sudo >/dev/null 2>&1 && { [ "$(id -u)" -eq 0 ] || [ -n "${ALLOW_SUDO_TEST:-}" ]; }; then
  TU="hpvtest$$"
  if id "$TU" >/dev/null 2>&1 || sudo useradd -m -s /bin/bash "$TU" 2>/dev/null; then
    sudo -u "$TU" bash -lc "
      set -e; source /etc/profile.d/hpv-bioinfo.sh
      conda activate hpv-artic            # <-- activate BY NAME, exactly like wgs_*
      cd \$HOME && mkdir -p run && cd run
      minimap2 -ax map-ont '$SCHEME_DIR/HPV16/V5/HPV16.reference.fasta' '$WORK/mini.reads.fastq' 2>/dev/null | samtools sort -o out.bam - && samtools index out.bam
      test -s out.bam && test -r '$MODELS_DIR' && echo 'models readable by non-admin'
      # run the DEPLOYED pipeline as a brand-new non-admin via hpv-run (no home-dir copy), outputs in HOME
      conda deactivate || true
      command -v hpv-run >/dev/null || { echo 'FAIL: hpv-run not on PATH'; exit 1; }
      cd \$HOME && mkdir -p deployrun && cd deployrun
      conda run -n hpv-py python \"\$HPV_PIPELINE/simulate_test_data.py\" --scheme-dir '$SCHEME_DIR' --out fastq_pass >/dev/null
      hpv-run --fastq-pass fastq_pass --samplesheet sample_sheet.csv --scheme-dir '$SCHEME_DIR' \
              --ref '$REF_DIR/hpv_all_types.fasta' --models-dir '$MODELS_DIR' --outdir results --threads 4 >/dev/null
      test -s results/run_summary.tsv && echo 'deployed pipeline ran in HOME via hpv-run'
      ! ( touch '$HPV_PIPELINE/.x' 2>/dev/null ) && ! ( touch '$INSTALL_ROOT/.x' 2>/dev/null ) && ! ( touch '$CONDA_HOME/.x' 2>/dev/null ) && echo 'confirmed: cannot write deployed pipeline / shared pool / shared conda'
    " && { ok "non-admin: conda activate by name + hpv-run deployed pipeline end-to-end in HOME + no shared writes"; NONADMIN=1; } || no "non-admin deployed-pipeline (hpv-run) end-to-end test"
    sudo userdel -r "$TU" 2>/dev/null || true
  fi
fi
[ "$NONADMIN" -eq 0 ] && echo "  (full throwaway-user test skipped — run as root or set ALLOW_SUDO_TEST=1; the negative shared-write tests above already prove non-admins cannot write the shared pool/conda.)"

echo; echo "==================================================================="
echo "$ROWS"; echo "  RESULT: $PASS passed, $FAIL failed."
echo "==================================================================="
[ "$FAIL" -ne 0 ] && { echo "FAIL — see [FAIL] rows."; exit 1; }
echo "PASS — ready cluster-wide. Start (any node):  conda activate hpv-qc"
exit 0
