#!/usr/bin/env bash
# =============================================================================
# setup_hpv_nanopore_env.sh
# Shared, ALL-USERS, CPU-ONLY Nanopore HPV tiling-amplicon analysis environment
# for an HPC cluster. Idempotent, sudo-aware, logged. Basecalling/demux happen on
# the instrument (GridION / PromethION P2 Solo) — this installs the ANALYSIS env only.
#
# Follows the lab convention: hpv-* envs are NAMED envs in the shared cluster conda
# (/hpc/opt/conda/envs), usable cluster-wide via `conda activate hpv-<name>` exactly like
# the existing wgs_* envs. Reference data / models / schemes / Nextflow live on the shared
# pool (/hpc/opt/hpv-bioinfo) so they are identical from every node.
#
# Usage:
#   sudo ./setup_hpv_nanopore_env.sh --bundle /path/to/hpv-analysis-bundle [--check]
#   ./setup_hpv_nanopore_env.sh --check --bundle ...        # dry-run, no changes
#   sudo ./setup_hpv_nanopore_env.sh --bundle ... --force-recreate   # rebuild hpv-* envs ONLY
#
# After install, ANY user on ANY node:  conda activate hpv-qc   (or: source the profile +
#   hpv-activate qc). Never reinstalls conda; never touches base or non-hpv- envs (wgs_* etc).
# =============================================================================
set -euo pipefail

# ---------- CONFIG (override via environment if desired) ----------
INSTALL_ROOT="${INSTALL_ROOT:-/hpc/opt/hpv-bioinfo}"   # shared pool (NOT node-local /opt)
SHARED_GROUP="${SHARED_GROUP:-bioinfo}"
CONDA_BASE="${CONDA_BASE:-}"                            # optional explicit shared-conda base
CONDA_ROOT="${CONDA_ROOT:-$INSTALL_ROOT/miniforge3}"   # only used if NO conda found anywhere
EXPECT_CONDA="${EXPECT_CONDA:-/hpc/opt/conda}"          # expected shared cluster conda
OWN_CONDARC="$INSTALL_ROOT/condarc"                     # project-local condarc (base untouched)
PKG_CACHE="$INSTALL_ROOT/pkgs"                          # shared, world-readable pkgs cache
HPV_ENV_FILE="$INSTALL_ROOT/etc/hpv-bioinfo.env"        # resolved paths for profile/verify
SCHEME_DIR="${SCHEME_DIR:-$INSTALL_ROOT/primer-schemes}"
REF_DIR="${REF_DIR:-$INSTALL_ROOT/references}"
MODELS_DIR="${MODELS_DIR:-$INSTALL_ROOT/clair3_models}"
BIN_DIR="$INSTALL_ROOT/bin"; LOG_DIR="$INSTALL_ROOT/logs"; LOCK_DIR="$INSTALL_ROOT/envs_locks"
PIPELINE_DIR="$INSTALL_ROOT/pipeline"; DOCS_DIR="$INSTALL_ROOT/docs"   # deployed pipeline CODE + manual
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; ENV_YML_DIR="$SCRIPT_DIR/envs"
PIPELINE_SRC="$SCRIPT_DIR/pipeline"   # source tree of the pipeline code
BUNDLE_DIR=""; DRYRUN=0; FORCE_RECREATE=0
HPV_ENVS="hpv-qc hpv-artic hpv-asm hpv-py hpv-typing hpv-smk"   # the ONLY envs we may touch
MAMBA_MIN_MAJ=23; MAMBA_MIN_MIN=10
MF_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
MODEL_PREFIX="${MODEL_PREFIX:-r1041_e82_400bps}"; MODEL_VARIANTS="${MODEL_VARIANTS:-fast hac sup}"

usage(){ sed -n '2,21p' "$0"; exit "${1:-0}"; }
while [ $# -gt 0 ]; do case "$1" in
  --bundle) BUNDLE_DIR="${2:-}"; shift 2;;
  --check|--dry-run) DRYRUN=1; shift;;
  --force-recreate) FORCE_RECREATE=1; shift;;
  -h|--help) usage 0;;
  *) echo "unknown arg: $1" >&2; usage 2;;
esac; done

TS="$(date +%Y%m%d_%H%M%S)"
if [ "$DRYRUN" -eq 0 ]; then TMP_LOG="$(mktemp /tmp/hpv_setup_${TS}.XXXX.log)"; exec > >(tee -a "$TMP_LOG") 2>&1; fi
log(){ printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
sect(){ printf '\n========== %s ==========\n' "$*"; }
warnx(){ printf '[WARN] %s\n' "$*"; }
die(){ printf '[FATAL] %s\n' "$*" >&2; exit 1; }
if [ "$(id -u)" -eq 0 ]; then SUDO=""; else SUDO="sudo"; fi
run(){ if [ "$DRYRUN" -eq 1 ]; then echo "  DRYRUN> $*"; else eval "$@"; fi; }
have(){ command -v "$1" >/dev/null 2>&1; }
umask 022

# ---------- PREFLIGHT ----------
sect "PREFLIGHT"
log "Mode: $([ $DRYRUN -eq 1 ] && echo 'DRY-RUN (no changes)' || echo 'INSTALL')$([ $FORCE_RECREATE -eq 1 ] && echo ' +force-recreate(hpv- envs only)')"
log "CPU-ONLY: no GPU required/checked/used (Clair3 on CPU). Cluster: shared pool /hpc."
[ -r /etc/os-release ] && { . /etc/os-release; log "OS: ${PRETTY_NAME:-unknown}"; }
case "${VERSION_ID:-}" in 20.04|22.04|24.04) log "Supported Ubuntu LTS (${VERSION_ID}).";;
  *) warnx "Not a known Ubuntu LTS ('${VERSION_ID:-?}'). Continuing, untested.";; esac
[ "$DRYRUN" -eq 0 ] && [ -n "$SUDO" ] && { sudo -v || die "sudo required for system steps."; }
[ -n "$BUNDLE_DIR" ] || die "--bundle /path/to/hpv-analysis-bundle is required."
[ -d "$BUNDLE_DIR" ] || die "bundle dir not found: $BUNDLE_DIR"
# confirm the shared pool parent of INSTALL_ROOT is present + writable (sudo)
SHARED_PARENT="$(dirname "$INSTALL_ROOT")"
if [ "$DRYRUN" -eq 0 ]; then
  [ -d "$SHARED_PARENT" ] || die "shared pool $SHARED_PARENT not found (is /hpc mounted on this node?)."
  $SUDO test -w "$SHARED_PARENT" || $SUDO mkdir -p "$INSTALL_ROOT" 2>/dev/null || die "cannot write $SHARED_PARENT with sudo."
  log "shared pool $SHARED_PARENT present and writable (sudo)."
else echo "  DRYRUN> verify shared pool $SHARED_PARENT exists and is sudo-writable"; fi
log "CPU threads: $(nproc 2>/dev/null||echo '?'); RAM: $(free -h 2>/dev/null|awk '/Mem:/{print $2}'||echo '?')"

sect "BUNDLE INTEGRITY (SHA256)"
[ -f "$BUNDLE_DIR/SHA256SUMS" ] || die "No SHA256SUMS in bundle — refusing unverified inputs."
( cd "$BUNDLE_DIR" && sha256sum -c SHA256SUMS >/tmp/hpv_sha_${TS}.out 2>&1 ) \
  && log "Bundle checksums OK ($(wc -l <"$BUNDLE_DIR/SHA256SUMS") files)." \
  || { cat /tmp/hpv_sha_${TS}.out; die "Bundle SHA256 mismatch — aborting before ingest."; }

# ---------- GROUP + INSTALL ROOT on shared pool (setgid 2775) ----------
sect "GROUP + INSTALL ROOT (shared pool)"
getent group "$SHARED_GROUP" >/dev/null && log "group '$SHARED_GROUP' exists." \
  || { log "creating group '$SHARED_GROUP'"; run "$SUDO groupadd '$SHARED_GROUP'"; }
[ -d "$INSTALL_ROOT" ] && log "$INSTALL_ROOT exists." \
  || { log "creating $INSTALL_ROOT (root:$SHARED_GROUP 2775) on shared pool"; run "$SUDO mkdir -p '$INSTALL_ROOT'"; }
run "$SUDO chown root:'$SHARED_GROUP' '$INSTALL_ROOT'"; run "$SUDO chmod 2775 '$INSTALL_ROOT'"
for d in "$BIN_DIR" "$LOG_DIR" "$LOCK_DIR" "$SCHEME_DIR" "$REF_DIR" "$MODELS_DIR" \
         "$PKG_CACHE" "$INSTALL_ROOT/etc"; do
  run "$SUDO install -d -o root -g '$SHARED_GROUP' -m 2775 '$d'"; done
if [ "$DRYRUN" -eq 0 ]; then run "$SUDO cp '$TMP_LOG' '$LOG_DIR/setup_${TS}.log'"; run "$SUDO chmod 0644 '$LOG_DIR/setup_${TS}.log'"; log "Logging to $LOG_DIR/setup_${TS}.log"; fi

# ---------- C. SYSTEM PACKAGES (apt) ----------
sect "APT PACKAGES"
APT_PKGS="build-essential git curl wget unzip pigz tabix bzip2 libbz2-dev liblzma-dev \
zlib1g-dev libcurl4-openssl-dev libssl-dev default-jre parallel datamash fonts-dejavu"
if [ "$DRYRUN" -eq 0 ]; then run "$SUDO apt-get update -y"; run "$SUDO DEBIAN_FRONTEND=noninteractive apt-get install -y $APT_PKGS";
else echo "  DRYRUN> apt-get install -y $APT_PKGS"; fi
if have apptainer; then log "apptainer present ($(apptainer --version 2>/dev/null))."
else log "installing apptainer (PPA; continue on failure)"
  if [ "$DRYRUN" -eq 0 ]; then ( run "$SUDO add-apt-repository -y ppa:apptainer/ppa" && run "$SUDO apt-get update -y" && run "$SUDO apt-get install -y apptainer" ) || warnx "apptainer install failed — continuing.";
  else echo "  DRYRUN> add-apt-repository ppa:apptainer/ppa && apt-get install -y apptainer"; fi
fi

# ============================================================================
# DETECT & REUSE the shared cluster conda (never reinstall; never touch base/others)
# ============================================================================
sect "CONDA: DETECT & REUSE shared cluster conda"
detect_conda(){
  local d b
  if [ -n "$CONDA_BASE" ] && [ -x "$CONDA_BASE/bin/conda" ]; then echo "$CONDA_BASE"; return 0; fi
  if [ -x "$EXPECT_CONDA/bin/conda" ]; then echo "$EXPECT_CONDA"; return 0; fi
  if have conda; then b="$(conda info --base 2>/dev/null || true)"; [ -n "$b" ] && [ -x "$b/bin/conda" ] && { echo "$b"; return 0; }; fi
  if [ -n "${CONDA_PREFIX:-}" ]; then b="${CONDA_PREFIX%%/envs/*}"; [ -x "$b/bin/conda" ] && { echo "$b"; return 0; }; fi
  for d in /hpc/opt/conda /opt/miniforge3 /opt/conda /opt/miniconda3 /usr/local/miniforge3 "${HOME:-/root}/miniforge3" "$CONDA_ROOT"; do
    [ -x "$d/bin/conda" ] && { echo "$d"; return 0; }; done
  return 1
}
conda_ge_min(){ local v maj min; v="$("$1/bin/conda" --version 2>/dev/null|awk '{print $2}')"; [ -n "$v" ]||return 1
  maj="${v%%.*}"; min="${v#*.}"; min="${min%%.*}"
  [ "${maj:-0}" -gt "$MAMBA_MIN_MAJ" ] 2>/dev/null && return 0
  [ "${maj:-0}" -eq "$MAMBA_MIN_MAJ" ] 2>/dev/null && [ "${min:-0}" -ge "$MAMBA_MIN_MIN" ] 2>/dev/null && return 0
  return 1; }

CONDA_HOME=""; FOUND_BASE="$(detect_conda || true)"
if [ -n "$FOUND_BASE" ]; then
  CONDA_HOME="$FOUND_BASE"; cver="$("$CONDA_HOME/bin/conda" --version 2>/dev/null||echo '?')"
  log "REUSING shared conda at $CONDA_HOME  ($cver)"
  [ "$DRYRUN" -eq 1 ] && echo "  DRYRUN> REUSING shared conda at $CONDA_HOME (NO Miniforge download)"
  [ "$CONDA_HOME" = "$EXPECT_CONDA" ] || warnx "shared conda is $CONDA_HOME (expected $EXPECT_CONDA) — proceeding with detected base."
  if [ "$DRYRUN" -eq 0 ] && ! conda_ge_min "$CONDA_HOME"; then
    die "Shared conda ($cver) < ${MAMBA_MIN_MAJ}.${MAMBA_MIN_MIN} (libmamba). Refusing to upgrade the
         cluster install. Ask the cluster admin to update it, or set CONDA_BASE."
  fi
else
  log "no conda found anywhere — would install Miniforge to $CONDA_ROOT (shared pool)"
  [ "$DRYRUN" -eq 1 ] && echo "  DRYRUN> no conda found — would install Miniforge to $CONDA_ROOT (download-then-run; batch -b -p)"
  if [ -x "$CONDA_ROOT/bin/conda" ]; then log "$CONDA_ROOT already has conda — reuse it."
  elif [ -e "$CONDA_ROOT" ]; then die "$CONDA_ROOT exists but no conda binary — refusing to clobber."
  else
    if [ "$DRYRUN" -eq 0 ]; then
      curl -fsSL -o /tmp/Miniforge3_${TS}.sh "$MF_URL"
      $SUDO bash /tmp/Miniforge3_${TS}.sh -b -p "$CONDA_ROOT"   # -b/-p are INSTALLER flags
      rm -f /tmp/Miniforge3_${TS}.sh
    fi
  fi
  CONDA_HOME="$CONDA_ROOT"
fi
CONDA="$CONDA_HOME/bin/conda"; MAMBA="$CONDA_HOME/bin/mamba"; [ -x "$MAMBA" ] || MAMBA="$CONDA"
# named envs land in the shared conda's envs dir (like wgs_*)
SHARED_ENVS_DIR="$CONDA_HOME/envs"
log "hpv-* will be NAMED envs in the shared conda's envs dir: $SHARED_ENVS_DIR"
log "  -> activatable cluster-wide as 'conda activate hpv-<name>', exactly like wgs_* envs."
[ "$DRYRUN" -eq 1 ] && echo "  DRYRUN> will create NAMED envs in $SHARED_ENVS_DIR ; base env and all non-hpv- envs (wgs_* etc.) are NEVER modified/updated/removed/solved-into."

# ---------- OUR condarc (channels/pkgs only; never edit the shared base condarc) ----------
sect "CONDARC (project-local; shared base condarc untouched)"
if [ "$DRYRUN" -eq 0 ]; then
  $SUDO tee "$OWN_CONDARC" >/dev/null <<EOF
channels: [conda-forge, bioconda]
channel_priority: strict
pkgs_dirs: [$PKG_CACHE]
auto_activate_base: false
always_yes: true
EOF
  run "$SUDO chmod 0644 '$OWN_CONDARC'"
else echo "  DRYRUN> write $OWN_CONDARC (channels+pkgs only; shared base condarc untouched; used via CONDARC= for hpv- creates only)"; fi
# run conda/mamba with OUR condarc, as root (install-time writes into the shared conda envs dir)
CR(){ if [ "$DRYRUN" -eq 1 ]; then echo "  DRYRUN> CONDARC=$OWN_CONDARC ${*}"; else $SUDO env CONDARC="$OWN_CONDARC" "$@"; fi; }
ENVRUN(){ local e="$1"; shift; CR "$CONDA" run -n "$e" "$@"; }   # run tool in NAMED env

# helper: does a NAMED env exist (activatable by name)?
env_exists(){ "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -Fxq "$1"; }

# ============================================================================
# CREATE hpv-* AS NAMED ENVS (name-scoped safety; never touch base/wgs_*/non-hpv-)
# ============================================================================
sect "NAMED hpv-* ENVS in $SHARED_ENVS_DIR"
export_lock(){ # $1=name
  [ "$DRYRUN" -eq 1 ] && return 0
  env_exists "$1" || return 0
  $SUDO env CONDARC="$OWN_CONDARC" "$CONDA" env export -n "$1" --no-builds > "$LOCK_DIR/$1.lock.yml" 2>/dev/null || true
  $SUDO env CONDARC="$OWN_CONDARC" "$CONDA" env export -n "$1"            > "$LOCK_DIR/$1.fullbuild.lock.yml" 2>/dev/null || true
  cp "$LOCK_DIR/$1.lock.yml" "$ENV_YML_DIR/$1.lock.yml" 2>/dev/null || true
  log "  exported $LOCK_DIR/$1.lock.yml"
}
create_env(){ # idempotent NAMED create; ONLY ever the six hpv-* names
  local name="$1" yml="$ENV_YML_DIR/$1.yml"
  case " $HPV_ENVS " in *" $name "*) ;; *) die "refusing: '$name' is not one of the allowed hpv-* envs";; esac
  [ -f "$yml" ] || { warnx "missing $yml — skipping $name"; return 0; }
  if env_exists "$name"; then
    if [ "$FORCE_RECREATE" -eq 1 ]; then
      case "$name" in hpv-*) ;; *) die "refusing to remove non-hpv- env '$name'";; esac   # belt-and-braces
      log "--force-recreate: removing ONLY named env $name then recreating"
      run "$SUDO env CONDARC='$OWN_CONDARC' '$CONDA' env remove -y -n '$name'"
    else
      log "env $name already exists in shared conda — SKIP (never auto-overwrite; --force-recreate to rebuild)"
      export_lock "$name"; return 0
    fi
  fi
  log "creating NAMED env $name in $SHARED_ENVS_DIR  (conda activate $name)"
  if [ "$DRYRUN" -eq 0 ]; then
    $SUDO env CONDARC="$OWN_CONDARC" "$MAMBA" env create -y -n "$name" -f "$yml"
    # all-users readability for THIS env only (name-scoped; never chmod other envs)
    [ -d "$SHARED_ENVS_DIR/$name" ] && $SUDO chmod -R a+rX "$SHARED_ENVS_DIR/$name"
  else echo "  DRYRUN> CONDARC=$OWN_CONDARC mamba env create -n $name -f $yml   (-> $SHARED_ENVS_DIR/$name)"; fi
  export_lock "$name"
}
[ "$DRYRUN" -eq 1 ] && echo "  DRYRUN> NO existing environment is touched; only these NEW named envs are created: $HPV_ENVS"
for e in $HPV_ENVS; do create_env "$e"; done

# ---------- ARTIC recipe cross-check (clair3 2.x) ----------
sect "ARTIC RECIPE CROSS-CHECK"
OFFICIAL_YML_URL="https://raw.githubusercontent.com/artic-network/fieldbioinformatics/master/environment.yml"
if [ "$DRYRUN" -eq 0 ]; then
  if curl -fsSL --max-time 30 -o /tmp/artic_official_${TS}.yml "$OFFICIAL_YML_URL" 2>/dev/null; then
    log "official ARTIC deps that hpv-artic.yml omits:"
    comm -23 <(grep -oE '^[[:space:]]*-[[:space:]]*[a-zA-Z0-9_.-]+' /tmp/artic_official_${TS}.yml|sed 's/[- ]//g'|sort -u) \
             <(grep -oE '^[[:space:]]*-[[:space:]]*[a-zA-Z0-9_.-]+' "$ENV_YML_DIR/hpv-artic.yml"|sed 's/[- ]//g'|sort -u) \
      | sed 's/^/    (official-only) /' || true
  else warnx "could not download official ARTIC environment.yml — skipped cross-check."; fi
else echo "  DRYRUN> curl official ARTIC environment.yml and diff deps vs hpv-artic.yml"; fi
grep -qE 'clair3[[:space:]]*>=[[:space:]]*2' "$ENV_YML_DIR/hpv-artic.yml" \
  && log "  confirmed: hpv-artic.yml pins clair3 >=2.0.1 (no <2.0 ceiling)." || warnx "clair3 ceiling missing in hpv-artic.yml!"

# ============================================================================
# CLAIR3 2.x MODELS for BOTH GridION and PromethION P2 Solo (UNCHANGED logic)
# ============================================================================
sect "CLAIR3 2.x MODELS (r10.4.1 fast/hac/sup; GridION + P2 Solo) -> $MODELS_DIR"
list_models(){ find "$MODELS_DIR" -maxdepth 2 -mindepth 1 -type d 2>/dev/null|sed 's#^#    #'||true; }
is_2x_dir(){ local d="$1"; [ -d "$d" ]||return 1
  find "$d" -maxdepth 1 -name '*.index' 2>/dev/null|grep -q . && return 1
  find "$d" -maxdepth 1 \( -name '*.pt' -o -name 'pytorch_model*' -o -name '*.pkl' -o -name '*.json' \) 2>/dev/null|grep -q .; }
have_variant(){ local v="$1" d; while IFS= read -r d; do [ -r "$d" ] && is_2x_dir "$d" && return 0; done \
  < <(find "$MODELS_DIR" -maxdepth 2 -type d -iname "*${MODEL_PREFIX}*${v}*" 2>/dev/null); return 1; }
models_complete(){ local v; for v in $MODEL_VARIANTS; do have_variant "$v"||return 1; done; return 0; }
if [ "$DRYRUN" -eq 1 ]; then
  echo "  DRYRUN> artic_get_models --model-dir $MODELS_DIR  (r10.4.1 Clair3 2.x set)"
  echo "  DRYRUN> LIST $MODELS_DIR and VERIFY 2.x (PyTorch) format for: $MODEL_VARIANTS; Rerio fallback; FAIL LOUDLY if incomplete"
elif models_complete; then log "Clair3 2.x models [$MODEL_VARIANTS] already present & readable — skip."; list_models
else
  log "fetching Clair3 2.x r10.4.1 model set via artic_get_models"
  ENVRUN hpv-artic artic_get_models --model-dir "$MODELS_DIR" 2>/dev/null || warnx "artic_get_models incomplete; will verify/fallback."
  log "models after artic_get_models:"; list_models
  if ! models_complete; then
    warnx "incomplete/wrong-format; trying ONT Rerio r10.4.1 (Clair3 2.x)."
    tmpr="$(mktemp -d /tmp/rerio_${TS}.XXXX)"
    ( git clone --depth 1 https://github.com/nanoporetech/rerio "$tmpr" \
      && ENVRUN hpv-artic python "$tmpr/download_model.py" --clair3 \
      && $SUDO bash -c "cp -r '$tmpr'/clair3_models/*${MODEL_PREFIX}* '$MODELS_DIR'/ 2>/dev/null || cp -r '$tmpr'/clair3_models/* '$MODELS_DIR'/" ) || warnx "Rerio fallback failed."
    rm -rf "$tmpr"; log "models after Rerio:"; list_models
  fi
  run "$SUDO chown -R root:'$SHARED_GROUP' '$MODELS_DIR'"
  run "$SUDO find '$MODELS_DIR' -type d -exec chmod 2775 {} +"
  run "$SUDO find '$MODELS_DIR' -type f -exec chmod 0444 {} +"
  if [ "$DRYRUN" -eq 0 ]; then
    miss=""; for v in $MODEL_VARIANTS; do have_variant "$v"||miss="$miss $v"; done
    [ -n "$miss" ] && die "Required Clair3 2.x r10.4.1 variant(s) missing/legacy/unreadable:$miss . Both GridION and P2 Solo need fast/hac/sup."
    for v in $MODEL_VARIANTS; do d="$(find "$MODELS_DIR" -maxdepth 2 -type d -iname "*${MODEL_PREFIX}*${v}*" 2>/dev/null|head -1)"; log "  [OK] $v model: $d (2.x PyTorch, readable)"; done
  fi
fi

# ---------- NEXTFLOW (shared pool bin) ----------
sect "NEXTFLOW"
if [ -x "$BIN_DIR/nextflow" ]; then log "nextflow present in $BIN_DIR."
else log "installing Nextflow -> $BIN_DIR"
  if [ "$DRYRUN" -eq 0 ]; then tmpd="$(mktemp -d /tmp/nf.XXXX)"
    ( cd "$tmpd" && curl -s https://get.nextflow.io | bash ) \
      && run "$SUDO install -o root -g '$SHARED_GROUP' -m 0755 '$tmpd/nextflow' '$BIN_DIR/nextflow'" || warnx "Nextflow install failed."; rm -rf "$tmpd"
  else echo "  DRYRUN> curl https://get.nextflow.io | bash; install to $BIN_DIR/nextflow"; fi
fi

# ---------- INGEST THE UPLOADED BUNDLE (shared pool) ----------
sect "INGEST BUNDLE -> $INSTALL_ROOT"
run "$SUDO cp -r '$BUNDLE_DIR/primer-schemes/.' '$SCHEME_DIR/'"
run "$SUDO cp -r '$BUNDLE_DIR/references/.' '$REF_DIR/'"
run "$SUDO cp '$BUNDLE_DIR/metadata/scheme_summary.tsv' '$INSTALL_ROOT/' 2>/dev/null || true"
if [ -f "$REF_DIR/hpv_all_types.mmi" ]; then log "hpv_all_types.mmi present."
else log "building minimap2 index (map-ont)"
  if [ "$DRYRUN" -eq 0 ]; then ENVRUN hpv-qc minimap2 -x map-ont -d "$REF_DIR/hpv_all_types.mmi" "$REF_DIR/hpv_all_types.fasta" || warnx "minimap2 index failed."; ENVRUN hpv-qc samtools faidx "$REF_DIR/hpv_all_types.fasta" || true
  else echo "  DRYRUN> minimap2 -x map-ont -d hpv_all_types.mmi hpv_all_types.fasta"; fi
fi
log "re-validating one scheme (HPV16/V5) post-ingest"
if [ "$DRYRUN" -eq 0 ]; then
  ENVRUN hpv-artic python - "$SCHEME_DIR/HPV16/V5/HPV16.primer.bed" "$SCHEME_DIR/HPV16/V5/HPV16.reference.fasta" <<'PYV' || warnx "scheme validation failed."
import sys; bed,ref=sys.argv[1],sys.argv[2]
try:
    from primalbedtools.bedfiles import BedLineParser; BedLineParser.from_file(bed); print("primalbedtools parsed",bed)
except Exception:
    L=sum(len(l.strip()) for l in open(ref) if not l.startswith(">")); n=0
    for ln in open(bed):
        f=ln.split("\t"); assert len(f)>=6 and f[5] in "+-"; s,e=int(f[1]),int(f[2]); assert 0<=s<e<=L; n+=1
    print(f"structural OK: {n} rows, ref {L} bp")
PYV
fi

# ============================================================================
# DEPLOY PIPELINE CODE to the shared pool so EVERY user can run it (not just the
# person whose HOME holds server-setup/). Idempotent: re-running REFRESHES the copy.
# HARD RULE preserved: this is read/execute-only for users; run_pipeline.py writes only
# to its --outdir (user CWD/HOME). Nothing here writes back into /hpc/opt at run time.
# ============================================================================
sect "DEPLOY PIPELINE CODE -> $PIPELINE_DIR"
run "$SUDO install -d -o root -g '$SHARED_GROUP' -m 2775 '$PIPELINE_DIR'"
run "$SUDO install -d -o root -g '$SHARED_GROUP' -m 2775 '$PIPELINE_DIR/testdata'"
# resolve a source file from the pipeline/ subdir first, else the server-setup/ root
pipe_src(){ if [ -f "$PIPELINE_SRC/$1" ]; then echo "$PIPELINE_SRC/$1";
            elif [ -f "$SCRIPT_DIR/$1" ]; then echo "$SCRIPT_DIR/$1"; else echo ""; fi; }
PIPE_SCRIPTS="run_pipeline.py hpvpipe.py qc_report.py qc_lib.py simulate_test_data.py"
PIPE_FILES="Snakefile main.nf nextflow.config PIPELINE_STAGES.md"
for f in $PIPE_SCRIPTS; do s="$(pipe_src "$f")"
  [ -n "$s" ] || { warnx "pipeline source missing: $f — skipping"; continue; }
  run "$SUDO install -o root -g '$SHARED_GROUP' -m 0755 '$s' '$PIPELINE_DIR/$f'"; done   # a+rx
for f in $PIPE_FILES; do s="$(pipe_src "$f")"
  [ -n "$s" ] || { warnx "pipeline source missing: $f — skipping"; continue; }
  run "$SUDO install -o root -g '$SHARED_GROUP' -m 0644 '$s' '$PIPELINE_DIR/$f'"; done    # a+r
[ -d "$PIPELINE_SRC/testdata" ] && run "$SUDO cp -r '$PIPELINE_SRC/testdata/.' '$PIPELINE_DIR/testdata/'"
# run_pipeline.py / Snakefile call ../qc_report.py with PYTHONPATH=../ (i.e. HERE.parent).
# With HERE=$PIPELINE_DIR that parent is $INSTALL_ROOT, so qc_report.py + qc_lib.py must ALSO
# live there for the unchanged relative path to resolve (pipeline logic is NOT modified).
for f in qc_report.py qc_lib.py; do s="$(pipe_src "$f")"; [ -n "$s" ] || continue
  run "$SUDO install -o root -g '$SHARED_GROUP' -m 0755 '$s' '$INSTALL_ROOT/$f'"; done
log "pipeline code deployed to $PIPELINE_DIR (qc_report.py/qc_lib.py also at $INSTALL_ROOT for ../ resolution)."

# ---------- DEPLOY DOCS (manual alongside the code) ----------
sect "DEPLOY DOCS -> $DOCS_DIR"
if [ -d "$SCRIPT_DIR/docs" ]; then
  run "$SUDO install -d -o root -g '$SHARED_GROUP' -m 2775 '$DOCS_DIR'"
  run "$SUDO cp -r '$SCRIPT_DIR/docs/.' '$DOCS_DIR/'"
  log "docs deployed to $DOCS_DIR (manual: HPV_Nanopore_Pipeline_Manual.{md,pdf,docx})"
else warnx "no docs/ in $SCRIPT_DIR — skipping docs deploy."; fi

# ---------- hpv-run LAUNCHER (on PATH for ALL users, any node) ----------
sect "hpv-run LAUNCHER -> $BIN_DIR/hpv-run"
if [ "$DRYRUN" -eq 0 ]; then
  $SUDO tee "$BIN_DIR/hpv-run" >/dev/null <<EOF
#!/usr/bin/env bash
# hpv-run — run the shared HPV pipeline from any node WITHOUT knowing absolute paths or
# pre-activating an env. Generated by setup_hpv_nanopore_env.sh. Runs run_pipeline.py inside
# the hpv-py env via 'conda run', and sets HPV_ARTIC_RUN so --use-artic reaches hpv-artic.
# All outputs go to the user's --outdir (default ./results); nothing is written to /hpc/opt.
#   usage:  hpv-run --fastq-pass fastq_pass --samplesheet sample_sheet.csv \\
#                   --scheme-dir "\$SCHEME_DIR" --ref "\$REF_DIR/hpv_all_types.fasta" \\
#                   --models-dir "\$MODELS_DIR" --outdir results --threads 16 --use-artic --require-artic
set -euo pipefail
HPV_INSTALL_ROOT="${INSTALL_ROOT}"
[ -r "\$HPV_INSTALL_ROOT/etc/hpv-bioinfo.env" ] && . "\$HPV_INSTALL_ROOT/etc/hpv-bioinfo.env"
: "\${HPV_CONDA_HOME:=${CONDA_HOME}}"
: "\${HPV_PIPELINE:=${PIPELINE_DIR}}"
export SCHEME_DIR="\${SCHEME_DIR:-${SCHEME_DIR}}" REF_DIR="\${REF_DIR:-${REF_DIR}}" MODELS_DIR="\${MODELS_DIR:-${MODELS_DIR}}"
export HPV_ARTIC_RUN="\${HPV_ARTIC_RUN:-\$HPV_CONDA_HOME/bin/conda run --no-capture-output -n hpv-artic}"
exec "\$HPV_CONDA_HOME/bin/conda" run --no-capture-output -n hpv-py \\
     python "\$HPV_PIPELINE/run_pipeline.py" "\$@"
EOF
  run "$SUDO chown root:'$SHARED_GROUP' '$BIN_DIR/hpv-run'"; run "$SUDO chmod 0755 '$BIN_DIR/hpv-run'"
  log "installed launcher: hpv-run (on PATH via $BIN_DIR) -> conda run -n hpv-py python \$HPV_PIPELINE/run_pipeline.py"
else echo "  DRYRUN> write $BIN_DIR/hpv-run (conda run -n hpv-py python \$HPV_PIPELINE/run_pipeline.py \"\$@\")"; fi

# ---------- ALL-USERS ACTIVATION (named envs) ----------
sect "ALL-USERS ACTIVATION (named envs)"
if [ "$DRYRUN" -eq 0 ]; then
  $SUDO tee "$HPV_ENV_FILE" >/dev/null <<EOF
# generated by setup_hpv_nanopore_env.sh ($TS) — shared conda + named envs
HPV_CONDA_HOME="$CONDA_HOME"
HPV_CONDA_SH="$CONDA_HOME/etc/profile.d/conda.sh"
HPV_ENV_MODE="named"
HPV_PIPELINE="$PIPELINE_DIR"
SCHEME_DIR="$SCHEME_DIR"
REF_DIR="$REF_DIR"
MODELS_DIR="$MODELS_DIR"
EOF
  run "$SUDO chmod 0644 '$HPV_ENV_FILE'"
else echo "  DRYRUN> write $HPV_ENV_FILE (HPV_CONDA_HOME=$CONDA_HOME, HPV_ENV_MODE=named)"; fi
run "$SUDO install -o root -g root -m 0644 '$SCRIPT_DIR/etc/profile.d/hpv-bioinfo.sh' /etc/profile.d/hpv-bioinfo.sh"
log "installed /etc/profile.d/hpv-bioinfo.sh (root:root 0644): AUTO-LOADS vars + PATH at login,"
log "  SILENT (no banner, safe in non-interactive shells), no conda init/activation."

# ---------- PERMISSIONS HARDENING (all users read/exec; only admins write) ----------
sect "PERMISSIONS HARDENING (all-users read/execute)"
if [ "$DRYRUN" -eq 0 ]; then
  run "$SUDO chgrp -R '$SHARED_GROUP' '$INSTALL_ROOT'"
  run "$SUDO find '$INSTALL_ROOT' -type d -exec chmod 2775 {} +"
  run "$SUDO find '$INSTALL_ROOT' -type f -exec chmod a+r {} +"
  run "$SUDO find '$INSTALL_ROOT' -type f -perm -u+x -exec chmod a+rx {} +"
  run "$SUDO find '$MODELS_DIR' -type f -exec chmod 0444 {} + 2>/dev/null || true"
  [ -d "$PKG_CACHE" ] && run "$SUDO chmod -R a+rX '$PKG_CACHE'"
  # name-scoped: ensure OUR named envs are world-readable; do NOT touch other envs
  for e in $HPV_ENVS; do [ -d "$SHARED_ENVS_DIR/$e" ] && run "$SUDO chmod -R a+rX '$SHARED_ENVS_DIR/$e'"; done
else echo "  DRYRUN> chgrp/chmod $INSTALL_ROOT (2775 dirs, world-read files, models 0444); chmod a+rX ONLY $SHARED_ENVS_DIR/hpv-* (other envs untouched)"; fi
log "hpv-* are NAMED envs in the shared conda ($SHARED_ENVS_DIR) — base and wgs_*/non-hpv- envs NOT touched."
log "HARD RULE: no analysis runtime step writes to $INSTALL_ROOT or the shared conda; caches/work/outputs -> user HOME/CWD."

# ---------- SUMMARY ----------
sect "SUMMARY"
PASS=0; FAIL=0
check(){ if eval "$2" >/dev/null 2>&1; then echo "  [PASS] $1"; PASS=$((PASS+1)); else echo "  [FAIL] $1"; FAIL=$((FAIL+1)); fi; }
if [ "$DRYRUN" -eq 0 ]; then
  check "shared conda ($CONDA_HOME)" "[ -x '$CONDA' ]"
  for e in $HPV_ENVS; do check "named env $e (conda activate $e)" "'$CONDA' env list | awk '{print \$1}' | grep -Fxq '$e'"; done
  for v in $MODEL_VARIANTS; do check "clair3 2.x model: $v" "find '$MODELS_DIR' -maxdepth 2 -type d -iname '*${MODEL_PREFIX}*${v}*' | grep -q ."; done
  check "Nextflow" "[ -x '$BIN_DIR/nextflow' ]"
  check "schemes ingested" "[ -f '$SCHEME_DIR/HPV16/V5/HPV16.primer.bed' ]"
  check "reference index" "[ -f '$REF_DIR/hpv_all_types.mmi' ]"
  check "pipeline code deployed" "[ -f '$PIPELINE_DIR/run_pipeline.py' ] && [ -f '$PIPELINE_DIR/hpvpipe.py' ]"
  check "qc_report at \$INSTALL_ROOT (../ resolution)" "[ -f '$INSTALL_ROOT/qc_report.py' ] && [ -f '$INSTALL_ROOT/qc_lib.py' ]"
  check "hpv-run launcher on PATH" "[ -x '$BIN_DIR/hpv-run' ]"
  check "docs deployed" "[ -f '$DOCS_DIR/HPV_Nanopore_Pipeline_Manual.md' ]"
  check "profile.d + env file" "[ -f /etc/profile.d/hpv-bioinfo.sh ] && [ -f '$HPV_ENV_FILE' ]"
  check "INSTALL_ROOT on shared pool, setgid" "[ \"\$(stat -c '%a' '$INSTALL_ROOT')\" = '2775' ]"
  echo; echo "  RESULT: $PASS passed, $FAIL failed."
else echo "  DRY-RUN complete — no changes made."; fi
echo; echo "Next: run  ./verify_env.sh   (any user, any node) to prove conda activate hpv-* by name, the deployed pipeline runs from \$HPV_PIPELINE, + no shared-pool writes + models readable."
echo "Any user on any node — the profile AUTO-LOADS at login (vars + PATH only, silent):"
echo "  hpv-run --help                          # run the pipeline (SCHEME_DIR/REF_DIR/MODELS_DIR/HPV_PIPELINE already set)"
echo "  conda activate hpv-qc                   # activate an env by name (or: hpv-activate qc)"
echo "  # SLURM / non-login batch jobs: add 'source /etc/profile.d/hpv-bioinfo.sh' at the top of the job script."
[ "$DRYRUN" -eq 0 ] && [ "${FAIL:-0}" -ne 0 ] && exit 1 || exit 0
