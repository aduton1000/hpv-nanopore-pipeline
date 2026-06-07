# /etc/profile.d/hpv-bioinfo.sh — system-wide auto-load for the shared HPV Nanopore
# analysis environment. Installed root:root 0644 by setup_hpv_nanopore_env.sh so every
# login shell sources it automatically (no manual `source` needed for interactive users).
#
# BEHAVIOUR AT LOAD: VARIABLES + PATH ONLY, and completely SILENT.
#   - exports HPV_PIPELINE, SCHEME_DIR, REF_DIR, MODELS_DIR, HPV_ARTIC_RUN (+ user-space
#     cache redirects) and prepends the shared bin/ (Nextflow + the hpv-run launcher) to PATH.
#   - prints NOTHING (no banner, no warnings, no errors) in any shell — interactive or not —
#     so it is safe for scp / sftp / rsync / non-interactive ssh.
#   - does NOT run conda init, does NOT `conda activate`, does NOT activate any env, and runs
#     no slow/heavy command. `conda activate hpv-<name>` / `hpv-activate` rely on the user's
#     own conda initialisation (exactly as for the lab's wgs_* envs); `hpv-run` needs no init
#     (it calls conda by absolute path).
#   - idempotent: re-sourcing (login + a manual source) never duplicates PATH or errors.

# --- resolved install root + data dirs (installer-written env file is VAR=val only; silent) ---
HPV_INSTALL_ROOT="/hpc/opt/hpv-bioinfo"
if [ -r "${HPV_INSTALL_ROOT}/etc/hpv-bioinfo.env" ]; then
  . "${HPV_INSTALL_ROOT}/etc/hpv-bioinfo.env"
fi
: "${SCHEME_DIR:=${HPV_INSTALL_ROOT}/primer-schemes}"
: "${REF_DIR:=${HPV_INSTALL_ROOT}/references}"
: "${MODELS_DIR:=${HPV_INSTALL_ROOT}/clair3_models}"
: "${HPV_PIPELINE:=${HPV_INSTALL_ROOT}/pipeline}"
export HPV_INSTALL_ROOT SCHEME_DIR REF_DIR MODELS_DIR HPV_PIPELINE

# --- shared bin (Nextflow + hpv-run launcher) on PATH, de-duplicated (idempotent) ---
case ":${PATH}:" in
  *":${HPV_INSTALL_ROOT}/bin:"*) ;;
  *) [ -d "${HPV_INSTALL_ROOT}/bin" ] && PATH="${HPV_INSTALL_ROOT}/bin:${PATH}" && export PATH ;;
esac

# --- user-space cache/work redirects (pure vars; uphold the HARD no-write-to-shared-pool rule) ---
export NXF_HOME="${NXF_HOME:-$HOME/.nextflow}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$HOME/.apptainer/cache}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$HOME/.conda/pkgs}"
export CLAIR3_MODEL_PATH="${MODELS_DIR}"

# --- default so `run_pipeline.py --use-artic` reaches hpv-artic (var only; user-overridable) ---
if [ -z "${HPV_ARTIC_RUN:-}" ]; then
  export HPV_ARTIC_RUN="${HPV_CONDA_HOME:-/hpc/opt/conda}/bin/conda run --no-capture-output -n hpv-artic"
fi

# --- convenience functions: DEFINED only (nothing is executed or printed at load) ---
# `hpv-activate <suffix>` == `conda activate hpv-<suffix>` (needs the user's normal conda init).
hpv-activate() {
  if [ -z "${1:-}" ]; then echo "usage: hpv-activate <qc|artic|asm|py|typing|smk>" >&2; return 2; fi
  conda activate "hpv-$1"
}
# Fallback `hpv-run` function — only defined if the launcher in $HPV_INSTALL_ROOT/bin is not on
# PATH; runs run_pipeline.py inside hpv-py via conda (absolute path; no init needed).
if ! command -v hpv-run >/dev/null 2>&1; then
  hpv-run() {
    "${HPV_CONDA_HOME:-/hpc/opt/conda}/bin/conda" run --no-capture-output -n hpv-py \
      python "${HPV_PIPELINE}/run_pipeline.py" "$@"
  }
fi

# Always end with success so a non-zero status here can never break a login shell.
:
