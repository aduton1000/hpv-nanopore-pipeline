#!/usr/bin/env python
"""
qc_report.py — publication-quality, amplicon-aware QC + visualisation for HPV
Nanopore tiling-amplicon runs. Run in the `hpv-py` env. All figures are written
(SVG+PDF+300dpi PNG via qc_lib) to a figures/ dir in the USER's working area
(default ./results/figures), never to /opt.

Inputs are the conventional outputs of PIPELINE_STAGES.md (all optional — each figure
is skipped with a message if its input is absent):
  results/qc/read_stats.tsv            barcode, sample_id, n_reads, total_bases, mean_len, mean_qual
  results/qc/lengths/<sample>.lengths.tsv      one read length (and optional qual) per line
  results/coverage/<sample>.<type>.amplicon.regions.bed.gz   mosdepth -b amplicon.bed output
  results/coverage/<sample>.<type>.per-base.bed.gz           mosdepth per-base output
  results/typing/type_counts.tsv       sample_id, type, n_reads
  results/consensus/consensus_qc.tsv   sample_id, type, pct_genome, pct_N
  (optional) a sequencing_summary*.txt for pycoQC run metrics

Usage:
  python qc_report.py --results results [--outdir results/figures] \
                      --scheme-dir "$SCHEME_DIR" [--sequencing-summary seq_summary.txt]
  python qc_report.py --demo --outdir /tmp/hpv_qc_demo   # synthesize tiny data, prove figures
"""
from __future__ import annotations
import argparse, glob, gzip, os, subprocess, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import qc_lib as Q
Q.set_theme()

FOUND, SKIPPED = [], []
def did(name, paths): FOUND.append(name); print(f"  [fig] {name} -> {', '.join(Path(p).name for p in paths)}")
def skip(name, why): SKIPPED.append(name); print(f"  [skip] {name}: {why}")


# ----------------------------- loaders -----------------------------
def _read_tsv(p):
    return pd.read_csv(p, sep="\t") if Path(p).exists() else None

def load_amplicon_depth(results):
    """Build a long df: sample, type, amplicon, mean_depth, pool from mosdepth region beds."""
    rows = []
    for f in glob.glob(f"{results}/coverage/*.amplicon.regions.bed.gz"):
        base = Path(f).name.replace(".amplicon.regions.bed.gz", "")
        if "." not in base: continue
        sample, htype = base.split(".", 1)
        op = gzip.open(f, "rt")
        for ln in op:
            c = ln.rstrip("\n").split("\t")
            # mosdepth -b: chrom start end name mean
            if len(c) < 5: continue
            name, mean = c[3], float(c[4])
            amp = name.split("_")[1] if "_" in name else name
            rows.append(dict(sample=sample, type=htype, amplicon=amp, mean_depth=mean))
    return pd.DataFrame(rows) if rows else None


# ----------------------------- figures -----------------------------
def fig_read_counts(read_stats, outdir):
    if read_stats is None: return skip("read_counts_yield", "results/qc/read_stats.tsv absent")
    df = read_stats.sort_values("barcode")
    fig, axes = plt.subplots(1, 2, figsize=(Q.FULL_W, 2.6))
    axes[0].bar(df["barcode"], df["n_reads"], color=Q.OKABE_ITO[0])
    axes[0].set(title="Reads per barcode", xlabel="barcode", ylabel="read count")
    axes[1].bar(df["barcode"], df["total_bases"]/1e6, color=Q.OKABE_ITO[1])
    axes[1].set(title="Yield per barcode", xlabel="barcode", ylabel="yield (Mb)")
    for ax in axes: ax.tick_params(axis="x", rotation=90)
    Q.footer(fig, "HPV Nanopore QC — read_stats.tsv")
    did("read_counts_yield", Q.save_fig(fig, outdir, "01_read_counts_yield"))

def fig_length_quality(results, outdir):
    files = glob.glob(f"{results}/qc/lengths/*.lengths.tsv")
    if not files: return skip("length_quality_dist", "results/qc/lengths/*.lengths.tsv absent")
    allL, allQ = [], []
    fig, axes = plt.subplots(1, 3, figsize=(Q.FULL_W, 2.6))
    for f in files:
        d = pd.read_csv(f, sep="\t", header=None)
        L = d[0].values; allL.append(L)
        axes[0].hist(L, bins=80, histtype="step", lw=0.8, alpha=0.7)
        if d.shape[1] > 1: allQ.append(d[1].values)
    allLc = np.concatenate(allL)
    axes[0].set(title="Read-length distribution (per barcode)", xlabel="read length (bp)", ylabel="reads")
    Q.annotate_length_window(axes[0]); axes[0].legend(fontsize=5.5, loc="upper right")
    axes[0].set_xlim(0, min(6000, allLc.max()))
    if allQ:
        allQc = np.concatenate(allQ)
        axes[1].hist(allQc, bins=60, color=Q.OKABE_ITO[2]); axes[1].set(title="Quality distribution (combined)", xlabel="mean read Q", ylabel="reads")
        h = axes[2].hist2d(allLc[:len(allQc)], allQc, bins=[80,50], cmap="magma")
        axes[2].set(title="Length vs quality", xlabel="read length (bp)", ylabel="mean read Q")
        axes[2].axvspan(Q.AMPLICON_MIN_BP, Q.AMPLICON_MAX_BP, color="w", alpha=0.0)
        axes[2].axvline(Q.HPV52_GAPFILL_BP, color="#00E0FF", ls="--", lw=0.8)
        axes[2].set_xlim(0, min(6000, allLc.max()))
    else:
        for a in (axes[1], axes[2]): a.set_visible(False)
    Q.footer(fig, "length window: 2–4 kb shaded; HPV52 1.1 kb gap-fill marked")
    did("length_quality_dist", Q.save_fig(fig, outdir, "02_length_quality"))

def fig_amplicon_heatmap(dep, outdir):
    if dep is None: return skip("per_amplicon_coverage_heatmap", "no mosdepth amplicon beds")
    for htype, sub in dep.groupby("type"):
        mat = sub.pivot_table(index="sample", columns="amplicon", values="mean_depth", aggfunc="mean")
        mat = mat.reindex(sorted(mat.columns, key=lambda x: (len(x), x)), axis=1)
        fig, ax = plt.subplots(figsize=(max(Q.COL_W, 0.5*mat.shape[1]+1.5),
                                        max(2.0, 0.3*mat.shape[0]+1.2)))
        im = ax.imshow(np.log10(mat.values+1), aspect="auto", cmap=Q.DEPTH_CMAP)
        ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns)
        ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index)
        ax.set(title=f"{htype}: per-amplicon mean depth (log10) — dropout at a glance",
               xlabel="amplicon", ylabel="sample")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("log10(depth+1)")
        # mark zero/low-depth (dropout) cells
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat.values[i, j] < 20:
                    ax.text(j, i, "·", ha="center", va="center", color="red", fontsize=8)
        Q.footer(fig, "red dot = mean depth <20x (amplicon dropout)")
        did(f"amplicon_heatmap_{htype}", Q.save_fig(fig, outdir, f"03_amplicon_heatmap_{htype}"))

def fig_pool_balance(dep, scheme_dir, outdir):
    if dep is None: return skip("pool_balance", "no amplicon depth")
    pool_of = {}
    for f in glob.glob(f"{scheme_dir}/*/V5/*.primer.bed") if scheme_dir else []:
        for ln in open(f):
            c = ln.split("\t")
            if len(c) >= 5 and "_" in c[3]:
                pool_of[(c[0], c[3].split("_")[1])] = c[4]
    if not pool_of: return skip("pool_balance", "scheme primer.bed pools not found (need --scheme-dir)")
    dep = dep.copy()
    dep["pool"] = [pool_of.get((f"{r.type}", r.amplicon)) for r in dep.itertuples()]
    g = dep.dropna(subset=["pool"]).groupby(["sample","pool"])["mean_depth"].mean().unstack("pool")
    if g is None or g.empty or g.shape[1] < 2: return skip("pool_balance", "need both pools")
    fig, ax = plt.subplots(figsize=(Q.COL_W+0.5, Q.COL_W))
    ax.scatter(g.iloc[:,0], g.iloc[:,1], color=Q.OKABE_ITO[0])
    m = float(np.nanmax(g.values))*1.05 if np.isfinite(np.nanmax(g.values)) else 1
    ax.plot([0,m],[0,m],"--",color="#999999",lw=0.8)
    for s,row in g.iterrows(): ax.annotate(s, (row.iloc[0],row.iloc[1]), fontsize=5)
    ax.set(title="Pool balance (mean amplicon depth)", xlabel=f"Pool {g.columns[0]} mean depth",
           ylabel=f"Pool {g.columns[1]} mean depth", xlim=(0,m), ylim=(0,m))
    Q.footer(fig, "off-diagonal points flag pool imbalance/failure")
    did("pool_balance", Q.save_fig(fig, outdir, "05_pool_balance"))

def fig_typing(type_counts, outdir):
    if type_counts is None: return skip("type_assignment", "results/typing/type_counts.tsv absent")
    mat = type_counts.pivot_table(index="sample_id", columns="type", values="n_reads",
                                  aggfunc="sum", fill_value=0)
    mat = mat.reindex([t for t in Q.HPV_TYPES if t in mat.columns], axis=1)
    fig, ax = plt.subplots(figsize=(Q.FULL_W, max(2.2, 0.3*mat.shape[0]+1.2)))
    bottom = np.zeros(mat.shape[0])
    for t in mat.columns:
        ax.barh(mat.index, mat[t], left=bottom, color=Q.TYPE_COLORS.get(t,"#888"), label=t)
        bottom += mat[t].values
    ax.set(title="Reads assigned per HPV type per sample (co-infection view)",
           xlabel="reads", ylabel="sample")
    ax.legend(ncol=6, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    Q.footer(fig, "from competitive minimap2 / kraken2 binning")
    did("type_assignment", Q.save_fig(fig, outdir, "06_type_assignment"))

def fig_consensus_qc(cons, outdir):
    if cons is None: return skip("consensus_completeness", "results/consensus/consensus_qc.tsv absent")
    mat = cons.pivot_table(index="sample_id", columns="type", values="pct_genome", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(max(Q.COL_W,0.5*mat.shape[1]+1.5), max(2.0,0.3*mat.shape[0]+1.2)))
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index)
    ax.set(title="Consensus completeness: % genome recovered", xlabel="type", ylabel="sample")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("% genome")
    Q.footer(fig, "pair with %N table; <90% = incomplete consensus")
    did("consensus_completeness", Q.save_fig(fig, outdir, "07_consensus_completeness"))

def fig_depth_tracks(results, scheme_dir, outdir):
    files = glob.glob(f"{results}/coverage/*.per-base.bed.gz")
    if not files: return skip("depth_tracks", "results/coverage/*.per-base.bed.gz absent")
    # one example track per (sample,type) found (full set in pipeline run)
    for f in files[:12]:
        base = Path(f).name.replace(".per-base.bed.gz",""); sample,htype = base.split(".",1)
        pos, dep = [], []
        for ln in gzip.open(f,"rt"):
            c=ln.split("\t"); s,e,d=int(c[1]),int(c[2]),float(c[3])
            pos.append((s+e)//2); dep.append(d)
        fig, ax = plt.subplots(figsize=(Q.FULL_W, 2.0))
        ax.fill_between(pos, dep, color=Q.OKABE_ITO[0], alpha=0.5, lw=0)
        # overlay amplicon boundaries + pool colours
        bed=f"{scheme_dir}/{htype}/V5/{htype}.amplicon.bed" if scheme_dir else None
        if bed and Path(bed).exists():
            for ln in open(bed):
                c=ln.split("\t"); a0,a1=int(c[1]),int(c[2]); pool=c[4]
                ax.axvspan(a0,a1,color=Q.POOL_COLORS.get(pool,"#ccc"),alpha=0.06)
                if "_part" in c[3]: ax.axvline(a0,color="red",ls=":",lw=0.6)  # origin-junction edge
        ax.set(title=f"{sample} × {htype}: depth across genome (pool-coloured; red = origin-junction)",
               xlabel="canonical position (bp)", ylabel="depth")
        ax.set_yscale("symlog")
        Q.footer(fig, "amplicon spans shaded by pool; origin-junction amplicon split shown")
        did(f"depth_track_{sample}_{htype}", Q.save_fig(fig, outdir, f"04_depth_{sample}_{htype}"))

def run_pycoqc(seqsum, outdir):
    if not seqsum:
        return skip("pycoqc", "no sequencing_summary uploaded (ASK user to upload it from the GridION "
                              "— FASTQ alone cannot produce channel/yield-over-time metrics)")
    out = Path(outdir).parent / "pycoQC.html"
    try:
        subprocess.run(["pycoQC","-f",seqsum,"-o",str(out)], check=True)
        did("pycoqc", [str(out)])
    except Exception as e:
        skip("pycoqc", f"pycoQC failed ({e})")

def run_multiqc(results, outdir):
    try:
        subprocess.run(["multiqc","-f","-o",str(Path(outdir).parent),results], check=True)
        did("multiqc", [str(Path(outdir).parent/'multiqc_report.html')])
    except Exception as e:
        skip("multiqc", f"multiqc not run ({e})")


# ----------------------------- demo data -----------------------------
def make_demo(results):
    rng = np.random.default_rng(0)
    R = Path(results); (R/"qc/lengths").mkdir(parents=True, exist_ok=True)
    (R/"coverage").mkdir(parents=True, exist_ok=True); (R/"typing").mkdir(parents=True, exist_ok=True)
    (R/"consensus").mkdir(parents=True, exist_ok=True)
    samples = ["SAMPLE_001","SAMPLE_002","SAMPLE_003"]; barcodes=["barcode01","barcode02","barcode03"]
    rs=[]
    for s,b in zip(samples,barcodes):
        n=int(rng.integers(4000,15000)); L=rng.normal(2600,500,n).clip(800,4800)
        rs.append(dict(barcode=b,sample_id=s,n_reads=n,total_bases=int(L.sum()),
                       mean_len=float(L.mean()),mean_qual=float(rng.normal(14,1))))
        q=rng.normal(14,2,n).clip(5,25)
        pd.DataFrame({0:L.astype(int),1:q.round(1)}).to_csv(R/f"qc/lengths/{s}.lengths.tsv",sep="\t",header=False,index=False)
    pd.DataFrame(rs).to_csv(R/"qc/read_stats.tsv",sep="\t",index=False)
    # amplicon depth + per-base for HPV16/HPV58, with one dropout
    tc=[]; cq=[]
    for s in samples:
        for ht,glen in (("HPV16",7921),("HPV58",7883)):
            rows=[]
            for amp in range(1,7):
                d = 0 if (s=="SAMPLE_003" and ht=="HPV16" and amp==4) else int(rng.integers(40,400))
                rows.append(f"{ht}\t{amp*100}\t{amp*100+1500}\t{ht}_{amp}\t{d}\n")
            with gzip.open(R/f"coverage/{s}.{ht}.amplicon.regions.bed.gz","wt") as fh: fh.writelines(rows)
            with gzip.open(R/f"coverage/{s}.{ht}.per-base.bed.gz","wt") as fh:
                for p in range(0,glen,200):
                    fh.write(f"{ht}\t{p}\t{p+200}\t{int(rng.integers(20,300))}\n")
            tc.append(dict(sample_id=s,type=ht,n_reads=int(rng.integers(1000,8000))))
            cq.append(dict(sample_id=s,type=ht,pct_genome=float(rng.uniform(85,100)),pct_N=float(rng.uniform(0,8))))
    pd.DataFrame(tc).to_csv(R/"typing/type_counts.tsv",sep="\t",index=False)
    pd.DataFrame(cq).to_csv(R/"consensus/consensus_qc.tsv",sep="\t",index=False)
    print(f"  demo data written under {R}")


def main():
    ap = argparse.ArgumentParser(description="HPV Nanopore publication-quality QC report")
    ap.add_argument("--results", default="results", help="pipeline results dir (read inputs)")
    ap.add_argument("--outdir", default=None, help="figures dir (default <results>/figures, in CWD)")
    ap.add_argument("--scheme-dir", default=os.environ.get("SCHEME_DIR"))
    ap.add_argument("--sequencing-summary", default=None)
    ap.add_argument("--demo", action="store_true", help="synthesize tiny demo data then report")
    ap.add_argument("--no-multiqc", action="store_true")
    a = ap.parse_args()
    results = Path(a.results); outdir = Path(a.outdir) if a.outdir else results/"figures"
    if a.demo:
        results.mkdir(parents=True, exist_ok=True); make_demo(results)
    print(f"== HPV QC report ==  results={results}  figures->{outdir}")
    print("  (figures: SVG+PDF+300dpi PNG, colourblind-safe, embedded fonts)")

    read_stats = _read_tsv(results/"qc"/"read_stats.tsv")
    dep        = load_amplicon_depth(results)
    type_counts= _read_tsv(results/"typing"/"type_counts.tsv")
    cons       = _read_tsv(results/"consensus"/"consensus_qc.tsv")

    fig_read_counts(read_stats, outdir)
    fig_length_quality(results, outdir)
    fig_amplicon_heatmap(dep, outdir)          # headline amplicon-QC figure
    fig_depth_tracks(results, a.scheme_dir, outdir)
    fig_pool_balance(dep, a.scheme_dir, outdir)
    fig_typing(type_counts, outdir)
    fig_consensus_qc(cons, outdir)
    run_pycoqc(a.sequencing_summary, outdir)
    if not a.no_multiqc: run_multiqc(str(results), outdir)

    print(f"\n== done ==  {len(FOUND)} figure groups written, {len(SKIPPED)} skipped.")
    print(f"   figures in: {outdir}  (open the MultiQC report for the combined view)")
    if not FOUND:
        print("   no inputs found — run with --demo to prove the module, or point --results at a real run.")
        sys.exit(1)


if __name__ == "__main__":
    main()
