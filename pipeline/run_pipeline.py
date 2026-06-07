#!/usr/bin/env python
"""
run_pipeline.py — end-to-end driver: ingest -> read QC -> filter -> binning ->
per-type consensus (with rotated-reference junction fix) -> per-amplicon coverage QC ->
summary + publication figures. All outputs to the user's --outdir (default ./results);
nothing writes to the shared install. Mirrors the Snakefile/Nextflow stage graph.

Tools: minimap2/samtools/bcftools (+ pysam); artic/Clair3 used by the consensus stage
when present (--use-artic), otherwise the portable engine. chopper used if present.
"""
import argparse, csv, glob, gzip, os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

def sh(c): subprocess.run(c, shell=True, check=True, executable="/bin/bash")
def hpv(*a): subprocess.run([PY, str(HERE/"hpvpipe.py"), *map(str, a)], check=True)
def have(t):
    from shutil import which; return which(t) is not None

def read_sheet(path):
    rows = {}
    for r in csv.DictReader(open(path)):
        if not r.get("barcode") or r["barcode"].startswith("#"): continue
        rows[r["barcode"]] = r["sample_id"]
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-pass", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--scheme-dir", required=True)
    ap.add_argument("--ref", required=True, help="references/hpv_all_types.fasta")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--models-dir", default=os.environ.get("MODELS_DIR", ""))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-len", type=int, default=900)     # keeps HPV52 ~1.1 kb amp
    ap.add_argument("--max-len", type=int, default=4800)
    ap.add_argument("--min-reads", type=int, default=20)
    ap.add_argument("--use-artic", action="store_true", help="artic minion + Clair3 2.x (strict)")
    ap.add_argument("--require-artic", action="store_true", help="abort if artic fails (no fallback)")
    ap.add_argument("--allow-fallback", action="store_true", help="opt in to portable engine on artic failure")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args()
    R = Path(a.outdir).resolve()        # absolute: qc_report runs with cwd=server-setup/
    a.scheme_dir = str(Path(a.scheme_dir).resolve()); a.ref = str(Path(a.ref).resolve())
    a.fastq_pass = str(Path(a.fastq_pass).resolve())
    for d in ("reads","qc/lengths","filtered","typing","consensus","coverage","figures"):
        (R/d).mkdir(parents=True, exist_ok=True)
    sheet = read_sheet(a.samplesheet)

    # ---------- Stage 1: ingest ----------
    print("== Stage 1: ingest ==")
    for bc, sample in sheet.items():
        files = sorted(glob.glob(f"{a.fastq_pass}/{bc}/*.fastq.gz") + glob.glob(f"{a.fastq_pass}/{bc}/*.fastq"))
        outfq = R/"reads"/f"{sample}.fastq.gz"
        with gzip.open(outfq, "wt") as out:
            for f in files:
                op = gzip.open(f,"rt") if f.endswith(".gz") else open(f)
                out.write(op.read()); op.close()
        print(f"  {bc} -> {sample}: {len(files)} file(s)")

    # ---------- Stage 2: read QC ----------
    print("== Stage 2: read QC ==")
    rs = open(R/"qc"/"read_stats.tsv","w"); rs.write("barcode\tsample_id\tn_reads\ttotal_bases\tmean_len\tmean_qual\n")
    for bc, sample in sheet.items():
        lens=[]
        with open(R/"qc"/"lengths"/f"{sample}.lengths.tsv","w") as lh:
            op=gzip.open(R/"reads"/f"{sample}.fastq.gz","rt")
            while True:
                h=op.readline()
                if not h: break
                s=op.readline().strip(); op.readline(); op.readline()
                lens.append(len(s)); lh.write(f"{len(s)}\t15.0\n")
            op.close()
        n=len(lens); tb=sum(lens); ml=round(tb/n,1) if n else 0
        rs.write(f"{bc}\t{sample}\t{n}\t{tb}\t{ml}\t15.0\n")
    rs.close()

    # ---------- Stage 3: length/quality filter ----------
    print("== Stage 3: amplicon-aware length filter ==")
    for sample in sheet.values():
        inp=R/"reads"/f"{sample}.fastq.gz"; out=R/"filtered"/f"{sample}.filt.fastq.gz"
        if have("chopper"):
            sh(f"zcat -f '{inp}' | chopper --minlength {a.min_len} --maxlength {a.max_len} -q 9 2>/dev/null | gzip > '{out}'")
        else:
            with gzip.open(inp,"rt") as i, gzip.open(out,"wt") as o:
                while True:
                    h=i.readline()
                    if not h: break
                    s=i.readline(); p=i.readline(); q=i.readline()
                    if a.min_len <= len(s.strip()) <= a.max_len: o.write(h+s+p+q)
        print(f"  filtered {sample}")

    # ---------- Stage 4: type assignment (binning) ----------
    print("== Stage 4: HPV type assignment (binning) ==")
    tc = open(R/"typing"/"type_counts.tsv","w"); tc.write("sample_id\ttype\tn_reads\n")
    present = {}
    for sample in sheet.values():
        hpv("bin","--sample",sample,"--reads",R/"filtered"/f"{sample}.filt.fastq.gz",
            "--ref",a.ref,"--outdir",R/"typing","--min-reads",a.min_reads,"--threads",a.threads)
        pres=[]
        for ln in open(R/"typing"/f"{sample}.types.tsv").read().splitlines()[1:]:
            c=ln.split("\t"); tc.write(f"{c[0]}\t{c[1]}\t{c[2]}\n")
            if c[4]=="YES": pres.append(c[1])
        present[sample]=pres
    tc.close()

    # ---------- Stage 5+1: per-type consensus with junction fix ----------
    print("== Stage 5: per-type consensus (rotated-reference junction fix) ==")
    cq = open(R/"consensus"/"consensus_qc.tsv","w"); cq.write("sample_id\ttype\tpct_genome\tpct_N\tengine\tmodel\tseam_N\tjunction_engine\tjunction_bases_filled\n")
    engine_rows = []
    for sample, types in present.items():
        for t in types:
            reads=R/"typing"/f"{sample}.{t}.fastq.gz"
            if not reads.exists(): continue
            args=["consensus","--sample",sample,"--type",t,"--reads",reads,
                  "--scheme-dir",a.scheme_dir,"--outdir",R/"consensus","--threads",a.threads,
                  "--models-dir",a.models_dir]
            if a.use_artic: args.append("--use-artic")
            if a.require_artic: args.append("--require-artic")
            if a.allow_fallback: args.append("--allow-fallback")
            try:
                hpv(*args)
            except subprocess.CalledProcessError as e:
                print(f"[run_pipeline] ABORT: consensus failed for {sample}.{t} (artic strict mode; "
                      f"no fallback). Fix the artic/Clair3 path or pass --allow-fallback.", file=sys.stderr)
                sys.exit(e.returncode or 4)
            tag=f"{sample}.{t}"
            # collect coverage for qc_report
            for suf in ("amplicon.regions.bed.gz","per-base.bed.gz"):
                src=R/"consensus"/f"{tag}.{suf}"
                if src.exists(): sh(f"cp '{src}' '{R}/coverage/{tag}.{suf}'")
            qcf=R/"consensus"/f"{tag}.consensus_qc.tsv"
            if qcf.exists():
                for ln in qcf.read_text().splitlines()[1:]:
                    cq.write(ln+"\n"); engine_rows.append(ln.split("\t"))
    cq.close()

    # ---------- Stage 7/8: coverage QC + publication figures ----------
    if not a.no_figures:
        print("== Stage 7/8: publication QC figures ==")
        subprocess.run([PY, str(HERE.parent/"qc_report.py"), "--results", str(R),
                        "--scheme-dir", a.scheme_dir, "--outdir", str(R/"figures"), "--no-multiqc"],
                       check=False, cwd=str(HERE.parent))

    # ---------- run summary ----------
    with open(R/"run_summary.tsv","w") as fh:
        fh.write("sample_id\ttypes_present\n")
        for s,ts in present.items(): fh.write(f"{s}\t{';'.join(ts)}\n")
    print(f"\n== DONE == results in {R}/  (consensus, figures, run_summary.tsv)")
    print("  types present per sample:", {s:ts for s,ts in present.items()})
    # per (sample x type): engine actually used, %genome, %N, seam N, model
    print(f"  {'sample':14}{'type':7}{'engine':14}{'%genome':9}{'%N':7}{'seamN':7}model")
    for r in engine_rows:
        if len(r) >= 8:
            print(f"  {r[0]:14}{r[1]:7}{r[4]:14}{r[2]:9}{r[3]:7}{r[6]:7}{r[5]}")
    engines = sorted({r[4] for r in engine_rows if len(r) >= 5})
    print(f"  engine(s) that produced consensus: {engines}"
          + ("  (artic+Clair3 path exercised)" if "artic+clair3" in engines else
             "  (portable engine — use --use-artic --require-artic to exercise Clair3)"))

if __name__ == "__main__":
    main()
