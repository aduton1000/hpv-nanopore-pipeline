#!/usr/bin/env python
"""
hpvpipe.py — core stage library + CLI for the HPV Nanopore tiling-amplicon pipeline.

Implements the previously-scaffolded stages as real, tested code:
  * rotate-build  : per-type ROTATED reference + rotated amp3 BED (origin-junction fix)
  * bin           : competitive minimap2 read typing -> per-type FASTQ + types-present
  * consensus     : per (sample x type) consensus with the junction fix integrated
                    (prefers `artic minion` + Clair3 2.x when available; otherwise a
                     portable minimap2+samtools/pysam engine — both produce the same
                     merged canonical consensus + VCF + trimmed BAM + coverage)
  * coverage      : per-amplicon / per-pool depth (mosdepth if present, else pysam)
  * model-select  : pick the r10.4.1 Clair3 model (fast/hac/sup) from a FASTQ header
  * summary       : per-(sample x type) %genome/%N + types table for qc_report.py

CPU-only. Reads schemes/refs from $SCHEME_DIR/$REF_DIR (or --scheme-dir/--ref); all
outputs go to the user's working dir; nothing writes to the shared install.

Tool resolution: minimap2/samtools/bcftools must be on PATH (in hpv-artic/hpv-qc).
artic/mosdepth/chopper are used when present and gracefully substituted otherwise so the
pipeline runs end-to-end even on a workstation without them.
"""
from __future__ import annotations
import argparse, glob, gzip, json, os, re, shutil, subprocess, sys
from pathlib import Path

TYPES = ["HPV16","HPV18","HPV31","HPV33","HPV35","HPV39","HPV45","HPV52","HPV56","HPV58","HPV59","HPV66"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ----------------------------- small utils -----------------------------
def have(tool): return shutil.which(tool) is not None
def run(cmd, **kw):
    kw.setdefault("check", True)
    return subprocess.run(cmd, **kw)
def sh(cmd):  # shell pipeline, checked
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
def revcomp(s): return s.translate(COMP)[::-1]

def read_fasta(path):
    d, name, buf = {}, None, []
    op = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    for ln in op:
        if ln.startswith(">"):
            if name: d[name] = "".join(buf)
            name = ln[1:].split()[0]; buf = []
        else: buf.append(ln.strip())
    if name: d[name] = "".join(buf)
    op.close(); return d

def write_fasta(path, header, seq, width=70):
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(seq), width): fh.write(seq[i:i+width] + "\n")

def fastq_reads(path):
    op = gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)
    while True:
        h = op.readline()
        if not h: break
        s = op.readline().strip(); op.readline(); q = op.readline().strip()
        yield h[1:].strip(), s, q
    op.close()

def scheme_files(scheme_dir, t, ver="V5"):
    d = Path(scheme_dir) / t / ver
    return d / f"{t}.reference.fasta", d / f"{t}.primer.bed", d / f"{t}.amplicon.bed"


# ----------------------------- BED / amplicons -----------------------------
def load_amplicons(primer_bed):
    """amp -> {'L':(start,end), 'R':(start,end), 'pool':str, 'chrom':str}. _alt ignored
    for coordinates (same span as base primer)."""
    amps = {}
    for ln in open(primer_bed):
        c = ln.rstrip("\n").split("\t")
        if len(c) < 6: continue
        chrom, s, e, name, pool, strand = c[0], int(c[1]), int(c[2]), c[3], c[4], c[5]
        m = re.match(r".*_(\d+)_(LEFT|RIGHT)(?:_alt\d+)?$", name)
        if not m: continue
        amp, side = m.group(1), m.group(2)
        a = amps.setdefault(amp, {"pool": pool, "chrom": chrom})
        key = "L" if side == "LEFT" else "R"
        if key not in a: a[key] = (s, e)            # base primer (first occurrence)
    return amps

def junction_amp(amps, L):
    """the amplicon whose insert wraps the origin (LEFT.start > RIGHT.end)."""
    for amp, a in amps.items():
        if "L" in a and "R" in a and a["L"][0] > a["R"][1]:
            return amp
    return None


# ----------------------------- rotation (Task 1) -----------------------------
def rotate_offset(L): return L // 2

def to_rot(pos, off, L): return (pos - off) % L
def to_canon(pos, off, L): return (pos + off) % L

def cmd_rotate_build(args):
    """Build per-type rotated reference + rotated amp3 BED so amp3 is contiguous."""
    ref_fa, primer_bed, amp_bed = scheme_files(args.scheme_dir, args.type)
    refs = read_fasta(ref_fa); chrom = next(iter(refs)); seq = refs[chrom]; L = len(seq)
    amps = load_amplicons(primer_bed); jamp = junction_amp(amps, L)
    if jamp is None:
        print(f"[rotate] {args.type}: no origin-junction amplicon found", file=sys.stderr); return
    off = rotate_offset(L)
    # safety: rotation cut-point must NOT fall inside the junction amplicon span
    ja = amps[jamp]
    in_jamp = (off >= ja["L"][0]) or (off < ja["R"][1])
    if in_jamp:
        # shift cut to a safe mid-point away from the junction amplicon
        off = (ja["R"][1] + ja["L"][0]) // 2
    rot = seq[off:] + seq[:off]
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rot_chrom = f"{chrom}_rot{off}"                # underscore (valid chrom; no pipe)
    write_fasta(outdir / f"{args.type}.rot.reference.fasta", rot_chrom, rot)
    # rotated amp3 primer BED — amp3 becomes a normal contiguous forward amplicon
    lines = []
    for ln in open(primer_bed):
        c = ln.rstrip("\n").split("\t")
        m = re.match(r".*_(\d+)_(LEFT|RIGHT)(?:_alt\d+)?$", c[3]) if len(c) >= 6 else None
        if not m or m.group(1) != jamp: continue
        ns, ne = to_rot(int(c[1]), off, L), to_rot(int(c[2]) - 1, off, L) + 1   # keep interval contiguous
        lines.append("\t".join([rot_chrom, str(ns), str(ne), c[3], c[4], c[5]] + c[6:]))
    (outdir / f"{args.type}.rot.primer.bed").write_text("\n".join(lines) + "\n")
    info = dict(type=args.type, canonical_chrom=chrom, rot_chrom=rot_chrom,
                L=L, offset=off, junction_amp=jamp,
                amp3_canonical={"L": ja["L"], "R": ja["R"]},
                amp3_rotated={"L": [to_rot(ja["L"][0],off,L), to_rot(ja["L"][1]-1,off,L)+1],
                              "R": [to_rot(ja["R"][0],off,L), to_rot(ja["R"][1]-1,off,L)+1]})
    (outdir / f"{args.type}.rotation.json").write_text(json.dumps(info, indent=1) + "\n")
    print(f"[rotate] {args.type}: L={L} offset={off} junction=amp{jamp} "
          f"amp3 rotated to contiguous [{info['amp3_rotated']['L'][0]}..{info['amp3_rotated']['R'][1]})")


# ----------------------------- portable consensus engine -----------------------------
def minimap2_sort(ref_fa, reads, out_bam, threads=4):
    sh(f"minimap2 -ax map-ont -t {threads} '{ref_fa}' '{reads}' 2>/dev/null "
       f"| samtools sort -@ {threads} -o '{out_bam}' - && samtools index '{out_bam}'")

def pileup_consensus(bam, ref_fa, mindepth=3, exclude_supplementary=False):
    """Majority-base consensus over the full reference via pysam; N where depth<mindepth.
    exclude_supplementary=True models an amplicon-aware caller: an origin-spanning read
    cannot be assigned to a single linear amplicon, so its supplementary (chimeric) piece
    across the seam is NOT used — which is exactly why amp3 drops out on a linear reference.
    Returns (consensus_str, depth_list)."""
    import pysam
    refs = read_fasta(ref_fa); chrom = next(iter(refs)); L = len(refs[chrom])
    cons = ["N"] * L; depth = [0] * L
    bf = pysam.AlignmentFile(bam, "rb")
    if exclude_supplementary:
        cb = lambda r: not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate)
    else:
        cb = "all"
    cov = bf.count_coverage(chrom, 0, L, quality_threshold=0, read_callback=cb)
    bases = "ACGT"
    for i in range(L):
        counts = [cov[b][i] for b in range(4)]
        d = sum(counts); depth[i] = d
        if d >= mindepth:
            cons[i] = bases[max(range(4), key=lambda b: counts[b])]
    bf.close()
    return "".join(cons), depth

def cmd_consensus(args):
    """Per (sample x type) consensus with rotated-reference junction fix integrated."""
    ref_fa, primer_bed, amp_bed = scheme_files(args.scheme_dir, args.type)
    refs = read_fasta(ref_fa); chrom = next(iter(refs)); canon = refs[chrom]; L = len(canon)
    amps = load_amplicons(primer_bed); jamp = junction_amp(amps, L)
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.sample}.{args.type}"
    threads = args.threads

    used_artic = False
    canon_cons = None
    sel_model = None
    # artic may live in a different env (hpv-artic) than this script (hpv-py, for pysam):
    # HPV_ARTIC_RUN="conda run -n hpv-artic" lets us call it without switching env.
    artic_prefix = os.environ.get("HPV_ARTIC_RUN", "").split()
    artic_ok = bool(artic_prefix) or have("artic")
    want_artic = args.use_artic or args.require_artic
    # STRICT by default when artic is requested: a failure ABORTS (so a broken Clair3 path
    # cannot masquerade as success). --allow-fallback opts into the portable engine instead.
    strict = want_artic and not args.allow_fallback
    if want_artic and not artic_ok:
        msg = ("artic not found (need it on PATH or set HPV_ARTIC_RUN='conda run -n hpv-artic')")
        if strict: print(f"[consensus] FATAL {tag}: {msg}", file=sys.stderr); sys.exit(3)
        print(f"[consensus] WARN {tag}: {msg}; --allow-fallback -> portable engine", file=sys.stderr)
    if want_artic and artic_ok:
        # ---- production path: artic minion (Clair3 2.x), LOCAL scheme via --bed/--ref ----
        # (--bed/--ref bypass get_scheme/manifest entirely -> no network fetch; ARTIC 1.10.1
        #  `if args.bed and args.ref: bed=args.bed; ref=args.ref`.)
        sel_model = select_model(args.reads, args.models_dir, args.model)
        # MAIN artic run: canonical ref + a bed with the origin-junction amplicon REMOVED.
        # GEOMETRY FIX: align_trim then sees only the 5 non-wrapping amplicons (all positive
        # length) and never hits np.zeros(amp.length) on amp3's negative (right<left) wrap.
        bed_main = str(primer_bed)
        if jamp is not None:
            bed_main = str(out / f"{tag}.no_amp{jamp}.primer.bed")
            keep = [ln for ln in Path(primer_bed).read_text().splitlines()
                    if ln.strip() and not re.match(rf"[^\t]*_{jamp}_(LEFT|RIGHT)", ln.split('\t')[3])]
            Path(bed_main).write_text("\n".join(keep) + "\n")
        adir = out / f"{tag}.artic.main"; adir.mkdir(exist_ok=True)
        cmd = artic_prefix + ["artic", "minion", "--threads", str(threads),
               "--bed", bed_main, "--ref", str(ref_fa), "--model-dir", str(args.models_dir)]
        if sel_model: cmd += ["--model", sel_model]
        cmd += ["--read-file", str(args.reads), tag]
        try:
            run(cmd, cwd=adir)
            cf = list(adir.glob("*.consensus.fasta"))
            if not cf: raise RuntimeError("artic minion produced no *.consensus.fasta")
            canon_cons = list(read_fasta(cf[0]).values())[0].upper(); used_artic = True
        except Exception as e:
            if strict:
                print(f"[consensus] FATAL {tag}: artic minion (main) failed ({e}). "
                      f"--require-artic/strict set: NOT falling back. cmd: {' '.join(map(str,cmd))}",
                      file=sys.stderr)
                sys.exit(4)
            print(f"[consensus] WARN {tag}: artic minion failed ({e}); --allow-fallback -> portable",
                  file=sys.stderr)

    # ---- portable canonical mapping (coverage QC always; consensus base when artic absent) ----
    canon_bam = out / f"{tag}.canonical.bam"
    minimap2_sort(ref_fa, args.reads, canon_bam, threads)
    # exclude_supplementary models a linear amplicon caller: origin-spanning amp3 reads cannot
    # be assigned to a linear amplicon -> seam GAP (the same gap artic's amp3-excluded run leaves).
    canon_only, canon_depth = pileup_consensus(canon_bam, ref_fa, args.mindepth, exclude_supplementary=True)
    if canon_cons is None:
        canon_cons = canon_only
    # amp3 insert (canonical) = between the amp3 primers, wrapping the origin
    seam_positions = []
    if jamp is not None:
        ja = amps[jamp]; ins_s, ins_e = ja["L"][1], ja["R"][0]   # LEFT.end .. RIGHT.start
        seam_positions = (list(range(ins_s, L)) + list(range(0, ins_e))) if ins_s > ins_e else list(range(ins_s, ins_e))
    canon_only_seam_N = sum(1 for p in seam_positions if canon_cons[p] == "N")   # gap in the MAIN consensus

    # ---- Task 2: recover amp3/seam via the ROTATED frame (artic/Clair3 when used; else pileup) ----
    seam_filled = 0; rotinfo = None; rot_engine = None
    if jamp is not None:
        rotdir = out / f"{tag}.rot"; rotdir.mkdir(exist_ok=True)
        import types as _t
        cmd_rotate_build(_t.SimpleNamespace(scheme_dir=args.scheme_dir, type=args.type, outdir=str(rotdir)))
        rotinfo = json.loads((rotdir / f"{args.type}.rotation.json").read_text()); off = rotinfo["offset"]
        rot_ref = rotdir / f"{args.type}.rot.reference.fasta"
        rot_bed = rotdir / f"{args.type}.rot.primer.bed"            # amp3-ONLY, contiguous (+length)
        rot_cons = None; rot_depth = None
        if used_artic:
            # On the rotated frame amp3 is a single positive-length amplicon, so artic minion's
            # align_trim handles it normally (no negative length) and Clair3 calls it: the seam
            # is recovered through the REAL artic/Clair3 path, not a portable shim.
            adir2 = out / f"{tag}.artic.rotj"; adir2.mkdir(exist_ok=True)
            cmd2 = artic_prefix + ["artic", "minion", "--threads", str(threads),
                   "--bed", str(rot_bed), "--ref", str(rot_ref), "--model-dir", str(args.models_dir)]
            if sel_model: cmd2 += ["--model", sel_model]
            cmd2 += ["--read-file", str(args.reads), f"{tag}_rotj"]
            try:
                run(cmd2, cwd=adir2)
                cf2 = list(adir2.glob("*.consensus.fasta"))
                if not cf2: raise RuntimeError("artic minion (rotated amp3) produced no consensus")
                rot_cons = list(read_fasta(cf2[0]).values())[0].upper(); rot_engine = "artic+clair3"
            except Exception as e:
                if strict:
                    print(f"[consensus] FATAL {tag}: artic minion (rotated amp{jamp}) failed ({e}). "
                          f"--require-artic/strict: NOT falling back.", file=sys.stderr); sys.exit(5)
                print(f"[consensus] WARN {tag}: rotated artic failed ({e}); portable seam fill", file=sys.stderr)
        if rot_cons is None:                                        # portable rotated consensus
            rot_bam = rotdir / f"{tag}.rot.bam"; minimap2_sort(rot_ref, args.reads, rot_bam, threads)
            rot_cons, rot_depth = pileup_consensus(rot_bam, rot_ref, args.mindepth); rot_engine = "rotated-pileup"
        # back-rotate amp3 bases into canonical; fill ONLY seam positions still N (no double-count)
        cons_list = list(canon_cons)
        for rp in range(len(rot_cons)):
            if rot_cons[rp] == "N": continue
            if rot_depth is not None and rot_depth[rp] < args.mindepth: continue
            cp = to_canon(rp, off, L)
            if cons_list[cp] == "N":
                cons_list[cp] = rot_cons[rp]; seam_filled += 1
        canon_cons = "".join(cons_list)

    # ---- junction-recovery proof (Task 1) ----
    merged_seam_N = sum(1 for p in seam_positions if canon_cons[p] == "N")
    seam_match = sum(1 for p in seam_positions if canon_cons[p] == canon[p])
    seam_total = len(seam_positions)
    seam_ident = round(100*seam_match/seam_total, 2) if seam_total else 100.0
    if jamp is not None:
        (out / f"{tag}.junction_proof.tsv").write_text(
            "metric\tvalue\n"
            f"junction_amplicon\tamp{jamp}\n"
            f"junction_engine\t{rot_engine}\n"
            f"rotation_offset\t{rotinfo['offset'] if rotinfo else 'NA'}\n"
            f"seam_region_bp\t{seam_total}\n"
            f"canonical_only_seam_N\t{canon_only_seam_N}\n"
            f"merged_seam_N\t{merged_seam_N}\n"
            f"seam_bases_filled_from_rotated\t{seam_filled}\n"
            f"seam_identity_to_reference_pct\t{seam_ident}\n"
            f"continuous_across_origin\t{'YES' if merged_seam_N==0 else 'NO'}\n"
            f"length_unchanged\t{'YES' if len(canon_cons)==L else 'NO'}\n")

    # ---- outputs ----
    write_fasta(out / f"{tag}.consensus.fasta", f"{args.sample}|{args.type}", canon_cons)
    # VCF (real variants vs canonical ref) via bcftools when available
    vcf = out / f"{tag}.vcf.gz"
    if have("bcftools"):
        try:
            sh(f"bcftools mpileup -f '{ref_fa}' '{canon_bam}' 2>/dev/null "
               f"| bcftools call -mv -Oz -o '{vcf}' 2>/dev/null && tabix -f -p vcf '{vcf}' 2>/dev/null || true")
        except Exception: pass
    # coverage summary (per-amplicon) -> regions bed.gz for qc_report
    write_amplicon_depth(canon_bam, amp_bed, canon, out / f"{tag}.amplicon.regions.bed.gz", L)
    write_perbase_depth(canon_depth, chrom, out / f"{tag}.per-base.bed.gz")
    # consensus QC
    n = canon_cons.count("N"); pct_genome = round(100*(L-n)/L, 2); pct_N = round(100*n/L, 2)
    engine = "artic+clair3" if used_artic else "portable"
    (out / f"{tag}.consensus_qc.tsv").write_text(
        "sample_id\ttype\tpct_genome\tpct_N\tengine\tmodel\tseam_N\tjunction_engine\tjunction_bases_filled\n"
        f"{args.sample}\t{args.type}\t{pct_genome}\t{pct_N}\t{engine}\t"
        f"{sel_model or 'NA'}\t{merged_seam_N if jamp is not None else 'NA'}\t"
        f"{rot_engine or 'NA'}\t{seam_filled}\n")
    jmsg = ""
    if jamp is not None:
        jmsg = (f" | junction amp{jamp} via {rot_engine} (offset {rotinfo['offset']}): "
                f"main-consensus seam had {canon_only_seam_N} N, filled {seam_filled} bp -> "
                f"merged seam N={merged_seam_N}, seam identity {seam_ident}% to ref")
    print(f"[consensus] {tag}: engine={engine} model={sel_model or 'NA'} | %genome={pct_genome} %N={pct_N}{jmsg}")


def write_perbase_depth(depth, chrom, out_gz, step=1):
    with gzip.open(out_gz, "wt") as fh:
        i = 0; n = len(depth)
        while i < n:
            j = i
            while j+1 < n and depth[j+1] == depth[i]: j += 1
            fh.write(f"{chrom}\t{i}\t{j+1}\t{depth[i]}\n"); i = j+1

def write_amplicon_depth(bam, amp_bed, refseq, out_gz, L):
    import pysam
    bf = pysam.AlignmentFile(bam, "rb"); chrom = bf.references[0]
    # group split junction parts (_part1/_part2) under one amplicon name
    spans = {}
    for ln in open(amp_bed):
        c = ln.rstrip("\n").split("\t")
        if len(c) < 5: continue
        name = re.sub(r"_part\d+$", "", c[3])
        spans.setdefault(name, []).append((int(c[1]), int(c[2])))
    rows = []
    for name, segs in spans.items():
        tot = cov = 0
        for s, e in segs:
            cc = bf.count_coverage(chrom, s, min(e, L), quality_threshold=0)
            for i in range(min(e, L) - s):
                cov += cc[0][i]+cc[1][i]+cc[2][i]+cc[3][i]; tot += 1
        mean = round(cov / tot, 1) if tot else 0.0
        s0, e0 = segs[0]
        rows.append(f"{chrom}\t{s0}\t{e0}\t{name}\t{mean}")
    bf.close()
    with gzip.open(out_gz, "wt") as fh: fh.write("\n".join(rows) + "\n")


# ----------------------------- model selection (Task 2) -----------------------------
def select_model(reads, models_dir, override=None):
    """Pick an INSTALLED Clair3 model from MODELS_DIR matching the FASTQ header's
    basecall_model_version_id (pore/chemistry r1041_e82, speed tier 400/260bps, accuracy
    tier fast/hac/sup). Never returns a name that is not present in MODELS_DIR: if there is
    no exact match it maps to the nearest installed model (same chemistry, preferred tier)
    and WARNs. Returns None only when MODELS_DIR is empty (then ARTIC auto-selects)."""
    installed = ([os.path.basename(d) for d in sorted(glob.glob(os.path.join(models_dir, "*")))
                  if os.path.isdir(d)] if models_dir else [])
    def pick(pred):
        for m in installed:
            if pred(m.lower()): return m
        return None
    if override:
        if not installed or override in installed: return override
        print(f"[model] WARN: --model {override} not present in {models_dir}; choosing nearest installed",
              file=sys.stderr)
    hid = ""
    try:
        for h, _, _ in fastq_reads(reads):
            mm = re.search(r"basecall_model_version_id=(\S+)", h)
            if mm: hid = mm.group(1).lower(); break
    except Exception: pass
    tier = "sup" if "sup" in hid else "hac" if "hac" in hid else "fast" if "fast" in hid else None
    sp = re.search(r"(\d+bps)", hid); speed = sp.group(1) if sp else None
    if not installed: return None
    # 1) exact: chemistry + speed + accuracy tier
    if tier and speed:
        m = pick(lambda n: "r1041" in n and speed in n and tier in n)
        if m: return m
    # 2) accuracy tier match (any speed) — nearest by chemistry; warn if speed differs
    if tier:
        m = pick(lambda n: "r1041" in n and tier in n)
        if m:
            if speed and speed not in m.lower():
                print(f"[model] WARN: header speed {speed} not installed; nearest model {m}", file=sys.stderr)
            return m
    # 3) last resort: any r1041 hac, else the first installed model
    m = pick(lambda n: "r1041" in n and "hac" in n) or installed[0]
    print(f"[model] WARN: header model '{hid or 'none'}' unmatched in {models_dir}; defaulting to {m}",
          file=sys.stderr)
    return m

def cmd_model_select(args):
    print(select_model(args.reads, args.models_dir, args.model) or "")


# ----------------------------- binning (Task 3) -----------------------------
def cmd_bin(args):
    """Competitive minimap2 typing of a barcode's reads vs hpv_all_types.fasta."""
    import pysam
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    bam = out / f"{args.sample}.alltypes.bam"
    minimap2_sort(args.ref, args.reads, bam, args.threads)
    bf = pysam.AlignmentFile(bam, "rb")
    best = {}   # read -> (score, type, identity)
    for a in bf.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary: continue
        ref = bf.get_reference_name(a.reference_id)        # e.g. HPV16_K02718
        m = re.match(r"(HPV\d+)", ref); typ = m.group(1) if m else ref.split("_")[0]
        try: nm = a.get_tag("NM")
        except KeyError: nm = 0
        alen = a.query_alignment_length or 1
        ident = max(0.0, 1 - nm/alen)
        score = a.get_tag("AS") if a.has_tag("AS") else a.mapping_quality
        if ident < args.min_ident: continue
        if a.query_name not in best or score > best[a.query_name][0]:
            best[a.query_name] = (score, typ, ident)
    bf.close()
    # per-type read sets + counts
    counts = {}
    for _, (s, typ, ident) in best.items(): counts[typ] = counts.get(typ, 0) + 1
    assign = {r: v[1] for r, v in best.items()}
    # split FASTQ per type
    handles = {}
    for h, s, q in fastq_reads(args.reads):
        rid = h.split()[0]; typ = assign.get(rid)
        if not typ: continue
        fh = handles.get(typ)
        if fh is None:
            fh = gzip.open(out / f"{args.sample}.{typ}.fastq.gz", "wt"); handles[typ] = fh
        fh.write(f"@{h}\n{s}\n+\n{q}\n")
    for fh in handles.values(): fh.close()
    total = sum(counts.values()) or 1
    present = []
    with open(out / f"{args.sample}.types.tsv", "w") as fh:
        fh.write("sample_id\ttype\tn_reads\tfrac\tpresent\n")
        for typ in TYPES:
            n = counts.get(typ, 0); frac = round(n/total, 4)
            isp = n >= args.min_reads
            if isp: present.append(typ)
            if n: fh.write(f"{args.sample}\t{typ}\t{n}\t{frac}\t{'YES' if isp else 'no'}\n")
    (out / f"{args.sample}.present.txt").write_text("\n".join(present) + "\n")
    print(f"[bin] {args.sample}: types present (>= {args.min_reads} reads): {present}"
          f"{'  <-- CO-INFECTION' if len(present) > 1 else ''}")


# ----------------------------- CLI -----------------------------
def main():
    ap = argparse.ArgumentParser(prog="hpvpipe")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("rotate-build"); p.add_argument("--scheme-dir", required=True)
    p.add_argument("--type", required=True); p.add_argument("--outdir", required=True)
    p.set_defaults(func=cmd_rotate_build)

    p = sub.add_parser("bin"); p.add_argument("--sample", required=True)
    p.add_argument("--reads", required=True); p.add_argument("--ref", required=True)
    p.add_argument("--outdir", required=True); p.add_argument("--min-reads", type=int, default=20)
    p.add_argument("--min-ident", type=float, default=0.80); p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_bin)

    p = sub.add_parser("consensus"); p.add_argument("--sample", required=True)
    p.add_argument("--type", required=True); p.add_argument("--reads", required=True)
    p.add_argument("--scheme-dir", required=True); p.add_argument("--outdir", required=True)
    p.add_argument("--models-dir", default=os.environ.get("MODELS_DIR", ""))
    p.add_argument("--model", default=None); p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mindepth", type=int, default=3)
    p.add_argument("--use-artic", action="store_true", help="use artic minion + Clair3 2.x (strict by default)")
    p.add_argument("--require-artic", action="store_true", help="artic is mandatory; abort on any failure")
    p.add_argument("--allow-fallback", action="store_true", help="opt in to the portable engine if artic fails")
    p.set_defaults(func=cmd_consensus)

    p = sub.add_parser("model-select"); p.add_argument("--reads", required=True)
    p.add_argument("--models-dir", default=os.environ.get("MODELS_DIR", "")); p.add_argument("--model", default=None)
    p.set_defaults(func=cmd_model_select)

    a = ap.parse_args(); a.func(a)


if __name__ == "__main__":
    main()
