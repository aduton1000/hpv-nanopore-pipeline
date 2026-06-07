#!/usr/bin/env python
"""
convert_scheme_to_artic.py — upgrade an ALREADY-INSTALLED primer-schemes/ tree in place
to the ARTIC 1.10.1 layout, WITHOUT a full re-ingest. Idempotent.

For each <scheme-dir>/HPV{t}/V5/ it:
  1. fixes the chrom/header  HPV{t}|{acc} -> HPV{t}_{acc}  (pipe is invalid for
     primalbedtools/ARTIC chrom ^[a-zA-Z0-9_.]+$) in reference.fasta + *.bed;
  2. writes the legacy name get_scheme_legacy() expects: HPV{t}.scheme.bed;
  3. creates the ARTIC 1.10.1 new layout HPV{t}/<length>/<version>/{primer.bed,reference.fasta};
  4. (after all types) writes a local <scheme-dir>/index.json manifest so get_scheme()
     resolves OFFLINE (falls back to local index.json, then skips download because the
     files already exist).

The analysis pipeline calls `artic minion --bed/--ref` (no manifest needed); this converter
additionally enables the `--scheme-name/--scheme-length/--scheme-version` form offline.

Usage (admin, on the cluster — schemes live under the shared install):
  python convert_scheme_to_artic.py --scheme-dir /hpc/opt/hpv-bioinfo/primer-schemes
  python convert_scheme_to_artic.py --scheme-dir ... --check     # dry-run
"""
import argparse, glob, hashlib, json, os, re, sys
from pathlib import Path

ARTIC_LEN = "2000"; ARTIC_VER = "v5.0.0"; VERSION = "V5"
CHROM_RE = re.compile(r"^[A-Za-z0-9_.]+$")

def fix_chrom(name):  # HPV16|K02718 -> HPV16_K02718 (only touch a pipe)
    return name.replace("|", "_")

def read_fasta1(p):
    hdr = None; seq = []
    for ln in open(p):
        if ln.startswith(">"):
            hdr = ln[1:].split()[0]
        else: seq.append(ln.strip())
    return hdr, "".join(seq)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scheme-dir", required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    SD = Path(a.scheme_dir)
    if not SD.is_dir(): sys.exit(f"scheme-dir not found: {SD}")
    manifest = {"primerschemes": {}}
    changed = 0
    for vdir in sorted(SD.glob("HPV*/V5")):
        ty = vdir.parent.name
        ref_fa = vdir / f"{ty}.reference.fasta"
        primer_bed = vdir / f"{ty}.primer.bed"
        if not ref_fa.exists() or not primer_bed.exists():
            print(f"  [skip] {ty}: missing reference/primer.bed in {vdir}"); continue
        hdr, seq = read_fasta1(ref_fa); new_hdr = fix_chrom(hdr)
        # 1) fix chrom in reference + all beds
        for f in [ref_fa] + list(vdir.glob("*.bed")):
            txt = f.read_text()
            if "|" in txt:
                if f.suffix == ".fasta":
                    txt = re.sub(r"^>.*$", f">{new_hdr}", txt, count=1, flags=re.M)
                else:
                    txt = "\n".join(("\t".join([fix_chrom(c.split('\t')[0])] + c.split('\t')[1:]) if c.strip() else c)
                                    for c in txt.splitlines()) + ("\n" if txt.endswith("\n") else "")
                if not a.check: f.write_text(txt)
                changed += 1
                print(f"  [fix] {f.relative_to(SD)}: chrom -> {new_hdr}")
        if not CHROM_RE.match(new_hdr):
            print(f"  [WARN] {ty}: chrom '{new_hdr}' still invalid", file=sys.stderr)
        primer_txt = primer_bed.read_text()
        ref_txt = f">{new_hdr}\n" + "".join(seq[i:i+70] + "\n" for i in range(0, len(seq), 70))
        # 2) legacy .scheme.bed
        if not a.check: (vdir / f"{ty}.scheme.bed").write_text(primer_txt)
        # 3) ARTIC 1.10.1 new layout
        nd = SD / ty / ARTIC_LEN / ARTIC_VER
        if not a.check:
            nd.mkdir(parents=True, exist_ok=True)
            (nd / "primer.bed").write_text(primer_txt)
            (nd / "reference.fasta").write_text(ref_txt)
        print(f"  [layout] {ty}: -> {ty}/{ARTIC_LEN}/{ARTIC_VER}/{{primer.bed,reference.fasta}} (+ {ty}.scheme.bed)")
        manifest["primerschemes"].setdefault(ty, {}).setdefault(ARTIC_LEN, {})[ARTIC_VER] = {
            "primer_bed_url": f"https://local.invalid/{ty}/{ARTIC_LEN}/{ARTIC_VER}/primer.bed",
            "primer_bed_md5": hashlib.md5(primer_txt.encode()).hexdigest(),
            "reference_fasta_url": f"https://local.invalid/{ty}/{ARTIC_LEN}/{ARTIC_VER}/reference.fasta",
            "reference_fasta_md5": hashlib.md5(ref_txt.encode()).hexdigest(),
            "status": "draft", "schemeName": ty, "ampliconSize": int(ARTIC_LEN),
            "schemeVersion": ARTIC_VER}
    # 4) local manifest
    if not a.check:
        (SD / "index.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"\n{'DRY-RUN: ' if a.check else ''}upgraded {len(manifest['primerschemes'])} schemes; "
          f"wrote {SD/'index.json'} ({changed} chrom fixes).")
    print("Validate offline:  conda run -n hpv-artic artic minion --bed "
          f"{SD}/HPV16/V5/HPV16.primer.bed --ref {SD}/HPV16/V5/HPV16.reference.fasta "
          "--model-dir $MODELS_DIR --read-file <reads> TEST   # (no manifest fetch)")

if __name__ == "__main__":
    main()
