#!/usr/bin/env python3
"""
ortholog_finder.py
------------------
Automate ortholog discovery across Costus genomes:
  1) makeblastdb for each species
  2) blastn for candidate gene FASTAs (from C. lasius)
  3) parse tabular hits, pick best per query per species
  4) extract regions with gffread (-w and -y)

Requires:
  - Python 3.8+
  - PyYAML (pip install pyyaml)
  - External tools in PATH or passed via config:
      makeblastdb, blastn, gffread

Usage:
  python3 ortholog_finder.py --config configs/orthologs.yaml
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import yaml
from collections import defaultdict, namedtuple

Hit = namedtuple("Hit", "qseqid sseqid pident length evalue bits sstart send sstrand")

def run(cmd, cwd=None):
    print("[RUN]", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def parse_blast_table(tsv_path):
    """
    Expect outfmt 6 with fields:
      qseqid sseqid pident length evalue bitscore sstart send sstrand
    """
    hits = []
    with open(tsv_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                # silently skip malformed rows
                continue
            qseqid, sseqid, pident, length, evalue, bits, sstart, send, sstrand = parts[:9]
            try:
                hits.append(
                    Hit(
                        qseqid=qseqid,
                        sseqid=sseqid,
                        pident=float(pident),
                        length=int(length),
                        evalue=float(evalue),
                        bits=float(bits),
                        sstart=int(sstart),
                        send=int(send),
                        sstrand=sstrand
                    )
                )
            except ValueError:
                # skip if conversion fails
                continue
    return hits

def choose_best_hits(hits, min_pident, min_align_len, max_evalue):
    """
    Return best hit per query (lowest evalue, then highest bits) after filtering.
    """
    by_query = defaultdict(list)
    for h in hits:
        if h.pident >= min_pident and h.length >= min_align_len and h.evalue <= max_evalue:
            by_query[h.qseqid].append(h)

    best = {}
    for q, hs in by_query.items():
        hs.sort(key=lambda x: (x.evalue, -x.bits, -x.pident, -x.length))
        best[q] = hs[0]
    return best  # dict: query_id -> Hit

def safe_region(start, end, pad, chrom_len=None):
    lo = min(start, end) - pad
    hi = max(start, end) + pad
    lo = 1 if lo < 1 else lo
    if chrom_len is not None:
        hi = chrom_len if hi > chrom_len else hi
    return lo, hi

def extract_with_gffread(gffread_cmd, genome_fasta, gff3, chrom, start, end, out_prefix):
    # nucleotide
    nt_out = f"{out_prefix}_transcripts.fa"
    aa_out = f"{out_prefix}_translated.fa"
    region = f"{chrom}:{start}-{end}"

    run([gffread_cmd, "-w", nt_out, "-g", genome_fasta, "-r", region, gff3])
    run([gffread_cmd, "-y", aa_out, "-g", genome_fasta, "-r", region, gff3])
    return nt_out, aa_out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="YAML config file")
    args = ap.parse_args()

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    blastn_cmd     = cfg.get("tools", {}).get("blastn", "blastn")
    makeblastdb_cmd= cfg.get("tools", {}).get("makeblastdb", "makeblastdb")
    gffread_cmd    = cfg.get("tools", {}).get("gffread", "gffread")

    queries        = cfg["queries"]            # list of FASTA paths for candidate genes (lasius)
    species        = cfg["species"]            # dict of species with genome/gff/db_prefix
    outdir         = Path(cfg.get("outdir", "results/orthologs"))
    ensure_dir(outdir)

    # BLAST params
    evalue_thresh  = float(cfg.get("filters", {}).get("max_evalue", 1e-10))
    min_pident     = float(cfg.get("filters", {}).get("min_pident", 50.0))
    min_align_len  = int(cfg.get("filters", {}).get("min_align_len", 100))
    pad_bp         = int(cfg.get("region_padding_bp", 1000))
    keep_tables    = bool(cfg.get("keep_blast_tables", True))

    # 1) makeblastdb for each species
    for sp_name, meta in species.items():
        genome = meta["genome_fasta"]
        db_pref= meta.get("db_prefix", str(Path(genome).with_suffix("")).replace(".fa","").replace(".fna","") + "_db")
        meta["db_prefix"] = db_pref  # store back
        ensure_dir(Path(db_pref).parent if "/" in db_pref else Path("."))

        # create DB if not present (naive check on .nhr file)
        if not Path(db_pref + ".nhr").exists():
            run([makeblastdb_cmd, "-in", genome, "-dbtype", "nucl", "-out", db_pref])
        else:
            print(f"[SKIP] DB exists for {sp_name}: {db_pref}")

    # 2) blastn each query against each species DB
    # outfmt includes sstart/send/sstrand so we can compute regions
    outfmt = "6 qseqid sseqid pident length evalue bitscore sstart send sstrand"
    blast_tables = []

    for qf in queries:
        qname = Path(qf).stem
        for sp_name, meta in species.items():
            db = meta["db_prefix"]
            sp_tag = sp_name.replace(" ", "_")
            tbl = outdir / f"{qname}__vs__{sp_tag}.blast.tsv"
            blast_tables.append((sp_name, qname, tbl))
            run([
                blastn_cmd, "-query", qf, "-db", db,
                "-outfmt", outfmt, "-evalue", str(evalue_thresh),
                "-max_target_seqs", "10", "-num_threads", str(cfg.get("threads", 1)),
                "-out", str(tbl)
            ])

    # 3) pick best hits per query per species and 4) extract with gffread
    summary_rows = []
    for sp_name, qname, tbl_path in blast_tables:
        hits = parse_blast_table(tbl_path)
        best = choose_best_hits(hits, min_pident, min_align_len, evalue_thresh)
        meta = species[sp_name]
        genome = meta["genome_fasta"]
        gff3   = meta["gff3"]

        sp_outdir = outdir / sp_name.replace(" ", "_") / qname
        ensure_dir(sp_outdir)

        if not keep_tables:
            # move or delete tables if desired; for now leave them
            pass

        if len(best) == 0:
            print(f"[WARN] No passing hits for {qname} in {sp_name}")
            continue

        for qid, h in best.items():
            # Compute region (strand-aware; pad)
            start, end = safe_region(h.sstart, h.send, pad_bp)
            out_prefix = sp_outdir / f"{qid}__{h.sseqid}_{start}-{end}"
            nt_out, aa_out = extract_with_gffread(
                gffread_cmd=gffread_cmd,
                genome_fasta=genome,
                gff3=gff3,
                chrom=h.sseqid,
                start=start,
                end=end,
                out_prefix=str(out_prefix)
            )
            summary_rows.append({
                "species": sp_name,
                "query_id": qid,
                "subject": h.sseqid,
                "start": start,
                "end": end,
                "strand": h.sstrand,
                "pident": h.pident,
                "length": h.length,
                "evalue": h.evalue,
                "bitscore": h.bits,
                "nt_fasta": nt_out,
                "aa_fasta": aa_out,
                "blast_table": str(tbl_path)
            })

    # 5) write a summary TSV
    import csv
    summary_path = outdir / "ortholog_summary.tsv"
    if summary_rows:
        keys = list(summary_rows[0].keys())
        with open(summary_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
            w.writeheader()
            for row in summary_rows:
                w.writerow(row)
        print(f"[DONE] Summary: {summary_path}")
    else:
        print("[DONE] No hits written; summary not created.")

if __name__ == "__main__":
    main()

