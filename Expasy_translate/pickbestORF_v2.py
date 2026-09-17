#!/usr/bin/env python3
"""
Pick the best VP2 ORF from an ExPASy 6-frame translation FASTA.

Logic:
  1. From all 6 frames, collect candidate ORFs (M...stop or M...end).
  2. Keep only candidates >= 200 aa.
  3. BLAST each candidate against a pre-built VP2 protein BLAST database.
  4. Keep candidates with >= 75% identity to VP2.
  5. Among those, pick the longest as the "best" ORF.
  6. Write results + log.

Usage:
    python pick_best_orf.py my_output.fasta best_orfs.fasta vp2_db [min_len] [min_pident]

Note: vp2_db is the PATH PREFIX of a pre-built blastp database, e.g.:
    makeblastdb -in vp2_ref.fasta -dbtype prot -out refs/blastdb/vp2_db
    then pass vp2_db="refs/blastdb/vp2_db"
"""
import sys
import os
import subprocess
import tempfile
from collections import defaultdict

MIN_LEN_DEFAULT = 200
MIN_PIDENT_DEFAULT = 75.0


def parse_fasta(path):
    records = {}
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records[header] = "".join(seq_lines)
    return records


def all_candidate_orfs(seq):
    """Return all M-start ORFs (pieces bounded by stop '-' or seq end)."""
    candidates = []
    for piece in seq.split("-"):
        m_pos = piece.find("M")
        if m_pos == -1:
            continue
        candidates.append(piece[m_pos:])
    return candidates


def main(infile, outfile, vp2_db, min_len=MIN_LEN_DEFAULT, min_pident=MIN_PIDENT_DEFAULT):
    logfile = outfile + ".log.txt"
    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(msg)

    def write_log_and_exit():
        with open(logfile, "w") as lf:
            lf.write("\n".join(log_lines) + "\n")

    source_name = os.path.splitext(os.path.basename(infile))[0]
    if source_name.endswith("_translated"):
        source_name = source_name[: -len("_translated")]

    log("=" * 80)
    log(f"VP2 ORF picking run")
    log(f"Input file:        {infile}")
    log(f"Output file:       {outfile}")
    log(f"VP2 BLAST db:      {vp2_db}")
    log(f"Min ORF length:    {min_len} aa")
    log(f"Min identity:      {min_pident}%")
    log("=" * 80)
    log()

    # sanity check the db files exist before doing any work
    required_ext = [".phr", ".pin", ".psq"]
    missing = [ext for ext in required_ext if not os.path.exists(vp2_db + ext)]
    if missing:
        log(f"ERROR: BLAST db '{vp2_db}' looks incomplete or missing (no {missing} found).")
        log(f"Build it first with: makeblastdb -in <vp2_ref.fasta> -dbtype prot -out {vp2_db}")
        write_log_and_exit()
        return

    records = parse_fasta(infile)
    samples = defaultdict(list)
    for header, seq in records.items():
        sample_id = header.split(":")[0].strip()
        samples[sample_id].append((header, seq))

    # 1+2. collect candidate ORFs >= min_len across all frames, sorted by length desc
    candidates = {}  # cand_id -> (sample_id, header, orf)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf:
        query_path = qf.name
        for sample_id, frames in samples.items():
            all_cands = []
            for header, seq in frames:
                for i, orf in enumerate(all_candidate_orfs(seq)):
                    if len(orf) >= min_len:
                        all_cands.append((header, i, orf))
            all_cands.sort(key=lambda x: len(x[2]), reverse=True)
            for header, i, orf in all_cands:
                # Use a BLAST-safe ID (no spaces or special characters)
                cand_id = f"{sample_id}_candidate_{i}_{len(candidates)}"
                candidates[cand_id] = (sample_id, header, orf)
                qf.write(f">{cand_id}\n{orf}\n")

    if not candidates:
        log(f"No candidate ORFs >= {min_len} aa found in any sample. Nothing to BLAST.")
        os.remove(query_path)
        write_log_and_exit()
        return

    log(f"Total candidate ORFs (>= {min_len} aa) across all samples: {len(candidates)}")
    log()

    # 3. blast all candidates against the pre-built VP2 protein db
    blast_out = query_path + ".blast.tsv"
    result = subprocess.run([
        "blastp", "-query", query_path, "-db", vp2_db,
        "-word_size", "3", "-evalue", "1e-3",
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore qcovs",
        "-out", blast_out
    ], capture_output=True, text=True)

    if result.returncode != 0:
        log(f"ERROR: blastp failed (exit code {result.returncode})")
        log(f"stderr: {result.stderr.strip()}")
        os.remove(query_path)
        write_log_and_exit()
        return

    # best hit (highest pident) per candidate
    best_hit = {}  # cand_id -> (pident, bitscore, qcovs)
    with open(blast_out) as f:
        for line in f:
            fld = line.rstrip("\n").split("\t")
            cand_id, pident, bitscore, qcovs = fld[0], float(fld[2]), float(fld[5]), float(fld[6])
            if cand_id not in best_hit or pident > best_hit[cand_id][0]:
                best_hit[cand_id] = (pident, bitscore, qcovs)

    if not best_hit:
        log("WARNING: blastp ran successfully but returned ZERO hits for ALL candidates.")
        log("This usually means the db, query, or blastp parameters are misconfigured —")
        log("check the db path and try a manual blastp run to confirm.")
        log()

    # 4+5. filter by min_pident, pick longest among passing candidates, per sample
    best_per_sample = {}  # sample_id -> (length, cand_id, header, orf, pident, bitscore, qcovs)
    for cand_id, (sample_id, header, orf) in candidates.items():
        pident, bitscore, qcovs = best_hit.get(cand_id, (0.0, 0.0, 0.0))
        if pident < min_pident:
            continue
        length = len(orf)
        if sample_id not in best_per_sample or length > best_per_sample[sample_id][0]:
            best_per_sample[sample_id] = (length, cand_id, header, orf, pident, bitscore, qcovs)

    # 6. write output + full log
    n_written = 0
    n_skipped = 0
    with open(outfile, "w") as out:
        for sample_id, frames in samples.items():
            log(f"Sample: {sample_id}")
            sample_cands = {cid: v for cid, v in candidates.items() if v[0] == sample_id}
            if not sample_cands:
                log(f"    no ORF >= {min_len} aa in any of the 6 frames — skipped")
                n_skipped += 1
                log()
                continue
            for cid, (sid, header, orf) in sample_cands.items():
                pident, bitscore, qcovs = best_hit.get(cid, (0.0, 0.0, 0.0))
                passed = "PASS" if pident >= min_pident else "fail"
                log(f"    frame={header}  len={len(orf)}aa  pident={pident:.1f}%  "
                    f"bitscore={bitscore:.1f}  qcov={qcovs:.0f}%  [{passed}]")
            if sample_id not in best_per_sample:
                log(f"    WARNING: no candidate reached {min_pident}% identity to VP2 — skipped")
                n_skipped += 1
                log()
                continue
            length, cand_id, header, orf, pident, bitscore, qcovs = best_per_sample[sample_id]
            log(f"    --> BEST: frame={header}  len={length}aa  pident={pident:.1f}%")
            out.write(f">{source_name} | expasy_id={sample_id} | best_frame={header} | "
                      f"length={length}aa | pident={pident:.1f}% | bitscore={bitscore:.1f}\n")
            out.write(orf + "\n")
            n_written += 1
            log()

    log("=" * 80)
    log(f"Summary: {n_written} sequences written, {n_skipped} samples skipped, "
        f"{len(samples)} total samples processed")
    log(f"Output FASTA: {outfile}")
    log(f"Log file:     {logfile}")
    log("=" * 80)

    os.remove(query_path)
    os.remove(blast_out)

    write_log_and_exit()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python pick_best_orf.py <input_6frame.fasta> <output_best_orfs.fasta> "
              "<vp2_db_prefix> [min_len=200] [min_pident=75]")
        sys.exit(1)
    infile, outfile, vp2_db = sys.argv[1:4]
    min_len = int(sys.argv[4]) if len(sys.argv) > 4 else MIN_LEN_DEFAULT
    min_pident = float(sys.argv[5]) if len(sys.argv) > 5 else MIN_PIDENT_DEFAULT
    main(infile, outfile, vp2_db, min_len, min_pident)
    
    
    
