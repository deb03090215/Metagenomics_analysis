#!/usr/bin/env python3

import os
import sys
import time
import subprocess
import argparse
import csv
import shutil
import shlex

import pandas as pd
from Bio import SeqIO, Entrez


def rename_orfs(gtf_path="orfs.gtf",
                fasta_in="orfs.fasta",
                fasta_out="orfs_renamed.fasta"):
    """
    Rename ORF sequences based on contig IDs in the GTF.
    """
    orf_to_contig = {}
    with open(gtf_path) as gtf:
        for line in gtf:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            contig = cols[0]
            attrs = cols[8]
            parts = [x for x in attrs.split(";") if "orf_id" in x]
            if parts:
                orf_id = parts[0].split("orf_id")[1].replace('"', "").strip()
                orf_to_contig[orf_id] = contig

    records = []
    for rec in SeqIO.parse(fasta_in, "fasta"):
        if rec.id in orf_to_contig:
            rec.id = f"{rec.id}_{orf_to_contig[rec.id]}"
            rec.description = ""
            records.append(rec)
        else:
            print(f"[warn] {rec.id} not in GTF, skipping")

    SeqIO.write(records, fasta_out, "fasta")
    print(f"[info] Renamed ORFs saved to {fasta_out}")


def download_fastas(csv_path, outdir="NCBI_aa", api_key=None, delay=0.3):
    """
    Download one protein FASTA per accession in the CSV's 'sseqid' column.
    """
    Entrez.email = "none@example.com"
    if api_key:
        Entrez.api_key = api_key

    df = pd.read_csv(csv_path, dtype=str)
    accs = sorted(df.get("sseqid", pd.Series()).dropna().unique())
    os.makedirs(outdir, exist_ok=True)
    failed = []

    for acc in accs:
        fasta_path = os.path.join(outdir, f"{acc}.fasta")
        need_fetch = True
        if os.path.exists(fasta_path):
            try:
                recs = list(SeqIO.parse(fasta_path, "fasta"))
                if len(recs) == 1:
                    print(f"[skip] {acc} → {recs[0].id}")
                    need_fetch = False
                else:
                    print(f"[warn] {acc}.fasta has {len(recs)} records, re-downloading")
            except Exception as e:
                print(f"[warn] Cannot parse {acc}.fasta ({e}), re-downloading")
        if need_fetch:
            print(f"[fetch] Downloading {acc} ...")
            try:
                with Entrez.efetch(db="protein", id=acc,
                                  rettype="fasta", retmode="text") as handle:
                    data = handle.read()
                with open(fasta_path, "w") as fw:
                    fw.write(data)
                recs = list(SeqIO.parse(fasta_path, "fasta"))
                if len(recs) == 1:
                    print(f"[saved] {acc} → {recs[0].id}")
                else:
                    print(f"[error] {acc}.fasta has {len(recs)} records")
                    failed.append(acc)
            except Exception as e:
                print(f"[error] Failed to download {acc}: {e}")
                failed.append(acc)
            time.sleep(delay)

    if failed:
        print("\n[summary] The following accessions failed or need manual checking:")
        for a in failed:
            print(" ", a)
    else:
        print("\n[summary] All reference FASTAs OK.")


def run_hmmscan(pfam_hmm: str, fasta: str, cpu: int, workdir: str = None):
    tblout = "pfam_results.tblout"
    domtblout = "pfam.domtblout"
    alignment = "pfam_alignment.txt"

    cmd1 = ["hmmscan", f"--cpu={cpu}", "--tblout", tblout, pfam_hmm, fasta]
    subprocess.run(cmd1, check=True, cwd=workdir)

    cmd2 = ["hmmscan", f"--cpu={cpu}", "--domtblout", domtblout, pfam_hmm, fasta]
    with open(os.path.join(workdir or "", alignment), "w") as aln_fh:
        subprocess.run(cmd2, stdout=aln_fh, check=True, cwd=workdir)


def tblout_to_csv(tblout_path: str, csv_path: str = "pfam_results.csv", workdir: str = None):
    header = [
        "target name","accession","query name","accession",
        "E-value","score","bias","E-value","score","bias",
        "exp","reg","clu","ov","env","dom","rep","inc","description of target"
    ]
    tbl_full = os.path.join(workdir or "", tblout_path)
    csv_full = os.path.join(workdir or "", csv_path)
    with open(tbl_full, 'r') as fin, open(csv_full, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(header)
        for line in fin:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            writer.writerow(parts)


def process_diamond(csv_path, contig_fasta, pfam_hmm, cpu, paired1, paired2):
    """
    Full pipeline: DIAMOND parsing → ORF finder → MAFFT/hmmscan → mapping → coverage.
    """
    base_dir = os.getcwd()
    df = pd.read_csv(csv_path)
    contigs = SeqIO.index(contig_fasta, "fasta")

    for virus in df["sscinames"].dropna().unique():
        virus_dir = os.path.join(base_dir, virus.replace(" ", "_"))
        os.makedirs(virus_dir, exist_ok=True)
        subset = df[df["sscinames"] == virus]

        for acc in subset["sseqid"].dropna().unique():
            acc_safe = acc.replace("/", "_")
            acc_dir = os.path.join(virus_dir, acc_safe)
            os.makedirs(acc_dir, exist_ok=True)

            os.chdir(acc_dir)
            try:
                # 1) Separate forward/reverse hits
                hits = subset[subset["sseqid"] == acc]
                fw_ids, rev_ids = set(), set()
                for _, row in hits.iterrows():
                    qid = row.get("qseqid")
                    if pd.isna(qid):
                        continue
                    if row["qstart"] > row["qend"]:
                        rev_ids.add(qid)
                    else:
                        fw_ids.add(qid)

                with open("forward_contig_list.txt", "w") as f:
                    f.write("\n".join(sorted(fw_ids)))
                if rev_ids:
                    with open("reverse_contig_list.txt", "w") as f:
                        f.write("\n".join(sorted(rev_ids)))

                # 2) Extract contigs & combine
                SeqIO.write((contigs[q] for q in fw_ids if q in contigs),
                            "forward_contigs.fasta", "fasta")
                if rev_ids:
                    rev_recs = []
                    for q in rev_ids:
                        if q in contigs:
                            rec = contigs[q]
                            rec.seq = rec.seq.reverse_complement()
                            rec.id = rec.id + "_rev"
                            rec.description = ""
                            rev_recs.append(rec)
                    SeqIO.write(rev_recs, "reverse_contigs.fasta", "fasta")

                with open("combine.fasta", "w") as out:
                    if os.path.exists("forward_contigs.fasta"):
                        out.write(open("forward_contigs.fasta").read())
                    if rev_ids:
                        out.write(open("reverse_contigs.fasta").read())

                # 3) ORF finder & rename
                subprocess.run("orffinder-to-gtf -in combine.fasta -out orfs.gtf -orf_size 300",
                               shell=True, check=True)
                subprocess.run("orffinder-to-sequence -in combine.fasta -out orfs.fasta"
                               " -outtype protein -orf_size 300",
                               shell=True, check=True)
                rename_orfs()

                # 檢查 orfs_renamed.fasta
                orfs_file = "orfs_renamed.fasta"
                need_alignment = False
                if os.path.exists(orfs_file):
                    if os.path.getsize(orfs_file) > 0:
                        need_alignment = True
                    else:
                        print(f"[warn] {orfs_file} is empty → removing and skipping alignment/hmmscan")
                        os.remove(orfs_file)
                else:
                    print(f"[warn] {orfs_file} not found → skipping alignment/hmmscan")

                if need_alignment:
                    # 4) MAFFT alignment
                    ref_fa = os.path.join(base_dir, "NCBI_aa", f"{acc}.fasta")
                    mafft_cmd = (
                        f"mafft --addfull orfs_renamed.fasta --keeplength {shlex.quote(ref_fa)}"
                        f" > combine_aa_aln.fasta"
                    )
                    subprocess.run(mafft_cmd, shell=True, check=True)

                    # 5) hmmscan
                    run_hmmscan(pfam_hmm, "orfs_renamed.fasta", cpu, workdir=acc_dir)
                    tblout_to_csv("pfam_results.tblout", "pfam_results.csv", workdir=acc_dir)

                # 6) Read mapping
                mapping_dir = os.path.join(acc_dir, "mapping")
                os.makedirs(mapping_dir, exist_ok=True)
                for rec in SeqIO.parse("combine.fasta", "fasta"):
                    name = "_".join(rec.id.split("_")[:2])
                    single_fa = os.path.join(mapping_dir, f"{name}.fasta")
                    SeqIO.write(rec, single_fa, "fasta")

                    bam = os.path.join(mapping_dir, f"{name}_aln.sorted.bam")
                    cmd_map = (
                        f"minimap2 -ax sr {shlex.quote(single_fa)} {shlex.quote(paired1)} {shlex.quote(paired2)}"
                        f" | samtools view -bS - | samtools sort -o {shlex.quote(bam)}"
                    )
                    subprocess.run(cmd_map, shell=True, check=True)
                    subprocess.run(f"samtools index {shlex.quote(bam)}", shell=True, check=True)

                # 7) Coverage metrics
                metrics_file = os.path.join(mapping_dir, "mapping_metrics.csv")
                with open(metrics_file, "w", newline="") as mh:
                    writer = csv.writer(mh)
                    writer.writerow([
                        "accession", "rname", "startpos", "endpos",
                        "numreads", "covbases", "coverage",
                        "meandepth", "meanbaseq", "meanmapq"
                    ])
                    for bam_fn in sorted(os.listdir(mapping_dir)):
                        if not bam_fn.endswith(".bam"):
                            continue
                        bam_path = os.path.join(mapping_dir, bam_fn)
                        cov = subprocess.run(
                            f"samtools coverage {shlex.quote(bam_path)}",
                            shell=True, capture_output=True, text=True, check=True
                        )
                        for line in cov.stdout.splitlines():
                            if line and not line.startswith("#"):
                                cols = line.split()
                                writer.writerow([bam_fn.replace("_aln.sorted.bam","")] + cols)
                                break

                print(f"[done] Processed {virus}/{acc_safe}")

            except subprocess.CalledProcessError as e:
                print(f"[error] Step failed for {acc_safe}: {e}")
            except Exception as e:
                print(f"[error] Unexpected error for {acc_safe}: {e}")
            finally:
                os.chdir(base_dir)


def parse_args():
    p = argparse.ArgumentParser(
        description="Process DIAMOND output through ORF finder, MAFFT, hmmscan, mapping, and coverage."
    )
    p.add_argument("csv_file",    help="DIAMOND parsing CSV")
    p.add_argument("contig_file", help="Contigs FASTA")
    p.add_argument("paired1",     help="Paired-end FASTQ #1")
    p.add_argument("paired2",     help="Paired-end FASTQ #2")
    p.add_argument("--pfam-hmm",   required=True, help="Path to Pfam-A.hmm")
    p.add_argument("--cpu",        type=int, default=1, help="Number of CPUs for hmmscan")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not shutil.which("hmmscan"):
        print("Error: hmmscan not found in PATH.", file=sys.stderr)
        sys.exit(1)

    print("[info] Downloading and verifying reference FASTAs...")
    download_fastas(args.csv_file, outdir="NCBI_aa")

    print("[info] Running full pipeline...")
    process_diamond(
        args.csv_file,
        args.contig_file,
        args.pfam_hmm,
        args.cpu,
        args.paired1,
        args.paired2
    )

    print("All done.")
