#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import shlex
import csv
import subprocess
import pandas as pd
from Bio import Entrez, SeqIO

# 讀取 Entrez API key（若有設定）
api_key = os.getenv("ENTREZ_KEY", None)
Entrez.email = "none@example.com"
if api_key:
    Entrez.api_key = api_key

def extract_accession(sseqid: str) -> str:
    """從 'gi|...|...|ACCESSION' 擷取 ACCESSION，否則回 None"""
    if not isinstance(sseqid, str):
        return None
    parts = sseqid.split('|')
    return parts[3] if len(parts) >= 4 and parts[0] == 'gi' else None

def download_references(accessions, outdir):
    """下載所有 accession 的 fasta 到 outdir"""
    os.makedirs(outdir, exist_ok=True)
    for acc in accessions:
        out_fa = os.path.join(outdir, f"{acc}.fasta")
        if os.path.exists(out_fa):
            print(f"[SKIP] {acc}.fasta 已存在")
            continue
        print(f"[FETCH] 下載 {acc}")
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc,
                                   rettype="fasta", retmode="text")
            seq = SeqIO.read(handle, "fasta")
            handle.close()
            SeqIO.write(seq, out_fa, "fasta")
        except Exception as e:
            print(f"[ERROR] 下載 {acc} 失敗：{e}")

def map_and_index(accessions, read1, read2, refdir, mapping_dir):
    """對每支 accession 做 mapping、排序、index，並用 samtools coverage 計算 coverage & depth"""
    os.makedirs(mapping_dir, exist_ok=True)
    csv_path = "mapping_metrics.csv"
    # CSV 欄位：accession + samtools coverage 的九個欄位
    fieldnames = [
        "accession",
        "rname", "startpos", "endpos",
        "numreads", "covbases", "coverage",
        "meandepth", "meanbaseq", "meanmapq"
    ]
    # 寫入 CSV header
    with open(csv_path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

    for acc in accessions:
        ref_fa     = os.path.join(refdir, f"{acc}.fasta")
        bam_prefix = os.path.join(mapping_dir, acc)
        sorted_bam = f"{bam_prefix}_aln.sorted.bam"

        # quote paths 用於 shell 語句
        q_ref = shlex.quote(ref_fa)
        q_r1  = shlex.quote(read1)
        q_r2  = shlex.quote(read2)
        q_bam = shlex.quote(sorted_bam)

        # 1) Mapping + sort + index
        if not os.path.exists(sorted_bam):
            print(f"[MAP] 開始 {acc}")
            cmd_map = (
                f"minimap2 -ax sr {q_ref} {q_r1} {q_r2} | "
                f"samtools view -bS - | "
                f"samtools sort -o {q_bam}"
            )
            try:
                subprocess.run(cmd_map, shell=True, check=True)
                subprocess.run(f"samtools index {q_bam}", shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] {acc} mapping 失敗，跳過: {e}")
                continue
        else:
            print(f"[SKIP] 已存在 BAM：{sorted_bam}")

        # 2) 用 samtools coverage 取代原本的 depth/coverage 計算
        print(f"[COVERAGE] 計算 {acc} ...")
        try:
            result = subprocess.run(
                ["samtools", "coverage", sorted_bam],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] samtools coverage 執行失敗 for {acc}：{e.stderr.strip()}")
            continue

        lines = result.stdout.strip().splitlines()
        if len(lines) < 2:
            print(f"[WARNING] {acc} 無 coverage 輸出，跳過")
            continue

        # 解析 header & values
        header = lines[0].lstrip('#').split('\t')
        values = lines[1].split('\t')
        if len(header) != len(values):
            print(f"[ERROR] {acc} 欄位數量不符：取得 {len(values)} 欄，預期 {len(header)} 欄，跳過")
            continue

        # coverage 欄位是百分比，若為 0 則刪除 BAM/BAI 並跳過
        try:
            cov_pct = float(values[header.index("coverage")])
        except ValueError:
            cov_pct = 0.0
        if cov_pct == 0.0:
            print(f"[CLEAN] {acc} coverage=0，刪除 BAM/BAI")
            try:
                os.remove(sorted_bam)
                bai = sorted_bam + ".bai"
                if os.path.exists(bai):
                    os.remove(bai)
            except OSError as e:
                print(f"[WARN] 刪除 BAM/BAI 失敗: {e}")
            continue

        # 組成要寫入的 row
        row = {"accession": acc}
        for col_name, col_val in zip(header, values):
            row[col_name] = col_val

        # 寫入 CSV
        with open(csv_path, "a", newline="", encoding="utf-8") as fout:
            writer = csv.DictWriter(fout, fieldnames=fieldnames)
            writer.writerow(row)

        print(f"[WRITE] 已寫入 {acc} coverage metrics")

    print(f"[DONE] 完成所有 accession 處理，輸出檔案：{csv_path}")

def main():
    parser = argparse.ArgumentParser(
        description="自動下載參考、mapping、並用 samtools coverage 計算 coverage & depth"
    )
    parser.add_argument("blast_csv", help="BLAST_parsing.csv 檔案")
    parser.add_argument("read1",     help="read1.fastq(.gz)")
    parser.add_argument("read2",     help="read2.fastq(.gz)")
    parser.add_argument(
        "--refdir", default="NCBI_nt",
        help="下載序列儲存資料夾 (預設: NCBI_nt)"
    )
    parser.add_argument(
        "--mapdir", default="mapping",
        help="mapping 結果資料夾 (預設: mapping)"
    )
    args = parser.parse_args()

    # 絕對路徑化
    base_dir   = os.getcwd()
    refdir_abs = os.path.join(base_dir, args.refdir)
    mapdir_abs = os.path.join(base_dir, args.mapdir)

    # 讀 CSV、擷取 accession、去重並排序
    df   = pd.read_csv(args.blast_csv, dtype=str)
    accs = sorted({
        extract_accession(s)
        for s in df["sseqid"]
        if extract_accession(s)
    })
    print(f"[INFO] 共發現 {len(accs)} 個不重複的 accession")

    # 下載 references、做 mapping & coverage 計算
    download_references(accs, refdir_abs)
    map_and_index(accs, args.read1, args.read2,
                  refdir=refdir_abs, mapping_dir=mapdir_abs)

if __name__ == "__main__":
    main()
