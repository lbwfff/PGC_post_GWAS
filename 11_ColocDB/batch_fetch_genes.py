#!/usr/bin/env python3
"""
Batch fetch COLOC and SMR data from COLOCdb for a list of genes.

Reads Table1.csv, extracts unique genes with totalscore > 1,
then fetches both COLOC and SMR results for each gene.
Skips genes that already have output files locally.
"""

import requests
import pandas as pd
import time
import os
import argparse

COLOC_API_URL = "https://ngdc.cncb.ac.cn/colocdb/api/coloc"
SMR_API_URL = "https://ngdc.cncb.ac.cn/colocdb/api/smr"


def fetch_coloc_for_gene(gene_id, page_size=100):
    """Fetch COLOC results for a given gene with pagination."""
    all_rows = []
    page_index = 1

    while True:
        params = {
            "pageSize": page_size,
            "pageIndex": page_index,
            "gene_id": gene_id,
        }

        r = requests.get(COLOC_API_URL, params=params, timeout=30)
        r.raise_for_status()
        j = r.json()

        meta = j["meta"]
        total = meta["total"]
        rows = j["data"]

        print(f"  [COLOC] page {page_index}: got {len(rows)} / total {total}")

        if not rows:
            break

        all_rows.extend(rows)

        if page_index * page_size >= total:
            break

        page_index += 1
        time.sleep(0.3)

    return all_rows


def fetch_smr_for_gene(gene, page_size=100):
    """Fetch SMR results for a given gene with pagination."""
    all_records = []
    page_index = 1

    while True:
        params = {
            "pageSize": page_size,
            "pageIndex": page_index,
            "gene": gene,
        }

        resp = requests.get(SMR_API_URL, params=params, timeout=30)
        resp.raise_for_status()
        j = resp.json()

        meta = j.get("meta", {})
        total = meta.get("total", 0)
        records = j.get("data", [])

        print(f"  [SMR] page {page_index}, got {len(records)}, total={total}")

        if not records:
            break

        all_records.extend(records)

        if page_index * page_size >= total:
            break

        page_index += 1
        time.sleep(0.3)

    return all_records


def get_genes_from_table(table_csv, totalscore_threshold=1):
    """
    Read Table1.csv and return a sorted list of unique RMP genes
    where totalscore > totalscore_threshold.
    """
    df = pd.read_csv(table_csv)

    # Ensure columns exist
    if "RMP" not in df.columns or "totalscore" not in df.columns:
        raise ValueError(
            "Table1.csv must contain 'RMP' and 'totalscore' columns. "
            f"Found columns: {list(df.columns)}"
        )

    # Filter by totalscore > threshold
    filtered = df[df["totalscore"] > totalscore_threshold]

    # Get unique gene names, sorted
    genes = sorted(filtered["RMP"].unique())

    print(f"Read {len(df)} rows from {table_csv}")
    print(f"Genes with totalscore > {totalscore_threshold}: {len(genes)}")
    print(f"Gene list: {genes}")

    return genes


def fetch_gene_if_needed(gene, output_dir="."):
    """
    Fetch COLOC and SMR data for a gene, but skip if local files already exist.
    Returns True if any fetching was done, False if skipped entirely.
    """
    coloc_file = os.path.join(output_dir, f"COLOCdb_{gene}.csv")
    smr_file = os.path.join(output_dir, f"COLOCdb_SMR_{gene}.csv")

    coloc_exists = os.path.isfile(coloc_file)
    smr_exists = os.path.isfile(smr_file)

    if coloc_exists and smr_exists:
        print(f"  → Both files already exist, skipping {gene}")
        return False

    if coloc_exists:
        print(f"  → COLOC file exists, fetching SMR only for {gene}")
        _fetch_and_save_smr(gene, output_dir)
        return True

    if smr_exists:
        print(f"  → SMR file exists, fetching COLOC only for {gene}")
        _fetch_and_save_coloc(gene, output_dir)
        return True

    # Neither file exists, fetch both
    print(f"  → Fetching both COLOC and SMR for {gene}")
    _fetch_and_save_coloc(gene, output_dir)
    _fetch_and_save_smr(gene, output_dir)
    return True


def _fetch_and_save_coloc(gene, output_dir):
    """Fetch COLOC data and save to CSV."""
    coloc_file = os.path.join(output_dir, f"COLOCdb_{gene}.csv")
    rows = fetch_coloc_for_gene(gene)

    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(coloc_file, index=False)
        print(f"  ✓ Saved {len(rows)} COLOC rows to {coloc_file}")
    else:
        # Save empty CSV with standard columns
        df = pd.DataFrame(columns=[
            "Gene", "Tissue", "eQTL", "GWAS", "chr", "pos", "ref",
            "alt", "AF", "yue_alpha", "yue_alpha_se", "yue_pval",
            "coloc_abf_pp0", "coloc_abf_pp1", "coloc_abf_pp2",
            "coloc_abf_pp3", "coloc_abf_pp4", "coloc_summary"
        ])
        df.to_csv(coloc_file, index=False)
        print(f"  ✓ No COLOC data, saved empty file to {coloc_file}")


def _fetch_and_save_smr(gene, output_dir):
    """Fetch SMR data and save to CSV."""
    smr_file = os.path.join(output_dir, f"COLOCdb_SMR_{gene}.csv")
    records = fetch_smr_for_gene(gene)

    if records:
        df = pd.DataFrame(records)
        df.to_csv(smr_file, index=False)
        print(f"  ✓ Saved {len(records)} SMR rows to {smr_file}")
    else:
        # Save empty CSV with standard columns
        df = pd.DataFrame(columns=[
            "Gene", "eQTL", "GWAS", "ProbeChr", "ProbePos", "SNP",
            "SNPChr", "SNPPos", "A1", "A2", "Freq", "b", "se",
            "p", "n", "ld_r2", "Direction", "p_SMR", "p_HEIDI"
        ])
        df.to_csv(smr_file, index=False)
        print(f"  ✓ No SMR data, saved empty file to {smr_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch fetch COLOC and SMR data from COLOCdb for genes in Table1.csv"
    )
    parser.add_argument(
        "-t", "--table",
        type=str,
        default="Table1.csv",
        help="Path to Table1.csv (default: Table1.csv)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default=".",
        help="Output directory for CSV files (default: current directory)"
    )
    parser.add_argument(
        "--threshold",
        type=int,
        default=1,
        help="totalscore threshold (default: 1, i.e. totalscore > 1)"
    )

    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Read genes from table
    print("=" * 60)
    print("Step 1: Reading gene list from Table1.csv")
    print("=" * 60)
    genes = get_genes_from_table(args.table, args.threshold)

    if not genes:
        print("No genes found with totalscore > {}. Nothing to do.".format(args.threshold))
        return

    # Step 2: Batch fetch for each gene
    print("\n" + "=" * 60)
    print(f"Step 2: Fetching data for {len(genes)} genes")
    print("=" * 60)

    fetched_count = 0
    skipped_count = 0

    for i, gene in enumerate(genes, 1):
        print(f"\n[{i}/{len(genes)}] Processing gene: {gene}")
        print("-" * 40)

        did_fetch = fetch_gene_if_needed(gene, output_dir)

        if did_fetch:
            fetched_count += 1
        else:
            skipped_count += 1

        # Slight delay between genes to be respectful to the API
        if i < len(genes):
            time.sleep(0.5)

    # Summary
    print("\n" + "=" * 60)
    print("Done!")
    print(f"  Total genes processed: {len(genes)}")
    print(f"  Fetched (or partially fetched): {fetched_count}")
    print(f"  Skipped (all files existed):    {skipped_count}")
    print("=" * 60)


if __name__ == "__main__":
    main()