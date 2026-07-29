#!/usr/bin/env python3
"""
Extract all 132,404 prophage sequences from 26k genome assemblies.

Usage:
    python3 extract_prophages.py

Inputs:
    - /home/erikg/phind/26k_prophage1.csv  (132,404 records)
    - /mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects/  (26k genomes)

Output:
    - prophage_homology_survey/full_prophages.fa
"""

import csv
import os
import subprocess
import sys
import time
import glob
from collections import defaultdict

# Paths
CSV_PATH = "/home/erikg/phind/26k_prophage1.csv"
GENOMES_DIR = "/mnt/nvme3n1/erikg/phind-genome-work/ecoli26k/v1/26k/canonical_objects/"
OUTPUT_DIR = "prophage_homology_survey"
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "full_prophages.fa")
LOG_FILE = os.path.join(OUTPUT_DIR, "extraction_log.txt")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# We'll accumulate all FASTA output and write at the end
# to avoid O(132k) file appends
all_output_lines = []
total_attempted = 0
total_success = 0
total_failed = 0
genome_stats = defaultdict(lambda: {"attempted": 0, "success": 0, "failed": 0})
failed_records = []


def extract_prophage_batch(genome, records):
    """
    Extract all prophage sequences for one genome in a single samtools faidx call.
    records: list of (prophage_id, scaffold, begin, end) tuples
    Returns list of (prophage_id, sequence_bases_or_None) pairs.
    """
    fasta_path = os.path.join(GENOMES_DIR, genome, f"{genome}.pansn.fa.gz")
    
    if not os.path.exists(fasta_path):
        return [(pid, None, f"FASTA not found: {fasta_path}") for pid, _, _, _ in records]

    # Build region arguments for samtools
    # Format: {fasta_seqid}:{begin}-{end}
    # The fasta_seqid follows the pattern: {genome}#1#{scaffold}
    regions = []
    pid_map = {}  # region_string -> prophage_id
    for prophage_id, scaffold, begin, end in records:
        begin_int = int(float(begin))
        end_int = int(float(end))
        fasta_seqid = f"{genome}#1#{scaffold}"
        # Ensure begin < end for samtools
        if begin_int > end_int:
            begin_int, end_int = end_int, begin_int
        # Clamp to 1-based positive coordinates
        if begin_int < 1:
            begin_int = 1
        region = f"{fasta_seqid}:{begin_int}-{end_int}"
        regions.append(region)
        pid_map[region] = prophage_id

    if not regions:
        return [(pid, None, "No valid regions") for pid, _, _, _ in records]

    try:
        # Run samtools faidx once for all regions
        cmd = ["samtools", "faidx", fasta_path] + regions
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )
        
        if result.returncode != 0:
            stderr = result.stderr.strip()[:200]
            return [(pid, None, f"samtools error: {stderr}") for pid, _, _, _ in records]
        
        # Parse the multi-FASTA output
        output = result.stdout
        sequences = parse_fasta_output(output)
        
        # Map regions to sequences
        results = []
        for region in regions:
            prophage_id = pid_map[region]
            seq = sequences.get(region)
            if seq is None:
                # The region might have been in the stderr as a warning
                results.append((prophage_id, None, "No sequence in output"))
            else:
                results.append((prophage_id, seq, None))
        
        return results

    except subprocess.TimeoutExpired:
        return [(pid, None, "samtools timeout") for pid, _, _, _ in records]
    except Exception as e:
        return [(pid, None, f"Exception: {str(e)[:200]}") for pid, _, _, _ in records]


def parse_fasta_output(output):
    """
    Parse samtools faidx output into {region_header: sequence_string}.
    Headers look like: >GCF_000005845.2#1#NC_000913.3:564815-582583
    """
    sequences = {}
    current_header = None
    current_seq_parts = []
    
    for line in output.splitlines():
        if line.startswith(">"):
            # Save previous sequence
            if current_header is not None:
                sequences[current_header] = "".join(current_seq_parts)
            # Extract just the region part (remove the leading '>')
            current_header = line[1:].strip()
            current_seq_parts = []
        else:
            current_seq_parts.append(line.strip())
    
    # Don't forget the last sequence
    if current_header is not None:
        sequences[current_header] = "".join(current_seq_parts)
    
    return sequences


def main():
    global total_attempted, total_success, total_failed
    
    print(f"Reading prophage CSV: {CSV_PATH}")
    start_time = time.time()
    
    # Read CSV and group by genome
    genome_groups = defaultdict(list)
    header_checked = False
    
    with open(CSV_PATH, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not header_checked:
                header_checked = True
                # Verify expected columns
                expected = ['end', 'genome', 'scaffold', 'begin', 'transposable', 'taxonomy', 'prophage_id']
                if row != expected:
                    print(f"Warning: Unexpected header: {row}")
                continue
            
            if len(row) < 7:
                continue
            
            end, genome, scaffold, begin, transposable, taxonomy, prophage_id = row[:7]
            genome_groups[genome].append((prophage_id, scaffold, begin, end))
    
    total_genomes = len(genome_groups)
    total_records = sum(len(records) for records in genome_groups.values())
    print(f"Loaded {total_records} records across {total_genomes} genomes")
    
    # Process genomes one at a time
    processed = 0
    genome_success = 0
    genome_fail = 0
    
    # Write output incrementally to avoid memory issues with 132k sequences
    with open(OUTPUT_FASTA, 'w') as out_f:
        with open(LOG_FILE, 'w') as log_f:
            log_f.write(f"Extraction log for {total_records} prophage records\n")
            log_f.write(f"Started: {time.ctime()}\n")
            log_f.write("=" * 60 + "\n")
            
            for genome, records in sorted(genome_groups.items()):
                processed += 1
                genome_stats[genome]["attempted"] = len(records)
                total_attempted += len(records)
                
                if processed % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = processed / elapsed if elapsed > 0 else 0
                    print(f"  [{processed}/{total_genomes}] genomes processed ({processed/elapsed:.1f}/s) - "
                          f"{total_success}/{total_attempted} records OK")
                
                results = extract_prophage_batch(genome, records)
                
                for prophage_id, seq, error in results:
                    if seq is not None:
                        out_f.write(f">{prophage_id}\n{seq}\n")
                        total_success += 1
                        genome_stats[genome]["success"] += 1
                    else:
                        total_failed += 1
                        genome_stats[genome]["failed"] += 1
                        failed_records.append((prophage_id, genome, error))
                        log_f.write(f"FAIL: {prophage_id} | {genome} | {error}\n")
    
    elapsed = time.time() - start_time
    success_rate = (total_success / total_attempted * 100) if total_attempted > 0 else 0
    
    print(f"\n{'=' * 60}")
    print(f"Extraction complete!")
    print(f"  Total attempted: {total_attempted}")
    print(f"  Total success:   {total_success}")
    print(f"  Total failed:    {total_failed}")
    print(f"  Success rate:    {success_rate:.1f}%")
    print(f"  Time elapsed:    {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  Output:          {OUTPUT_FASTA}")
    print(f"  Log:             {LOG_FILE}")
    
    # Summary by genome success rate
    genomes_with_failures = [g for g, s in genome_stats.items() if s["failed"] > 0]
    if genomes_with_failures:
        print(f"\nGenomes with failures ({len(genomes_with_failures)}):")
        for g in sorted(genomes_with_failures)[:20]:
            s = genome_stats[g]
            print(f"  {g}: {s['success']}/{s['attempted']} OK")
        if len(genomes_with_failures) > 20:
            print(f"  ... and {len(genomes_with_failures) - 20} more")
    
    # Validation
    print(f"\n{'=' * 60}")
    print(f"Validation:")
    print(f"  All 132,404 records attempted: {'YES' if total_attempted == 132404 else f'NO ({total_attempted})'}")
    print(f"  Success rate > 90%: {'YES' if success_rate > 90 else 'NO'}")
    print(f"  FASTA lines written: {total_success}")
    
    # Check FASTA format
    fasta_lines = 0
    fasta_headers = 0
    with open(OUTPUT_FASTA, 'r') as f:
        for line in f:
            fasta_lines += 1
            if line.startswith('>'):
                fasta_headers += 1
                # Check single-line header
                if len(line.strip().split()) > 1:
                    print(f"  WARNING: Multi-token header: {line.strip()[:80]}")
    
    print(f"  FASTA total lines: {fasta_lines}")
    print(f"  FASTA headers: {fasta_headers}")
    print(f"  All single-line headers: {'YES' if fasta_headers == total_success else 'NO'}")


if __name__ == "__main__":
    main()