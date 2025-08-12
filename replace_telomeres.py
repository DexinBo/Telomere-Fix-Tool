#!/usr/bin/env python3
import os
import re
import subprocess
import tempfile
import shutil
import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import concurrent.futures
from typing import Dict, Optional, List
import argparse
import sys

# 配置日志系统
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def find_mummer_path() -> None:
    """Find the path to nucmer, delta-filter, and show-coords."""
    for tool in ['nucmer', 'delta-filter', 'show-coords']:
        if not shutil.which(tool):
            raise EnvironmentError(
                f"'{tool}' not found in PATH. Please ensure MUMmer is installed and accessible."
            )
    logger.info("MUMmer tools (nucmer, delta-filter, show-coords) found in PATH.")

def run_nucmer(ref_fa: str, qry_fa: str, output_prefix: str, min_identity: float = 80, min_length: int = 1000) -> str:
    """
    Runs nucmer alignment and filters the results.
    The ref_fa and qry_fa can now contain multiple sequences.
    Uses list form for subprocess.run for better security.
    """
    delta_file = f"{output_prefix}.delta"
    filtered_delta_file = f"{output_prefix}.filtered.delta"
    coords_file = f"{output_prefix}.coords"

    # nucmer command
    cmd_nucmer = ["nucmer", "--maxmatch", f"--prefix={output_prefix}", ref_fa, qry_fa]
    logger.debug(f"Running nucmer: {' '.join(cmd_nucmer)}")
    subprocess.run(cmd_nucmer, check=True, capture_output=True)

    # delta-filter command
    cmd_filter = ["delta-filter", "-i", str(min_identity), "-l", str(min_length), delta_file]
    logger.debug(f"Running delta-filter: {' '.join(cmd_filter)}")
    with open(filtered_delta_file, "w") as outfile:
        subprocess.run(cmd_filter, check=True, stdout=outfile)

    # show-coords command
    cmd_coords = ["show-coords", "-rclT", filtered_delta_file]
    logger.debug(f"Running show-coords: {' '.join(cmd_coords)}")
    with open(coords_file, "w") as outfile:
        subprocess.run(cmd_coords, check=True, stdout=outfile)

    return coords_file

def parse_coords_for_ends_optimized(coords_file: str) -> tuple[Optional[Dict], Optional[Dict]]:
    """
    Parses alignment results from nucmer run on extracted end sequences.
    Finds the best hit for both 5' and 3' ends based on a combined score
    of identity and alignment length. This function no longer filters by position.
    """
    best_5_end_hit = None
    best_3_end_hit = None
    
    best_5_end_score = -1.0
    best_3_end_score = -1.0

    with open(coords_file) as f:
        for line in f:
            if line.startswith('/') or line.startswith('[') or not line.strip():
                continue
            
            try:
                cols = line.strip().split('\t')
                if len(cols) < 13: continue

                ref_start, ref_end = int(cols[0]), int(cols[1])
                qry_start, qry_end = int(cols[2]), int(cols[3])
                identity = float(cols[6])
                ref_id = cols[11]
                qry_id = cols[12]
                
                alignment_length = ref_end - ref_start + 1
                current_score = identity * alignment_length
                
                current_hit = {
                    "ref_start": ref_start, "ref_end": ref_end,
                    "qry_start": qry_start, "qry_end": qry_end,
                    "identity": identity,
                    "qry_id": qry_id,
                    "ref_id": ref_id
                }

                if ref_id.endswith("_5end"):
                    if current_score > best_5_end_score:
                        best_5_end_score = current_score
                        best_5_end_hit = current_hit
                elif ref_id.endswith("_3end"):
                    if current_score > best_3_end_score:
                        best_3_end_score = current_score
                        best_3_end_hit = current_hit

            except (IndexError, ValueError) as e:
                logger.warning(f"Skipping malformed line in coords file: '{line.strip()}' due to error: {e}")
                continue

    return best_5_end_hit, best_3_end_hit

def process_single_contig(chr_rec: SeqRecord, consensus_fa: str, consensus_seqs: Dict[str, str], min_identity: float, end_region: int, start_threshold: int) -> SeqRecord:
    """
    Processes a single contig by extracting its ends, running nucmer on them,
    and replacing the ends if significant alignments are found.
    """
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            chr_length = len(chr_rec.seq)
            
            # 1. Create a temporary FASTA file containing only the 5' and 3' end sequences
            ends_records = []
            
            seq_5_end = chr_rec.seq[:min(end_region, chr_length)]
            rec_5_end = SeqRecord(seq_5_end, id=f"{chr_rec.id}_5end", description="")
            ends_records.append(rec_5_end)
            
            if chr_length > end_region:
                seq_3_end = chr_rec.seq[-end_region:]
                rec_3_end = SeqRecord(seq_3_end, id=f"{chr_rec.id}_3end", description="")
                ends_records.append(rec_3_end)
            
            temp_ends_fa = os.path.join(temp_dir, f"ends_{chr_rec.id}.fasta")
            with open(temp_ends_fa, "w") as f:
                SeqIO.write(ends_records, f, "fasta")
            
            output_prefix = os.path.join(temp_dir, f"align_{chr_rec.id}")
            
            # 2. Run nucmer using the extracted end sequences as reference
            coords_file = run_nucmer(temp_ends_fa, consensus_fa, output_prefix, min_identity)
            
            # 3. Parse results from the optimized coords file
            best_5_hit, best_3_hit = parse_coords_for_ends_optimized(coords_file)

            new_seq_list = list(str(chr_rec.seq))
            
            # 5' End Replacement
            if best_5_hit and best_5_hit['qry_id'] in consensus_seqs:
                # 检查比对是否在 contig 前 start_threshold bp内
                if best_5_hit['ref_start'] > start_threshold:
                    logger.info(f"Found 5' end hit, but alignment start position ({best_5_hit['ref_start']}) is not within the first {start_threshold}bp. Skipping replacement.")
                    best_5_hit = None
                else:
                    chr_align_start_5 = best_5_hit['ref_start']
                    chr_align_end_5 = best_5_hit['ref_end']
                    
                    replace_seq_5 = consensus_seqs[best_5_hit['qry_id']][0:best_5_hit['qry_end']]
                    
                    original_len_5 = len("".join(new_seq_list))
                    new_seq_list = list(replace_seq_5 + "".join(new_seq_list[best_5_hit['ref_end']:]))
                    new_len_5 = len("".join(new_seq_list))
                    
                    logger.info(f"Replaced 5'end hit for '{chr_rec.id}' from [{chr_align_start_5}-{chr_align_end_5}] (orig_len:{original_len_5}, new_len:{new_len_5}) on consensus '{best_5_hit['qry_id']}' with {best_5_hit['identity']:.2f}% identity.")
            else:
                logger.info(f"No significant alignment found at the 5' end for '{chr_rec.id}'.")

            # 3' End Replacement
            if best_3_hit and best_3_hit['qry_id'] in consensus_seqs:
                chr_align_start_3_abs = best_3_hit['ref_start'] + (chr_length - end_region)

                # 检查比对是否在 contig 末尾 start_threshold bp内
                if chr_align_start_3_abs < chr_length - start_threshold:
                    logger.info(f"Found 3' end hit, but alignment start position ({chr_align_start_3_abs}) is not within the last {start_threshold}bp of contig length ({chr_length}). Skipping replacement.")
                    best_3_hit = None
                else:
                    chr_align_end_3_abs = best_3_hit['ref_end'] + (chr_length - end_region)
                    current_chr_len = len("".join(new_seq_list))
                    len_change = current_chr_len - chr_length
                    
                    adjusted_ref_start = chr_align_start_3_abs + len_change
                    
                    replace_seq_3 = consensus_seqs[best_3_hit['qry_id']][best_3_hit['qry_start']-1:]
                    
                    original_len_3 = len("".join(new_seq_list))
                    new_seq_list = list("".join(new_seq_list[:adjusted_ref_start-1]) + replace_seq_3)
                    new_len_3 = len("".join(new_seq_list))
                    
                    logger.info(f"Replaced 3'end hit for '{chr_rec.id}' from [{chr_align_start_3_abs}-{chr_align_end_3_abs}] (orig_len:{original_len_3}, new_len:{new_len_3}) on consensus '{best_3_hit['qry_id']}' with {best_3_hit['identity']:.2f}% identity.")
            else:
                logger.info(f"No significant alignment found at the 3' end for '{chr_rec.id}'.")
            
            final_seq = "".join(new_seq_list)
            return SeqRecord(Seq(final_seq), id=f"{chr_rec.id}_telomeres_replaced", description="")
    
    except Exception as e:
        logger.error(f"Error processing contig '{chr_rec.id}': {e}", exc_info=True)
        return chr_rec

def replace_telomeres_multi(original_fa: str, consensus_fa: str, output_fa: str, min_identity: float = 80, end_region: int = 20000, num_workers: int = 4, start_threshold: int = 100) -> None:
    """
    Main function now uses a ProcessPoolExecutor for parallel processing
    and optimizes nucmer calls by only aligning contig ends.
    """
    find_mummer_path()

    consensus_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(consensus_fa, "fasta")}
    logger.info(f"Loaded {len(consensus_seqs)} consensus telomere sequences.")
    
    all_contigs = list(SeqIO.parse(original_fa, "fasta"))
    total_contigs = len(all_contigs)
    logger.info(f"Total contigs to process: {total_contigs}")

    modified_records = []
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                process_single_contig,
                chr_rec,
                consensus_fa,
                consensus_seqs,
                min_identity,
                end_region,
                start_threshold
            ): chr_rec.id
            for chr_rec in all_contigs
        }

        for future in concurrent.futures.as_completed(futures):
            modified_records.append(future.result())

    with open(output_fa, "w") as f:
        SeqIO.write(modified_records, f, "fasta")
    
    logger.info(f"\nProcess complete. All {len(modified_records)} sequences saved to {output_fa}")
    
####################

if __name__ == "__main__":
    
    full_command = "python3 " + " ".join(sys.argv)
    logger.info(f"Command run: {full_command}")

    parser = argparse.ArgumentParser(
        description='Replace chromosome ends (telomeres) with a consensus sequence using MUMmer. Now with parallel processing and optimized end-only alignment.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-r', '--reference', required=True, help='FASTA file of the original chromosome contigs.')
    parser.add_argument('-c', '--consensus', required=True, help='FASTA file of the consensus telomere sequences.')
    parser.add_argument('-o', '--output', required=True, help='Path for the output FASTA file with replaced ends.')
    parser.add_argument('-i', '--min_identity', type=float, default=80.0, help='Minimum alignment identity percentage (default: 80.0).')
    parser.add_argument('-e', '--end_region', type=int, default=20000, help='Size of the region (in bp) at each end to scan for telomeres (default: 20000).')
    parser.add_argument('-p', '--processes', type=int, default=8, help='Number of parallel processes to use (default: 8).')
    parser.add_argument('-s', '--start_threshold', type=int, default=100, help='Maximum start position (in bp) on the contig ends for a valid alignment (default: 100).')
    
    args = parser.parse_args()
    
    try:
        replace_telomeres_multi(
            args.reference,
            args.consensus,
            args.output,
            min_identity=args.min_identity,
            end_region=args.end_region,
            num_workers=args.processes,
            start_threshold=args.start_threshold
        )
    except (ValueError, EnvironmentError, subprocess.CalledProcessError) as e:
        logger.critical(f"\nAn error occurred: {e}", exc_info=True)
    except Exception as e:
        logger.critical(f"\nAn unexpected error occurred: {e}", exc_info=True)
