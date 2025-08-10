#!/usr/bin/env python3
import os
import re
import subprocess
import tempfile
import shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import concurrent.futures # 导入并行处理库

# 其他函数 (find_mummer_path, run_nucmer, parse_coords_for_ends) 保持不变...

def find_mummer_path():
    """Find the path to nucmer, delta-filter, and show-coords."""
    for tool in ['nucmer', 'delta-filter', 'show-coords']:
        if not shutil.which(tool):
            raise EnvironmentError(
                f"'{tool}' not found in PATH. Please ensure MUMmer is installed and accessible."
            )

def run_nucmer(consensus_fa, chromosome_fa, output_prefix, min_identity=80, min_length=1000):
    """
    Runs nucmer alignment and filters the results.
    The consensus_fa and chromosome_fa can now contain multiple sequences.
    """
    delta_file = f"{output_prefix}.delta"
    cmd_nucmer = f"nucmer --maxmatch --prefix={output_prefix} {chromosome_fa} {consensus_fa}"
    subprocess.run(cmd_nucmer, shell=True, check=True, capture_output=True)

    filtered_delta_file = f"{output_prefix}.filtered.delta"
    cmd_filter = f"delta-filter -i {min_identity} -l {min_length} {delta_file} > {filtered_delta_file}"
    subprocess.run(cmd_filter, shell=True, check=True, capture_output=True)

    coords_file = f"{output_prefix}.coords"
    cmd_coords = f"show-coords -rclT {filtered_delta_file} > {coords_file}"
    subprocess.run(cmd_coords, shell=True, check=True, capture_output=True)

    return coords_file
#
def parse_coords_for_ends(coords_file, chr_length, end_region_size=20000):
    """
    Parses alignment results to find the best hit for both 5' and 3' ends
    based on a combined score of identity and alignment length.
    """
    best_5_end_hit = None
    best_3_end_hit = None
    
    # 初始化最佳得分
    best_5_end_score = -1.0
    best_3_end_score = -1.0

    with open(coords_file) as f:
        lines = [line for line in f if not line.startswith('/') and not line.startswith('[')]
        
    for line in lines:
        if not line.strip():
            continue
            
        try:
            cols = line.strip().split('\t')
            if len(cols) < 13: continue

            ref_start, ref_end = int(cols[0]), int(cols[1])
            qry_start, qry_end = int(cols[2]), int(cols[3])
            identity = float(cols[6])
            
            # 计算比对长度
            alignment_length = ref_end - ref_start + 1
            
            # 计算综合得分
            current_score = identity * alignment_length
            
            current_hit = {
                "ref_start": ref_start, "ref_end": ref_end,
                "qry_start": qry_start, "qry_end": qry_end,
                "identity": identity,
                "qry_id": cols[12]
            }

        except (IndexError, ValueError):
            continue

        # 检查是否为5'端最佳比对
        if ref_start <= end_region_size:
            if current_score > best_5_end_score:
                best_5_end_score = current_score
                best_5_end_hit = current_hit

        # 检查是否为3'端最佳比对
        if ref_end >= (chr_length - end_region_size):
            if current_score > best_3_end_score:
                best_3_end_score = current_score
                best_3_end_hit = current_hit

    return best_5_end_hit, best_3_end_hit
#

# 新增的函数，用于封装单个contig的处理逻辑
def process_single_contig(chr_rec, consensus_fa, consensus_seqs, min_identity, end_region):
    """Processes a single contig and returns the modified SeqRecord."""
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_chr_fa = os.path.join(temp_dir, f"chr_{chr_rec.id}.fasta")
            with open(temp_chr_fa, "w") as f:
                SeqIO.write(chr_rec, f, "fasta")
            
            output_prefix = os.path.join(temp_dir, f"align_{chr_rec.id}")
            
            coords_file = run_nucmer(consensus_fa, temp_chr_fa, output_prefix, min_identity)
            
            chr_length = len(chr_rec.seq)
            best_5_hit, best_3_hit = parse_coords_for_ends(coords_file, chr_length, end_region)

            new_seq_list = list(str(chr_rec.seq))
            
            # 5' End Replacement
            if best_5_hit and best_5_hit['qry_id'] in consensus_seqs:
                print(f"  Found 5' end hit for '{chr_rec.id}' on consensus '{best_5_hit['qry_id']}' with {best_5_hit['identity']:.2f}% identity.")
                replace_seq_5 = consensus_seqs[best_5_hit['qry_id']][0:best_5_hit['qry_end']]
                new_seq_list = list(replace_seq_5 + "".join(new_seq_list[best_5_hit['ref_end']:]))
                print("  Replaced 5' end.")
            else:
                print(f"  No significant alignment found at the 5' end for '{chr_rec.id}'.")

            # 3' End Replacement
            if best_3_hit and best_3_hit['qry_id'] in consensus_seqs:
                current_chr_len = len("".join(new_seq_list))
                len_change = current_chr_len - chr_length
                adjusted_ref_start = best_3_hit['ref_start'] + len_change
                
                print(f"  Found 3' end hit for '{chr_rec.id}' on consensus '{best_3_hit['qry_id']}' with {best_3_hit['identity']:.2f}% identity.")
                replace_seq_3 = consensus_seqs[best_3_hit['qry_id']][best_3_hit['qry_start']-1:]
                new_seq_list = list("".join(new_seq_list[:adjusted_ref_start-1]) + replace_seq_3)
                print("  Replaced 3' end.")
            else:
                print(f"  No significant alignment found at the 3' end for '{chr_rec.id}'.")
            
            final_seq = "".join(new_seq_list)
            return SeqRecord(Seq(final_seq), id=f"{chr_rec.id}_telomeres_replaced", description="")
    
    except Exception as e:
        print(f"Error processing contig '{chr_rec.id}': {e}")
        return chr_rec # Return original record on failure


def replace_telomeres_multi(original_fa, consensus_fa, output_fa, min_identity=80, end_region=20000, num_workers=4):
    """
    Main function now uses a ProcessPoolExecutor for parallel processing.
    """
    find_mummer_path()

    consensus_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(consensus_fa, "fasta")}
    
    # 将所有contig记录读取到内存
    all_contigs = list(SeqIO.parse(original_fa, "fasta"))
    total_contigs = len(all_contigs)
    print(f"Total contigs to process: {total_contigs}")

    modified_records = []
    
    # 使用 ProcessPoolExecutor 进行并行处理
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
        # 提交所有contig的任务
        futures = [
            executor.submit(
                process_single_contig,
                chr_rec,
                consensus_fa,
                consensus_seqs,
                min_identity,
                end_region
            )
            for chr_rec in all_contigs
        ]

        # 收集结果
        for future in concurrent.futures.as_completed(futures):
            modified_records.append(future.result())
            print(f"Processed {len(modified_records)}/{total_contigs} contigs.")
    
    # 写入最终输出
    with open(output_fa, "w") as f:
        SeqIO.write(modified_records, f, "fasta")
    
    print(f"\nProcess complete. All {len(modified_records)} sequences saved to {output_fa}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Replace chromosome ends (telomeres) with a consensus sequence using MUMmer. Now with parallel processing.',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-r', '--reference', required=True, help='FASTA file of the original chromosome contigs.')
    parser.add_argument('-c', '--consensus', required=True, help='FASTA file of the consensus telomere sequences.')
    parser.add_argument('-o', '--output', required=True, help='Path for the output FASTA file with replaced ends.')
    parser.add_argument('-i', '--min_identity', type=float, default=80.0, help='Minimum alignment identity percentage (default: 80.0).')
    parser.add_argument('-e', '--end_region', type=int, default=20000, help='Size of the region (in bp) at each end to scan for telomeres (default: 20000).')
    parser.add_argument('-p', '--processes', type=int, default=8, help='Number of parallel processes to use (default: 8).')
    
    args = parser.parse_args()
    
    try:
        replace_telomeres_multi(
            args.reference,
            args.consensus,
            args.output,
            min_identity=args.min_identity,
            end_region=args.end_region,
            num_workers=args.processes
        )
    except (ValueError, EnvironmentError, subprocess.CalledProcessError) as e:
        print(f"\nAn error occurred: {e}")
