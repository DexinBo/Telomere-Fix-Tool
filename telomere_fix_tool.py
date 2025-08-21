#!/usr/bin/env python3

import argparse
import logging
import logging.handlers
import os
import subprocess
import sys
import tempfile
import multiprocessing
import threading
from collections import namedtuple
from typing import List, Optional, Dict, Tuple
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --- 1. Logging Queue Configuration ---
def log_listener_process(queue: multiprocessing.Queue):
    log_format = '%(asctime)s - %(processName)s - %(levelname)s - \n%(message)s\n'
    root = logging.getLogger()
    handler = logging.StreamHandler()
    formatter = logging.Formatter(log_format, datefmt='%Y-%m-%d %H:%M:%S')
    handler.setFormatter(formatter)
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    while True:
        try:
            record = queue.get()
            if record is None: break
            logger = logging.getLogger(record.name)
            logger.handle(record)
        except Exception:
            import traceback
            print('Error in log listener:', file=sys.stderr)
            traceback.print_exc(file=sys.stderr)

def worker_configurer(queue: multiprocessing.Queue):
    root = logging.getLogger()
    if not root.handlers:
        handler = logging.handlers.QueueHandler(queue)
        root.addHandler(handler)
        root.setLevel(logging.INFO)

# --- 2. Core Processing Logic ---
Alignment = namedtuple('Alignment', ['s_start','s_end','q_start','q_end','s_len','q_len','identity','s_strand','q_strand','s_name','q_name'])

def run_command(command: list):
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result.stdout
    except FileNotFoundError:
        logging.error(f"Command '{command[0]}' not found.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command '{' '.join(command)}' failed: {e.stderr.strip()}")

def parse_show_coords(coords_output: str) -> List[Alignment]:
    alignments = []
    lines = coords_output.strip().split('\n')
    data_lines = [line for line in lines if not line.startswith('[') and not line.startswith('=')]
    if data_lines and data_lines[0].startswith('S1'): data_lines = data_lines[1:]
    for line in data_lines:
        if not line.strip(): continue
        parts = line.split('\t')
        try:
            s_start,s_end,q_start,q_end,s_len,q_len,identity,s_name,q_name = int(parts[0]),int(parts[1]),int(parts[2]),int(parts[3]),int(parts[4]),int(parts[5]),float(parts[6]),parts[9],parts[10]
            q_strand = '-' if int(parts[8]) < 0 else '+'
            alignments.append(Alignment(s_start,s_end,q_start,q_end,s_len,q_len,identity,'+',q_strand,s_name,q_name))
        except (ValueError, IndexError): continue
    return alignments

def find_best_end_alignment(alignments: List[Alignment], ref_len: int, end_type: str, min_identity: float, min_length: int, end_distance: int) -> Tuple[Optional[Alignment], Optional[Alignment]]:
    """
    Filters for the best end alignment.
    Returns: (The best alignment that fully meets the criteria, The best candidate that meets positional criteria but may fail others)
    """
    best_valid_aln, best_valid_score = None, (-1,-1)
    best_candidate_aln, best_candidate_score = None, (-1,-1)
    for aln in alignments:
        is_geometric_candidate = False
        if end_type == '5prime':
            if aln.s_start <= end_distance: is_geometric_candidate = True
        elif end_type == '3prime':
            if aln.s_end >= (ref_len - end_distance): is_geometric_candidate = True
        
        if is_geometric_candidate:
            current_score = (aln.identity, aln.s_len)
            if current_score > best_candidate_score:
                best_candidate_score, best_candidate_aln = current_score, aln
            
            is_valid = aln.identity >= min_identity and aln.s_len >= min_length
            if is_valid and current_score > best_valid_score:
                best_valid_score, best_valid_aln = current_score, aln
                
    return best_valid_aln, best_candidate_aln

def align_fragment(fragment_record: SeqRecord, consensus_file_path: str, params: argparse.Namespace) -> List[Alignment]:
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as tmp_frag_file:
        SeqIO.write(fragment_record, tmp_frag_file, "fasta")
        tmp_frag_path = tmp_frag_file.name
    
    prefix = tempfile.mktemp()
    delta_path = f"{prefix}.delta"
    filtered_delta_path = f"{prefix}.delta.filter"

    try:
        nucmer_command = ["nucmer", "--maxmatch", "--threads", str(params.nucmer_threads), "-p", prefix, tmp_frag_path, consensus_file_path]
        run_command(nucmer_command)
        
        delta_filter_command = ["delta-filter", "-i", "80", delta_path]
        filtered_delta_output = run_command(delta_filter_command)
        
        with open(filtered_delta_path, 'w') as f:
            f.write(filtered_delta_output)

        coords_output = run_command(["show-coords", "-THrd", filtered_delta_path])

    finally:
        os.remove(tmp_frag_path)
        if os.path.exists(delta_path):
            os.remove(delta_path)
        if os.path.exists(filtered_delta_path): 
            os.remove(filtered_delta_path)
    
    return parse_show_coords(coords_output)

def process_ref_sequence(
        ref_record: SeqRecord,
        consensus_records_dict: Dict[str, SeqRecord],
        consensus_file_path: str,
        params: argparse.Namespace
) -> SeqRecord:
    """Core processing function, generates a detailed multi-line report."""
    report_lines = []
    ref_len = len(ref_record.seq)
    margin = params.search_margin
    report_lines.append(f"--- Processing Report for: {ref_record.id} (Mode: {params.mode}) (Length: {ref_len} bp) ---")
    
    best_5p_aln, best_3p_aln = None, None

    # 1. Process 5' end
    five_prime_len = min(ref_len, margin)
    alignments_5p = align_fragment(ref_record[:five_prime_len], consensus_file_path, params)
    valid_5p_aln, candidate_5p_aln = find_best_end_alignment(alignments_5p, five_prime_len, '5prime', params.min_identity, params.min_length, params.end_distance)
    
    if valid_5p_aln:
        best_5p_aln = valid_5p_aln
        consensus_len = len(consensus_records_dict[best_5p_aln.q_name].seq)
        report_lines.append(
            f"  [5' End] ✓ Best match found: Consensus='{best_5p_aln.q_name}'(Len:{consensus_len}), Identity={best_5p_aln.identity:.2f}%, "
            f"RefPos=[{best_5p_aln.s_start}:{best_5p_aln.s_end}], ConsensusPos=[{best_5p_aln.q_start}:{best_5p_aln.q_end}], "
            f"AlnLength={best_5p_aln.s_len}")
    elif candidate_5p_aln:
        consensus_len = len(consensus_records_dict[candidate_5p_aln.q_name].seq)
        report_lines.append(
            f"  [5' End] ✗ No qualified match. Closest candidate: Consensus='{candidate_5p_aln.q_name}'(Len:{consensus_len}), Identity={candidate_5p_aln.identity:.2f}%, "
            f"RefPos=[{candidate_5p_aln.s_start}:{candidate_5p_aln.s_end}], ConsensusPos=[{candidate_5p_aln.q_start}:{candidate_5p_aln.q_end}], "
            f"AlnLength={candidate_5p_aln.s_len} (Thresholds: Length>={params.min_length}, Identity>={params.min_identity}%)")
    else:
        report_lines.append("  [5' End] ✗ No alignments found in the terminal region.")

    # 2. Process 3' end
    if ref_len > margin:
        three_prime_len = min(ref_len, margin)
        alignments_3p = align_fragment(ref_record[-three_prime_len:], consensus_file_path, params)
        valid_3p_aln_frag, candidate_3p_aln_frag = find_best_end_alignment(alignments_3p, three_prime_len, '3prime', params.min_identity, params.min_length, params.end_distance)
        
        if valid_3p_aln_frag:
            offset = ref_len - three_prime_len
            best_3p_aln = valid_3p_aln_frag._replace(s_start=valid_3p_aln_frag.s_start + offset, s_end=valid_3p_aln_frag.s_end + offset)
            consensus_len = len(consensus_records_dict[best_3p_aln.q_name].seq)
            report_lines.append(
                f"  [3' End] ✓ Best match found: Consensus='{best_3p_aln.q_name}'(Len:{consensus_len}), Identity={best_3p_aln.identity:.2f}%, "
                f"RefPos=[{best_3p_aln.s_start}:{best_3p_aln.s_end}], ConsensusPos=[{best_3p_aln.q_start}:{best_3p_aln.q_end}], "
                f"AlnLength={best_3p_aln.s_len}")
        elif candidate_3p_aln_frag:
            consensus_len = len(consensus_records_dict[candidate_3p_aln_frag.q_name].seq)
            report_lines.append(
                f"  [3' End] ✗ No qualified match. Closest candidate: Consensus='{candidate_3p_aln_frag.q_name}'(Len:{consensus_len}), Identity={candidate_3p_aln_frag.identity:.2f}%, "
                f"RefPos(frag):[{candidate_3p_aln_frag.s_start}:{candidate_3p_aln_frag.s_end}], ConsensusPos=[{candidate_3p_aln_frag.q_start}:{candidate_3p_aln_frag.q_end}], "
                f"AlnLength={candidate_3p_aln_frag.s_len} (Thresholds: Length>={params.min_length}, Identity>={params.min_identity}%)")
        else:
            report_lines.append("  [3' End] ✗ No alignments found in the terminal region.")

    # 3. Modify sequence and generate final report
    original_len = len(ref_record.seq)
    modified_description = []
    
    prefix, suffix = Seq(""), Seq("")
    middle_start_idx, middle_end_idx = 0, original_len
    
    if params.mode == 'replace':
        if best_5p_aln:
            consensus = consensus_records_dict[best_5p_aln.q_name].seq
            if best_5p_aln.q_strand == '+':
                # Forward alignment (Ref head vs Con head): slice the entire head of the Consensus from the start to the end of the alignment (q_end)
                prefix = consensus[:best_5p_aln.q_end]
            else:
                # Reverse alignment (Ref head vs Con tail): q_end is smaller; slice the entire Consensus tail including the alignment region (from q_end to the end) and use its reverse complement to replace the Ref's head
                tail_to_replace_with = consensus[best_5p_aln.q_end - 1:]
                prefix = tail_to_replace_with.reverse_complement()

            middle_start_idx = best_5p_aln.s_end
            modified_description.append("5p")

        if best_3p_aln:
            consensus = consensus_records_dict[best_3p_aln.q_name].seq
            if best_3p_aln.q_strand == '+':
                # Forward alignment (Ref tail vs Con tail): slice the entire tail of the Consensus from the start of the alignment (q_start) to the end
                suffix = consensus[best_3p_aln.q_start - 1:]
            else:
                # Reverse alignment (Ref tail vs Con head): q_start is larger; slice the entire Consensus head including the alignment region (from the start to q_start) and use its reverse complement to replace the Ref's tail
                head_to_replace_with = consensus[:best_3p_aln.q_start]
                suffix = head_to_replace_with.reverse_complement()

            middle_end_idx = best_3p_aln.s_start - 1
            if "3p" not in modified_description: modified_description.append("3p")

    elif params.mode == 'stitch':
        if best_5p_aln:
            consensus = consensus_records_dict[best_5p_aln.q_name].seq
            
            if best_5p_aln.q_strand == '+':
                # Forward alignment (Ref head vs Con head): slice the prefix of the Consensus
                prefix = consensus[:best_5p_aln.q_start - 1]
            else:
                # Reverse alignment (Ref head vs Con tail): q_start is larger; slice the tail of the original Consensus (from q_start onwards) and use its reverse complement as the prefix for the Ref
                tail_to_add = consensus[best_5p_aln.q_start:]
                prefix = tail_to_add.reverse_complement()
            
            middle_start_idx = best_5p_aln.s_start - 1
            modified_description.append("5p")

        if best_3p_aln:
            consensus = consensus_records_dict[best_3p_aln.q_name].seq
            
            if best_3p_aln.q_strand == '+':
                # Forward alignment (Ref tail vs Con tail): slice the suffix of the Consensus
                suffix = consensus[best_3p_aln.q_end:]
            else:
                # Reverse alignment (Ref tail vs Con head): q_end is smaller; slice the head of the original Consensus (up to q_end) and use its reverse complement as the suffix for the Ref
                head_to_add = consensus[:best_3p_aln.q_end - 1]
                suffix = head_to_add.reverse_complement()

            middle_end_idx = best_3p_aln.s_end
            if "3p" not in modified_description: modified_description.append("3p")

    if modified_description:
        middle = ref_record.seq[middle_start_idx:middle_end_idx]
        modified_seq = prefix + middle + suffix
        result_summary = f"  [Result] Sequence MODIFIED (mode: {params.mode}, ends: {','.join(modified_description)}). Original Length: {original_len}, New Length: {len(modified_seq)}"
        new_record = SeqRecord(modified_seq, id=f"{ref_record.id}_ends_modified", description=result_summary)
    else:
        result_summary = f"  [Result] Sequence NOT modified."
        new_record = ref_record

    report_lines.append(result_summary)
    logging.info("\n".join(report_lines))

    return new_record


def main():
    parser = argparse.ArgumentParser(
        description="Modify chromosome ends (telomeres) with a consensus sequence using MUMmer.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Mode options (--mode):
  replace: In 'replace' mode, the entire end of the reference sequence containing the alignment is removed 
           and replaced with the corresponding end of the consensus sequence.
  stitch:  In 'stitch' mode, the aligned region on the reference is kept as an anchor, and the consensus 
           sequence is used to fill in the missing parts outside of this anchor.
"""
    )
    parser.add_argument("--mode", choices=['replace', 'stitch'], default='stitch', help="Choose the logic for end modification. (default: stitch)")
    parser.add_argument("-r", "--ref", required=True, help="Path to the input reference sequence FASTA file.")
    parser.add_argument("-c", "--consensus", required=True, help="Path to the consensus sequence FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Path for the output FASTA file with modified sequences.")
    parser.add_argument("-P", "--processes", type=int, default=8, help="Number of parallel processes.")
    parser.add_argument("-T", "--nucmer_threads", type=int, default=8, help="Number of threads for each nucmer task.")
    parser.add_argument("--min_identity", type=float, default=80.0, help="Minimum identity percentage for a valid alignment.")
    parser.add_argument("--min_length", type=int, default=2000, help="Minimum length for a valid alignment.")
    parser.add_argument("--search_margin", type=int, default=20000, help="Length of the fragment from each end to search for alignments.")
    parser.add_argument("--end_distance", type=int, default=1000, help="Maximum distance from the sequence end for an alignment to be considered.")
    args = parser.parse_args()
    
    log_queue = multiprocessing.Queue(-1)
    listener_thread = threading.Thread(target=log_listener_process, args=(log_queue,))
    listener_thread.start()
    
    worker_configurer(log_queue)
    
    logging.info(f"Script running in '{args.mode}' mode.")
    total_threads, cpu_cores = args.processes * args.nucmer_threads, multiprocessing.cpu_count()
    logging.info(f"Configuration: {args.processes} parallel processes, with {args.nucmer_threads} threads per nucmer task.")
    if total_threads > cpu_cores:
        logging.warning(f"Total threads ({total_threads}) > CPU cores ({cpu_cores}). This may lead to performance degradation!")

    try:
        consensus_records = {rec.id: rec for rec in SeqIO.parse(args.consensus, "fasta")}
        if not consensus_records:
            logging.critical(f"Consensus sequence file '{args.consensus}' is empty.")
            sys.exit(1)
        logging.info(f"Successfully read {len(consensus_records)} consensus sequences.")
    except FileNotFoundError:
        logging.critical(f"Consensus sequence file '{args.consensus}' not found.")
        sys.exit(1)
        
    try:
        ref_records = list(SeqIO.parse(args.ref, "fasta"))
        if not ref_records:
            logging.warning(f"Reference sequence file '{args.ref}' is empty.")
            sys.exit(0)
        logging.info(f"Found {len(ref_records)} reference sequences to process.\n")
    except FileNotFoundError:
        logging.critical(f"Reference sequence file '{args.ref}' not found.")
        sys.exit(1)

    modified_records = []
    with ProcessPoolExecutor(max_workers=args.processes, initializer=worker_configurer, initargs=(log_queue,)) as executor:
        try:
            results_iterator = executor.map(process_ref_sequence, ref_records, repeat(consensus_records), repeat(args.consensus), repeat(args))
            modified_records = list(results_iterator)
        except Exception as e:
            logging.critical(f"A critical error occurred during processing: {e}")
            executor.shutdown(wait=True)
            log_queue.put(None)
            listener_thread.join()
            sys.exit(1)

    changed_count = sum(1 for rec in modified_records if rec.id.endswith('_ends_modified'))
    logging.info(f"\nProcessing complete. {changed_count} out of {len(ref_records)} reference sequences were modified.")
    
    SeqIO.write(modified_records, args.output, "fasta")
    logging.info(f"Results saved to: {args.output}")

    log_queue.put(None)
    listener_thread.join()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
