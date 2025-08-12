This script is a comprehensive and optimized tool designed to identify and replace telomeric or repetitive regions at the ends of multiple contig sequences within a FASTA file. It leverages the powerful MUMmer suite for sequence alignment and Biopython for sequence manipulation, all while employing parallel processing for efficiency.

Key Features

1. Parallel Processing: It uses Python's concurrent.futures module to process multiple contigs simultaneously, dramatically speeding up execution on multi-core systems.

2. Intelligent Alignment Selection: Instead of just using percentage identity, it now employs a combined scoring strategy (identity multiplied by alignment length) to select the "best" alignment for replacement, ensuring more biologically meaningful results.

3. Custom Replacement Logic: It precisely implements your desired replacement strategy: replacing a contig's end with the entire relevant portion of the consensus sequence (from its start to alignment end for the 5' side, or from alignment start to its end for the 3' side).

4. Robust Dependency Checking: It verifies the presence of required external MUMmer tools before starting, preventing runtime errors.

5. Clean Operation: Uses temporary directories to manage intermediate alignment files, automatically cleaning them up upon completion.

6. Multi-Sequence Support: Handles FASTA files containing multiple contigs (reference) and multiple consensus sequences (query), iterating through each.


## Usage
```
python3 replace_telomeres.py [-h] -r REFERENCE -c CONSENSUS -o OUTPUT [-m {replace,splice}] [-i MIN_IDENTITY] [-l MIN_LENGTH]
                             [-e END_REGION] [-p PROCESSES] [-s START_THRESHOLD]

Modify chromosome ends (telomeres) with a consensus sequence using MUMmer. Now with parallel processing and two modification modes.

options:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        FASTA file of the original chromosome contigs.
  -c CONSENSUS, --consensus CONSENSUS
                        FASTA file of the consensus telomere sequences.
  -o OUTPUT, --output OUTPUT
                        Path for the output FASTA file with modified ends.
  -m {replace,splice}, --mode {replace,splice}
                        Modification mode: 
                        "replace": Replace the aligned contig end with the corresponding consensus sequence and its flanking region. (Default)
                        "splice": Prepend/append only the flanking region of the consensus sequence without replacing the original contig sequence.
  -i MIN_IDENTITY, --min_identity MIN_IDENTITY
                        Minimum alignment identity percentage (default: 80.0).
  -l MIN_LENGTH, --min_length MIN_LENGTH
                        Minimum alignment length (default: 2000).
  -e END_REGION, --end_region END_REGION
                        Size of the region (in bp) at each end to scan for telomeres (default: 20000).
  -p PROCESSES, --processes PROCESSES
                        Number of parallel processes to use (default: 8).
  -s START_THRESHOLD, --start_threshold START_THRESHOLD
                        Maximum start position (in bp) on the contig ends for a valid alignment (default: 1000).```
