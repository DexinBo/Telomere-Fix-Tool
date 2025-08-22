# Patching Of Sequence Telomeres (POST)

**POST**: Telomere end correction tool leveraging consensus-guided alignment.

##  Features
- **Two core modes**: `stitch` (extend) and `replace` (correct)
- **Reverse alignment handling**: smart parsing of `show-coords` including reverse-complement logic
- **High-performance workflow**: fast filtering via `delta-filter`, refined by custom Python logic
- **Parallel processing & robust logging**: uses `multiprocessing`, includes logging-safe queues and cleanup
- **Temporary file management**: automatic cleanup of intermediate files

##  Installation
```bash
git clone https://github.com/DexinBo/POST.git
cd POST
```

## Workflow

The script follows this general workflow for each input reference sequence:

`Input Contig` -\> `Extract End Fragment` -\> `nucmer` -\> `delta-filter` -\> `show-coords` -\> `Custom Python Filtering` -\> `Apply Stitch/Replace Logic` -\> `Output Modified Contig`

---

##  Requirements
- Python 3.7+
- [Biopython](https://biopython.org/): `pip install biopython`
- **MUMmer v4** tools (`nucmer`, `delta-filter`, `show-coords`) must be installed and in your `PATH`

---

## Installation

1.  Clone this repository:
    ```bash
    git clone https://github.com/DexinBo/Telomere_fix_tool.git
    cd Telomere_fix_tool
    ```
2.  Ensure all software listed in the **Requirements** section is installed.

## Usage

The script is run from the command line.

### Basic Command

```bash
python your_script_name.py --mode <stitch|replace> -r <reference.fasta> -c <consensus.fasta> -o <output.fasta> [OPTIONS]
```

### Arguments

| Argument | Short | Description | Default |
| :--- | :--- | :--- | :--- |
| `--mode` | | The core logic to use: `stitch` (extend) or `replace` (correct). | `stitch` |
| `--ref` | `-r` | **Required.** Path to the input reference FASTA file (e.g., your contigs). | |
| `--consensus` | `-c` | **Required.** Path to the consensus FASTA file (e.g., plasmids, trusted references). | |
| `--output` | `-o` | **Required.** Path for the output FASTA file with modified sequences. | |
| `--processes` | `-P` | Number of parallel processes to run for processing multiple reference sequences. | `8` |
| `--nucmer_threads`| `-T` | Number of threads for each individual `nucmer` alignment task. | `8` |
| `--min_identity` | | The minimum percent identity for an alignment to be considered valid for modification. | `80.0` |
| `--min_length` | | The minimum length (in bp) for an alignment to be considered valid for modification. | `2000` |
| `--search_margin`| | The length of the fragment (in bp) from each end to use for searching alignments. | `20000` |
| `--end_distance` | | An alignment must start (for 5' end) or end (for 3' end) within this distance of the fragment's edge to be considered. | `1000` |

### Example Command

To extend the ends of `my_contigs.fasta` using `complete_plasmids.fasta` as the reference, using 8 parallel processes, with each alignment using 4 threads:

```bash
python post.py --mode stitch \
    -r ref.fasta \
    -c consensus.fasta \
    -o output.fasta \
    -P 8 \
    -T 4 \
    --min_identity 80.0 \
    --min_length 2000 \
    --search_margin 20000 \
    --end_distance 1000
```

## Core Logic: Stitch vs. Replace

The key difference between the two modes is how they treat the aligned region, or `[anchor]`.

### Stitch Mode (Extension)

The goal of stitch mode is to extend an incomplete sequence. It assumes the aligned region `[anchor]` is correct and uses it as a bridge to add the missing adjacent sequence from the consensus.

  - **The `[anchor]` is preserved on the reference sequence.**
  - The script adds the part of the consensus that lies **outside** the `[anchor]`.

**Diagram:**

```
Ref:       ------------[anchor]
Consensus: ------------[anchor][BBBBB]

Stitched Result: ------------[anchor][BBBBB]
```

### Replace Mode (Correction)

The goal of replace mode is to correct a potentially erroneous end. It assumes the entire end of the reference sequence, including the `[anchor]`, is untrustworthy and should be replaced by the "golden standard" end from the consensus.

  - **The `[anchor]` is discarded** along with the rest of the reference sequence's end.
  - The script replaces it with the entire corresponding end from the consensus, which **also includes the `[anchor]`**.

**Diagram:**

```
Ref:       ------------[anchor][XXXXX]   <-- This entire block is removed
Consensus: ------------[anchor][YYYYY]   <-- This entire block is copied

Replaced Result: ------------[anchor][YYYYY]
```
