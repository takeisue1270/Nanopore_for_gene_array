# A Targeted Nanopore Sequencing Pipeline for Tandem Gene Array Engineering

*[日本語版 README はこちら](README.ja.md)*

## Overview

This repository provides a shell-based analysis pipeline for targeted Oxford Nanopore sequencing data, specifically designed for the analysis of **tandem gene arrays generated or modified by gene duplication-based approaches**.

The pipeline enables barcode-wise processing to quantify and characterize copy number variation and structural features within tandemly duplicated gene arrays. It performs read alignment, target and flanking sequence extraction, coverage calculation, and optional normalization using a `.fsa` reference genome.

The main workflow is implemented in `Alignment.sh`. Every step is controlled by command-line arguments or environment variables, so a run can be recorded and repeated exactly.

The normalization step (step 5) works **only for *Saccharomyces cerevisiae***: it reads chromosome lengths from a reference whose sequence names are of the `chrI` form. Every other step is organism-independent.

This is a research tool, provided as is and without warranty of any kind, as the [licence](LICENSE) states. It is not intended for diagnostic, therapeutic or clinical use. Work with engineered organisms is regulated — in Japan by the Cartagena Act — and your institution's biosafety committee, not this repository, decides what you may do.

---

## Pipeline Summary

The pipeline follows a step-wise workflow designed for the analysis of tandem gene arrays engineered by gene duplication. For each barcode, the following steps are performed:

1. Merge FASTQ files per barcode with minimum read length filtering
2. Align reads to a reference genome using **minimap2** and generate sorted BAM files
3. Generate WIG coverage tracks (optional)
4. Generate bedgraph coverage tracks (optional)
5. Normalize bedgraph coverage using a `.fsa` reference genome (optional)
6. Extract reads spanning the entire targeted tandem gene array (optional)
7. Extract target-hit reads only (optional)
8. Extract upstream and downstream flanking sequences around target hits (optional)
9. Align extracted sequences and generate YASS SVG visualizations (optional)
10. Perform BLASTN analysis for sequence validation or classification (optional)

Each step can be enabled or disabled independently using environment variables.

---

## Requirements

### Core tools

* minimap2
* samtools
* bedtools
* seqkit
* python >= 3.8

### Optional tools

* igvtools (for WIG generation)
* yass (for visualization)
* BLAST+ (for BLASTN analysis)

### Python packages

```bash
pip install mappy biopython numpy
```

* `mappy` — minimap2's Python bindings, used by the three extraction scripts
* `biopython` and `numpy` — used by `Bedgraph_normalize.py` only

### Python scripts (included in `./NSA/`)

* `Target_seq_extraction.py`
* `Extract_tgt_only_mappy.py`
* `Extract_tgt_flanks_mappy.py`
* `Bedgraph_normalize.py` — third-party; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
* `yass_yop_to_svg.py`

---

## Reference files you have to supply

The pipeline expects these in `./NSA/`. None of them is included — they depend on your target and your assembly, and sequence files of that size do not belong in a Git repository.

| File | What it is | How to make it |
| --- | --- | --- |
| `${REFBASE}.fa` | The alignment reference. For an engineered array this is usually the construct or the array unit, not a whole genome. | Export from your sequence editor as FASTA. The minimap2 index `.mmi` is built automatically on first use. |
| `${REFBASE}.UP1000.fa` | 1 kb immediately upstream of the array in the target locus. Used to find reads that span the whole array. | Cut from the locus sequence. The name must match `REFBASE`. |
| `${REFBASE}.DOWN1000.fa` | 1 kb immediately downstream, same purpose. | As above. |
| `ACT1.fa` | A single-copy control locus, used by `RUN_BLAST=1` to give a per-sample denominator. | The *ACT1* coding sequence of your strain; any single-copy locus serves the same purpose. |
| `S288C_reference_sequence_R64-2-1_20150113_0CUP1RU_1rDNARU_phiX.fsa` | The *S. cerevisiae* S288C R64-2-1 reference, modified in-house: the CUP1 repeat reduced to 0 units, the rDNA repeat to 1 unit, with phiX appended. Used by `RUN_FSA_NORM=1` and by the flanking-sequence mapping. | Start from the [SGD R64-2-1 genome release](https://www.yeastgenome.org/) and make the same edits for your own repeat units, or point `FSA_REF` at your own reference. Chromosome names must stay in `chrI` form for normalization to work. |

Collapsing a tandem repeat in the normalization reference is deliberate: leaving the repeat at its reference copy number makes the coverage ratio of an engineered array uninterpretable. Document whatever you do here, because the normalized output cannot be read without it.

---

## Directory Structure

```text
project_root/
├── Alignment.sh
├── NSA/
│   ├── Target_seq_extraction.py
│   ├── Extract_tgt_only_mappy.py
│   ├── Extract_tgt_flanks_mappy.py
│   ├── Bedgraph_normalize.py
│   └── yass_yop_to_svg.py
├── fastq_pass/
│   └── barcodeXX/
│       └── *.fastq(.gz)
└── README.md
```

---

## Usage

### Basic command

```bash
./Alignment.sh PREFIX START END REFBASE [THREADS] [MINLEN]
```

### Example

```bash
./Alignment.sh T975 1 96 S288C 14 200
```

### Arguments

| Argument | Description                                                   |
| -------- | ------------------------------------------------------------- |
| PREFIX   | Sample prefix (e.g. T975)                                     |
| START    | First barcode index (e.g. 1)                                  |
| END      | Last barcode index (e.g. 96)                                  |
| REFBASE  | Reference base name (FASTA expected at `./NSA/${REFBASE}.fa`) |
| THREADS  | Number of threads (default: 14)                               |
| MINLEN   | Minimum read length filter (default: 200)                     |

---

## Environment Variables

Key environment variables controlling optional steps:

| Variable     | Description                                     |
| ------------ | ----------------------------------------------- |
| MAKE_WIG     | Generate WIG file (0/1)                         |
| MAKE_BED     | Generate bedgraph (0/1)                         |
| RUN_FSA_NORM | Normalize bedgraph using `.fsa` reference (0/1) |
| RUN_TARGET   | Run target sequence extraction (0/1)            |
| RUN_TGT_ONLY | Extract target-hit reads only (0/1)             |
| RUN_FLANKS   | Extract flanking sequences (0/1)                |
| RUN_YASS     | Generate YASS SVGs (0/1)                        |
| RUN_BLAST    | Run BLASTN analysis (0/1)                       |

---

## Output Files

Typical output files include:

* `merged.<PREFIX>.<BARCODE>.fastq`
* `<PREFIX>.<BARCODE>.exp.sort.bam`
* `<PREFIX>.<BARCODE>.exp.<REFBASE>.bedgraph`
* `<PREFIX>.<BARCODE>.exp.<REFBASE>.norm.bedgraph` (if normalization enabled)
* `merged.<PREFIX>.<BARCODE>.tgt_flanks.fa` (if flanking extraction enabled)

All files are written to the current working directory.

---

## Notes on Reproducibility

* All parameters are explicitly controlled via command-line arguments or environment variables.
* Existing output files with identical names will be overwritten.
* FASTA headers are sanitized when necessary to avoid SAM parsing issues.

---

## Citation

Cite the Zenodo archive of this pipeline:

* Takesue H. A Targeted Nanopore Sequencing Pipeline for Tandem Gene Array Engineering. Zenodo. doi:[10.5281/zenodo.18440611](https://doi.org/10.5281/zenodo.18440611)

The metadata is in [CITATION.cff](CITATION.cff), so GitHub's **Cite this repository** button gives the same thing.

Cite the tools the pipeline runs as well — your results are theirs as much as this wrapper's:

* minimap2: Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191
* samtools / htslib: Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience*. 2021;10(2):giab008. doi:10.1093/gigascience/giab008
* bedtools: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010;26(6):841-842. doi:10.1093/bioinformatics/btq033
* SeqKit: Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*. 2016;11(10):e0163962. doi:10.1371/journal.pone.0163962
* BLAST+ (if `RUN_BLAST=1`): Camacho C, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421. doi:10.1186/1471-2105-10-421
* YASS (if `RUN_YASS=1`): Noé L, Kucherov G. YASS: enhancing the sensitivity of DNA similarity search. *Nucleic Acids Research*. 2005;33(Web Server issue):W540-W543. doi:10.1093/nar/gki478
* IGVtools (if `MAKE_WIG=1`): Robinson JT, et al. Integrative genomics viewer. *Nature Biotechnology*. 2011;29(1):24-26. doi:10.1038/nbt.1754

If you used `RUN_FSA_NORM=1`, cite the normalization script's own archive as well — it is other people's work:

* Okada S. poccopen/Bedgraph_norm_ratio. Zenodo. doi:[10.5281/zenodo.11515695](https://doi.org/10.5281/zenodo.11515695)

---

## Third-party code

`NSA/Bedgraph_normalize.py` is redistributed from the **Bedgraph_norm_ratio** project by Satoshi Okada ([@poccopen](https://github.com/poccopen)) under CC BY 4.0, and is **not** covered by this repository's MIT licence. Details, including what the licence requires of you, are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

---

## License

The code in this repository is licensed under the MIT License; see [LICENSE](LICENSE). The one exception is `NSA/Bedgraph_normalize.py`, covered by CC BY 4.0 as described above.

Oxford Nanopore Technologies, ONT, MinION, GridION and PromethION are trademarks of Oxford Nanopore Technologies plc. The other tools named here belong to their respective authors. This project is independent of all of them: not affiliated with, endorsed by, or supported by any of them, and it names them only to say what it reads and what it runs.

---

## Contact

Questions and bug reports: please open a GitHub issue.
