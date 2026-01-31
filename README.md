# A Targeted Nanopore Sequencing Pipeline for Tandem Gene Array Engineering

## Overview

This repository provides a shell-based analysis pipeline for targeted Oxford Nanopore sequencing data, specifically designed for the analysis of **tandem gene arrays generated or modified by gene duplication-based approaches**.

The pipeline enables barcode-wise processing to quantify and characterize copy number variation and structural features within tandemly duplicated gene arrays. It performs read alignment, target and flanking sequence extraction, coverage calculation, and optional normalization using a `.fsa` reference genome.

The main workflow is implemented in `Alignment.sh` and is designed to be reproducible, modular, and suitable for publication and long-term archiving (e.g. GitHub / Zenodo).

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

### Python scripts (included in `./NSA/`)

* `Target_seq_extraction.py`
* `Extract_tgt_only_mappy.py`
* `Extract_tgt_flanks_mappy.py`
* `Bedgraph_normalize.py`
* `yass_yop_to_svg.py`

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

If you use this pipeline in your research, please cite the Zenodo archive (DOI will be added after release).

In addition, this pipeline makes use of third-party software for coverage normalization (see below), which should be cited appropriately.

A `CITATION.cff` file is recommended for GitHub-based citation support.

---

## Third-party tools and attribution

### Bedgraph normalization

Coverage normalization is performed using `Bedgraph_normalize.py`, which is part of the **Bedgraph_norm_ratio** project developed by poccopen:

* GitHub repository: [https://github.com/poccopen/Bedgraph_norm_ratio](https://github.com/poccopen/Bedgraph_norm_ratio)

Please cite the original repository (and any associated publication, if applicable) when using normalized bedgraph outputs generated by this pipeline.

---

## License

MIT.

---

## Contact

For questions or issues, please use the GitHub Issues page after publication.
