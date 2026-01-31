#!/usr/bin/env bash
# Description:
#   This script performs a barcode-wise analysis of Nanopore sequencing data,
#   including read alignment, target/flanking sequence extraction, coverage
#   calculation, and optional normalization using a .fsa reference genome.
#
# Fast Alignment Pipeline: Minimap2 → BAM → Target_seq_extraction → (optional) tgt-only extract → (optional) YASS SVG → (optional) BLASTN
#
# Usage:
#   ./Alignment.sh PREFIX START END REFBASE [THREADS] [MINLEN]
# Example:
#   ./Alignment.sh T975 1 96 NxAlphoid 14 200
#
# Notes:
#   - REF is fixed: ./NSA/${REFBASE}.fa  (mmi auto-build/reuse)
#   - TARGET is also REFBASE: ./NSA/${REFBASE}.UP1000.fa ./NSA/${REFBASE}.DOWN1000.fa ./NSA/${REFBASE}.fa
#   - YASS output dir: <positive_dir>/${PREFIX}_${barcode}/
#
# Dependencies:
#   minimap2, samtools, seqkit, (igvtools optional), (bedtools optional), yass(optional), python3,
#   ./NSA/Target_seq_extraction.py, ./NSA/yass_yop_to_svg.py, ./NSA/Extract_tgt_only_mappy.py(optional)
#
# All output files are generated in the current working directory.
# Existing files with the same names will be overwritten.

set -euo pipefail
export LC_ALL=C

# ===== Configurable parameters (can be overridden by environment variables) =====
MAKE_WIG=${MAKE_WIG:-0}           # 1: run igvtools using BAM file
MAKE_BED=${MAKE_BED:-0}           # 1: run bedtools using BAM file
RUN_FSA_NORM=${RUN_FSA_NORM:-0}   # 1: normalize bedgraph using reference (Bedgraph_normalize.py)
RUN_TARGET=${RUN_TARGET:-0}       # 1: run Target_seq_extraction
RUN_TGT_ONLY=${RUN_TGT_ONLY:-0}   # 1: extract target-only reads from merged.fastq (mappy)
RUN_FLANKS=${RUN_FLANKS:-0}       # 1: extract upstream/downstream flanking sequences (mappy)
RUN_YASS=${RUN_YASS:-0}           # 1: generate YASS SVG (requires RUN_TARGET=1)
RUN_BLAST=${RUN_BLAST:-0}         # 1: run BLAST (optional)

MM2_CHUNK=${MM2_CHUNK:-2g}           # minimap2 -K
SORT_MEM=${SORT_MEM:-1G}             # samtools sort -m
TMPDIR=${TMPDIR:-/tmp}
NORM_PY=${NORM_PY:-./NSA/Bedgraph_normalize.py}

# ===== Command-line arguments =====
PREFIX=${1:?}
START=${2:?}
END=${3:?}
REFBASE=${4:?}
THREADS=${5:-14}
MINLEN=${6:-200}

# ===== Parameters for YASS =====
YASS_PREFIX="${PREFIX}"  
TARGET="${REFBASE}"             

# ===== Reference files =====
REF="./NSA/${REFBASE}.fa"
REF_MMI="${REF}.mmi"

UP="./NSA/${TARGET}.UP1000.fa"
DOWN="./NSA/${TARGET}.DOWN1000.fa"
TARGETFA="./NSA/${TARGET}.fa"

# ===== Mapping target-derived sequences to the .fsa reference =====
FSA_REF="./NSA/S288C_reference_sequence_R64-2-1_20150113_0CUP1RU_1rDNARU_phiX.fsa"
FSA_MMI="${FSA_REF}.mmi"

# ===== Parametes for BLAST =====
BLAST_TASK=${BLAST_TASK:-megablast}
BLAST_EVALUE=${BLAST_EVALUE:-1e-50}
BLAST_THREADS=${BLAST_THREADS:-$THREADS}
BLAST_MAX_TARGET_SEQS=${BLAST_MAX_TARGET_SEQS:-1}
BLAST_OUTFMT=${BLAST_OUTFMT:-6}
BLAST_SUMMARY=${BLAST_SUMMARY:-"blast_summary.${PREFIX}.csv"}
ACT1="./NSA/ACT1.fa"
DB_REF_PREFIX="./NSA/${REFBASE}.fa.db"
DB_ACT_PREFIX="./NSA/ACT1.fa.db"

# ===== Pre-processing =====
echo "[*]Decompressing input FASTQ files (if needed)..."
if command -v pigz >/dev/null 2>&1; then
  find ./fastq_pass -type f -name "*.gz" -print0 | xargs -0 -P "${THREADS}" -I{} pigz -d -p "${THREADS}" -f "{}" || true
else
  gzip -d -r ./fastq_pass/ 2>/dev/null || true
fi

# =====  minimap2 index ===== 
if [[ ! -f "${REF_MMI}" ]]; then
  echo "[*] Building minimap2 index: ${REF_MMI}"
  minimap2 -x map-ont -d "${REF_MMI}" "${REF}"
fi

# ===== minimap2 index for fsa ===== 
if [[ ! -f "${FSA_MMI}" ]]; then
  echo "[*] Building minimap2 index for FSA..."
  minimap2 -x map-ont -d "${FSA_MMI}" "${FSA_REF}"
fi

# ===== BLAST database preparation (optional) ===== 
if [[ "${RUN_BLAST}" -eq 1 ]]; then
  if [[ ! -f "${DB_REF_PREFIX}.nsq" ]]; then
    echo "[*] makeblastdb (REFBASE) ..."
    makeblastdb -in "${REF}" -dbtype nucl -parse_seqids -out "${DB_REF_PREFIX}"
  fi
  if [[ ! -f "${DB_ACT_PREFIX}.nsq" ]]; then
    echo "[*] makeblastdb (ACT1) ..."
    makeblastdb -in "${ACT1}" -dbtype nucl -parse_seqids -out "${DB_ACT_PREFIX}"
  fi
  [[ -f "${BLAST_SUMMARY}" ]] || echo "Sample,TotalReads,REF_HitReads,ACT1_HitReads,REF/Total,ACT1/Total,REF/ACT1" > "${BLAST_SUMMARY}"
fi

# ===== Main loop over barcodes =====
for i in $(seq -f "%02g" "${START}" "${END}"); do
  echo "[*] Processing barcode ${PREFIX}.${i}..."
  BARCODE="${i}"
  MERGED_FQ="merged.${PREFIX}.${BARCODE}.fastq"
  OUT_BAM="${PREFIX}.${BARCODE}.exp.sort.bam"
  OUT_WIG="${PREFIX}.${BARCODE}.exp.sort.wig"
  OUT_BED="${PREFIX}.${BARCODE}.exp.${REFBASE}.bedgraph"

  # 1) Generate merged.fastq if not present
  if [[ ! -f "${MERGED_FQ}" ]]; then
    seqkit seq -m "${MINLEN}" ./fastq_pass/barcode${BARCODE}/*.fastq -o "${MERGED_FQ}"
  fi

  # 2) Map reads with minimap2 and generate sorted BAM
  minimap2 -a -x map-ont -t "${THREADS}" -K "${MM2_CHUNK}" -s "${MINLEN}" --secondary=no "${REF_MMI}" "${MERGED_FQ}" \
    | samtools view -@ "${THREADS}" -b \
    | samtools sort -@ "${THREADS}" -m "${SORT_MEM}" -o "${OUT_BAM}"

  # 3) Generate WIG file (optional)
  if [[ "${MAKE_WIG}" -eq 1 ]]; then
    igvtools count -w 1 --bases "${OUT_BAM}" "${OUT_WIG}" "${REF}"
  fi

  # 4) Generate bedgraph (optional)
  if [[ "${MAKE_BED}" -eq 1 ]]; then
    bedtools genomecov -ibam "${OUT_BAM}" -bg > "${OUT_BED}"
 
  # 5) Normalize bedgraph using .fsa reference (optional)
  # Normalization is performed using genome length information from the .fsa reference.
    if [[ "${RUN_FSA_NORM}" -eq 1 ]]; then
      if [[ -f "${NORM_PY}" ]]; then
        OUT_BED_NORM="${OUT_BED%.bedgraph}.norm.bedgraph"
        rm -f "${OUT_BED_NORM}"  
        python3 "${NORM_PY}" "${FSA_REF}" "${OUT_BED}" "${OUT_BED_NORM}"
      else
        echo "[WARN] NORM_PY not found: ${NORM_PY} (skip normalization)"
      fi
    fi
  fi

  # 6) Target whole tandem gene array_extraction (optional)
  POSITIVE_OUT=""
  if [[ "${RUN_TARGET}" -eq 1 ]]; then
    python3 ./NSA/Target_seq_extraction.py "${UP}" "${DOWN}" "${TARGETFA}" "${MERGED_FQ}"
    POSITIVE_IN="${MERGED_FQ%.fastq}.target.positive"
    POSITIVE_OUT="RCC.${PREFIX}.${BARCODE}.${REFBASE}.positive"
    if [[ -f "${POSITIVE_IN}" ]]; then
      mv -f "${POSITIVE_IN}" "${POSITIVE_OUT}"
    else
      echo "[WARN] positive file not found: ${POSITIVE_IN} (skip downstream for ${PREFIX}.${BARCODE})"
      POSITIVE_OUT=""
    fi
  fi

  # 7) Extract target hits only (optional)
  if [[ "${RUN_TGT_ONLY}" -eq 1 ]]; then
    python3 ./NSA/Extract_tgt_only_mappy.py "${TARGETFA}" "${MERGED_FQ}" --threads "${THREADS}" --min-mapq 10
    # output: merged.<PREFIX>.<BARCODE>.tgt.fasta
  fi

  # 8) Extract flanking sequences around target hits (optional)
  if [[ "${RUN_FLANKS}" -eq 1 ]]; then
    python3 ./NSA/Extract_tgt_flanks_mappy.py "${TARGETFA}" "${MERGED_FQ}" --threads "${THREADS}" --min-mapq 10
    # Output: merged.<PREFIX>.<BARCODE>.tgt_flanks.fa

    FA="merged.${PREFIX}.${BARCODE}.tgt_flanks.fa"

    if [[ -f "${FA}" ]]; then
      # Remove leading '@' from FASTA headers to avoid misinterpretation of read names as SAM header lines by samtools
      FA_CLEAN="${FA%.fa}.noat.fa"
      sed -E '/^>/ s/^>@/>/' "${FA}" > "${FA_CLEAN}"
    
      BASENAME=$(basename "${FA_CLEAN}" .fa)
      BAM="${BASENAME}.vs_fsa.sort.bam"
      BEDGRAPH="${BASENAME}.vs_fsa.bedgraph"
    
      minimap2 -a -x map-ont -t "${THREADS}" "${FSA_MMI}" "${FA_CLEAN}" \
        | samtools view -@ "${THREADS}" -b \
        | samtools sort -@ "${THREADS}" -o "${BAM}"
    
      samtools index "${BAM}"
      bedtools genomecov -ibam "${BAM}" -bg > "${BEDGRAPH}"
    else
      echo "[WARN] 5-2 output not found: ${FA}"
    fi
  fi

  # 9) YASS alignment and SVG generation (optional)
  if [[ "${RUN_YASS}" -eq 1 ]]; then
    if [[ "${RUN_TARGET}" -ne 1 ]]; then
      echo "[WARN] RUN_YASS=1 requires RUN_TARGET=1; skipping YASS for ${PREFIX}.${BARCODE}"
    elif [[ -z "${POSITIVE_OUT}" || ! -f "${POSITIVE_OUT}" ]]; then
      echo "[WARN] positive not found; skip YASS for ${PREFIX}.${BARCODE}"
    else
      POS_DIR="$(dirname "$(readlink -f "${POSITIVE_OUT}")")"
      YDIR="${POS_DIR}/${YASS_PREFIX}_${BARCODE}"
      mkdir -p "${YDIR}/fa" "${YDIR}/yop" "${YDIR}/svg"
      MAP_TSV="${YDIR}/id_map.tsv"
      [[ -f "${MAP_TSV}" ]] || echo -e "short_id\tfull_id" > "${MAP_TSV}"

      # positive to FASTA
      POS_FASTA="${YDIR}/positive.fa"
      tail -n +2 "${POSITIVE_OUT}" | awk -F'\t' '
        {
          full=$2; seq=$3;
          sub(/\r$/, "", full); sub(/\r$/, "", seq);
          print ">" full "\n" seq
        }' > "${POS_FASTA}"

      # Split per sequence
      seqkit split -i -O "${YDIR}/fa" "${POS_FASTA}"

      find "${YDIR}/fa" -type f -name "*.fa" -print0 \
        | xargs -0 -P "${THREADS}" -I{} bash -c '
            fa="{}"
            full_id=$(head -n 1 "$fa" | sed "s/^>//" | tr -d "\r")
            h=$(printf "%s" "$full_id" | md5sum | awk "{print substr(\$1,1,8)}")
            prefix=$(printf "%s" "$full_id" | tr -cd "A-Za-z0-9._-@" | cut -c1-40)
            short="${prefix}_${h}"
            printf "%s\t%s\n" "$short" "$full_id" >> "'"${MAP_TSV}"'"

            yop="'"${YDIR}"'/yop/${short}.yop"
            svg="'"${YDIR}"'/svg/${short}.svg"

            yass "$fa" "'"${REF}"'" -o "$yop" >/dev/null 2>&1 || exit 0
            python3 ./NSA/yass_yop_to_svg.py "$yop" "$svg"
          '
    fi
  fi

  # 10) BLASTN (optional)
  if [[ "${RUN_BLAST}" -eq 1 ]]; then
    FA="${MERGED_FQ%.fastq}.fa"
    [[ -f "${FA}" ]] || seqkit fq2fa "${MERGED_FQ}" > "${FA}"

    REF_TSV="${PREFIX}.${BARCODE}.${REFBASE}.blast.tsv"
    blastn -task "${BLAST_TASK}" -db "${DB_REF_PREFIX}" -query "${FA}" \
           -num_threads "${BLAST_THREADS}" -evalue "${BLAST_EVALUE}" \
           -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" -outfmt "${BLAST_OUTFMT}" \
           > "${REF_TSV}"

    ACT_TSV="${PREFIX}.${BARCODE}.ACT1.blast.tsv"
    blastn -task "${BLAST_TASK}" -db "${DB_ACT_PREFIX}" -query "${FA}" \
           -num_threads "${BLAST_THREADS}" -evalue "${BLAST_EVALUE}" \
           -max_target_seqs "${BLAST_MAX_TARGET_SEQS}" -outfmt "${BLAST_OUTFMT}" \
           > "${ACT_TSV}"

    TOTAL=$(grep -c "^>" "${FA}" || echo 0)
    REF_HITS=$(awk "{print \$1}" "${REF_TSV}" | sort -u | wc -l | tr -d " ")
    ACT_HITS=$(awk "{print \$1}" "${ACT_TSV}" | sort -u | wc -l | tr -d " ")

    REF_TOTAL=$(awk -v r="${REF_HITS}" -v t="${TOTAL}" 'BEGIN{if(t>0) printf "%.6f", r/t; else print 0}')
    ACT_TOTAL=$(awk -v r="${ACT_HITS}" -v t="${TOTAL}" 'BEGIN{if(t>0) printf "%.6f", r/t; else print 0}')
    REF_ACT=$(awk -v r="${REF_HITS}" -v a="${ACT_HITS}" 'BEGIN{if(a>0) printf "%.6f", r/a; else print "Inf"}')

    echo "${PREFIX}.${BARCODE},${TOTAL},${REF_HITS},${ACT_HITS},${REF_TOTAL},${ACT_TOTAL},${REF_ACT}" >> "${BLAST_SUMMARY}"
  fi

  echo "[*] Finished processing ${PREFIX}.${BARCODE}"
done

echo "[*] All analyses completed successfully."
