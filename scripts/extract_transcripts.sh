#!/usr/bin/env bash
# extract_transcripts.sh
# ----------------------
# Extract translated (amino acid) transcripts from a genome + GFF3
# for a specified genomic interval using gffread.
#
# Usage:
#   ./extract_transcripts.sh \
#     --genome C.lasius_genome_NCBI.fasta \
#     --gff Lasius_annotation.renamed.gff3 \
#     --region Chrom6:3343495-4962302 \
#     --out-prefix chr6_region \
#     --outdir results/ \
#     [--gffread-cmd /path/to/gffread]
#
# Outputs:
#   <outdir>/<out-prefix>_aa_transcripts.fa  (amino-acid FASTA)
#
# Notes:
#   - Region format must match the reference naming in your GFF/FASTA

set -euo pipefail

# -------- Args --------
GENOME=""
GFF=""
REGION=""
OUT_PREFIX="transcripts"
OUTDIR="."
GFFREAD_CMD="gffread"

print_help() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
}

# Simple arg parser
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome) GENOME="${2:-}"; shift 2 ;;
    --gff) GFF="${2:-}"; shift 2 ;;
    --region) REGION="${2:-}"; shift 2 ;;
    --out-prefix) OUT_PREFIX="${2:-}"; shift 2 ;;
    --outdir) OUTDIR="${2:-}"; shift 2 ;;
    --gffread-cmd) GFFREAD_CMD="${2:-}"; shift 2 ;;
    -h|--help) print_help; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; print_help; exit 1 ;;
  esac
done

# -------- Validation --------
if ! command -v "$GFFREAD_CMD" >/dev/null 2>&1; then
  echo "ERROR: gffread not found (tried '$GFFREAD_CMD')." >&2
  exit 1
fi

if [[ -z "${GENOME}" || -z "${GFF}" || -z "${REGION}" ]]; then
  echo "ERROR: --genome, --gff, and --region are required." >&2
  print_help
  exit 1
fi

if [[ ! -f "${GENOME}" ]]; then
  echo "ERROR: Genome FASTA not found: ${GENOME}" >&2
  exit 1
fi

if [[ ! -f "${GFF}" ]]; then
  echo "ERROR: GFF3 file not found: ${GFF}" >&2
  exit 1
fi

# Normalize/prepare output directory
mkdir -p "${OUTDIR}"

# -------- Outputs --------
AA_OUT="${OUTDIR%/}/${OUT_PREFIX}_aa_transcripts.fa"

echo "[INFO] Genome:   ${GENOME}"
echo "[INFO] GFF3:     ${GFF}"
echo "[INFO] Region:   ${REGION}"
echo "[INFO] Outdir:   ${OUTDIR}"
echo "[INFO] Out aa:   ${AA_OUT}"
echo "[INFO] gffread:  ${GFFREAD_CMD}"

# -------- Extract translated (amino acid) transcripts --------
"$GFFREAD_CMD" -y "${AA_OUT}" -g "${GENOME}" -r "${REGION}" "${GFF}"

echo "[DONE] Wrote:"
echo "  - ${AA_OUT}"
