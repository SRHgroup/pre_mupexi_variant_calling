#!/usr/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash research/test_single_vcf.sh \
    -i /path/to/phased.vcf.gz \
    -p SAMPLE_ID \
    [-o OUTDIR] [--max-distance N] [--min-cluster-size N] [--rna-label LABEL]

Example:
  bash research/test_single_vcf.sh \
    -i /home/projects/SRHgroup/projects/SingelCell_Bladder/data/rna/rnadnavar/variant_calling/mutect2/RNA_TUMOUR_vs_DNA_NORMAL/DNAt_DNAn_RNAt_merged_genotyped_TUMOR_phased.vcf.gz \
    -p SingelCell_Test \
    -o research/output_single
USAGE
}

outdir="research/output_single"
max_distance=500
min_cluster_size=2
rna_label="RNA_TUMOUR"

while :; do
  case ${1:-} in
    -i|--input) vcf="$2"; shift 2 ;;
    -p|--patient) patient="$2"; shift 2 ;;
    -o|--outdir) outdir="$2"; shift 2 ;;
    --max-distance) max_distance="$2"; shift 2 ;;
    --min-cluster-size) min_cluster_size="$2"; shift 2 ;;
    --rna-label) rna_label="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) break ;;
  esac
done

[ -n "${vcf:-}" ] || { usage; exit 1; }
[ -n "${patient:-}" ] || { usage; exit 1; }
[ -f "$vcf" ] || { echo "ERROR: input VCF not found: $vcf" >&2; exit 1; }

mkdir -p "$outdir"

python3 "$(dirname "$0")/extract_rna_editing_clusters.py" \
  --input "${patient}=${vcf}" \
  --out-variants "${outdir}/rna_edit_variants.tsv" \
  --out-clusters "${outdir}/rna_edit_clusters.tsv" \
  --max-distance "$max_distance" \
  --min-cluster-size "$min_cluster_size" \
  --rna-label "$rna_label"

python3 "$(dirname "$0")/plot_rna_editing_clusters.py" \
  --variants "${outdir}/rna_edit_variants.tsv" \
  --clusters "${outdir}/rna_edit_clusters.tsv" \
  --out-prefix "${outdir}/rna_edit_cluster_landscape"

echo "[done] one-sample test complete in: $outdir"
