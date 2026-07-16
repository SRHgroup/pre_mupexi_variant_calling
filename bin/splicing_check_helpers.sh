#!/usr/bin/bash

# Shared output and scheduler-marker resolution for splicing status commands.

splicing_check_star_root() {
  if [ -n "${splicing_sjdir:-}" ]; then
    printf '%s\n' "$splicing_sjdir"
  elif [ -n "${splicing_stardir:-}" ]; then
    printf '%s\n' "$splicing_stardir"
  elif [ -n "${stardir:-}" ]; then
    printf '%s\n' "$stardir"
  elif [ -n "${bamdir:-}" ]; then
    printf '%s/star\n' "$(dirname "$bamdir")"
  else
    return 1
  fi
}

splicing_check_output_root() {
  if [ -n "${splicing_outdir:-}" ]; then
    printf '%s\n' "$splicing_outdir"
  elif [ -n "${datadir:-}" ]; then
    printf '%s/splicing\n' "${datadir%/}"
  else
    return 1
  fi
}

splicing_check_sample_id() {
  local patient="$1"
  printf '%s_RNA_%s\n' "$patient" "${tumor_tag:-TUMOR}"
}

splicing_check_is_rna_sample_id() {
  local sample_id="$1"
  case "$sample_id" in
    *"_${rna_tumor_label:-RNA_TUMOR}"|*"_${out_rna_tumor_label:-${rna_tumor_label:-RNA_TUMOR}}"|*_RNA_TUMOR|*_RNA_TUMOUR) return 0 ;;
    *) return 1 ;;
  esac
}

splicing_check_step_prefix() {
  case "$1" in
    spl1) printf '%s\n' "splicing_spl1" ;;
    spl2) printf '%s\n' "splicing_spl2" ;;
    spl3) printf '%s\n' "splicing_spl3" ;;
    spl3.5) printf '%s\n' "splicing_spl3_5" ;;
    spl4) printf '%s\n' "splicing_spl4" ;;
    *) return 1 ;;
  esac
}

splicing_check_step_outputs() {
  local patient="$1"
  local step="$2"
  local sample_id star_root output_root plain gz
  sample_id="$(splicing_check_sample_id "$patient")"

  case "$step" in
    spl1)
      star_root="$(splicing_check_star_root)" || return 1
      plain="${star_root%/}/${patient}/${sample_id}/${sample_id}.SJ.out.tab"
      gz="${plain}.gz"
      if [ -s "$plain" ]; then
        printf '%s\n' "$plain"
      elif [ -s "$gz" ]; then
        printf '%s\n' "$gz"
      elif [ -e "$plain" ]; then
        printf '%s\n' "$plain"
      elif [ -e "$gz" ]; then
        printf '%s\n' "$gz"
      else
        printf '%s\n' "$plain"
      fi
      ;;
    spl2)
      output_root="$(splicing_check_output_root)" || return 1
      printf '%s/%s/%s.spl2.novel_junctions.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      ;;
    spl3)
      output_root="$(splicing_check_output_root)" || return 1
      printf '%s/%s/%s.spl3.event_annotated.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      ;;
    spl3.5)
      output_root="$(splicing_check_output_root)" || return 1
      printf '%s/%s/%s.spl3.5.cancer_unique.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      printf '%s/%s/%s.spl3.5.normal_present.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      printf '%s/%s/%s.spl3.5.normal_filter_summary.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      ;;
    spl4)
      output_root="$(splicing_check_output_root)" || return 1
      printf '%s/%s/%s.spl4.neojunctions.tsv\n' "${output_root%/}" "$sample_id" "$sample_id"
      printf '%s/%s/%s.spl4.neojunctions.nt.fa\n' "${output_root%/}" "$sample_id" "$sample_id"
      printf '%s/%s/%s.spl4.neojunctions.aa.fa\n' "${output_root%/}" "$sample_id" "$sample_id"
      ;;
    *) return 1 ;;
  esac
}

splicing_check_output_requires_content() {
  local step="$1"
  local path="$2"
  if [ "$step" = "spl4" ]; then
    case "$path" in
      *.fa|*.fasta) return 1 ;;
    esac
  fi
  return 0
}

splicing_check_step_completion() {
  local patient="$1"
  local step="$2"
  local outputs out detail="" missing="" empty=""
  outputs="$(splicing_check_step_outputs "$patient" "$step")" || {
    printf 'CONFIG_ERROR\tunable to resolve %s paths\n' "$step"
    return 1
  }

  while IFS= read -r out; do
    [ -n "$out" ] || continue
    if [ ! -e "$out" ]; then
      [ -n "$missing" ] && missing="${missing} ; "
      missing="${missing}${out}"
    elif splicing_check_output_requires_content "$step" "$out" && [ ! -s "$out" ]; then
      [ -n "$empty" ] && empty="${empty} ; "
      empty="${empty}${out}"
    fi
    [ -n "$detail" ] && detail="${detail} ; "
    detail="${detail}${out}"
  done <<< "$outputs"

  if [ -n "$missing" ]; then
    printf 'MISSING\t%s\n' "$missing"
    return 1
  fi
  if [ -n "$empty" ]; then
    printf 'EMPTY\t%s\n' "$empty"
    return 1
  fi
  printf 'DONE\t%s\n' "$detail"
}

splicing_check_spl4_input() {
  local patient="$1"
  local sample_id output_root suffix
  sample_id="$(splicing_check_sample_id "$patient")"
  output_root="$(splicing_check_output_root)" || return 1
  suffix="${splicing_spl4_input_suffix:-.spl3.event_annotated.tsv}"
  printf '%s/%s/%s%s\n' "${output_root%/}" "$sample_id" "$sample_id" "$suffix"
}

splicing_check_step_input_status() {
  local patient="$1"
  local step="$2"
  local input_info input_state input_detail star_root patient_root shard spl4_input
  case "$step" in
    spl1)
      star_root="$(splicing_check_star_root)" || {
        printf 'NO_INPUT\tCONFIG needs splicing_sjdir, splicing_stardir, stardir, or bamdir\n'
        return 1
      }
      patient_root="${star_root%/}/${patient}"
      shard=""
      if [ -d "$patient_root" ]; then
        shard="$(find "$patient_root" -type f -name '*.[0-9][0-9][0-9][0-9].SJ.out.tab*' -size +0c -print -quit 2>/dev/null || true)"
      fi
      if [ -n "$shard" ]; then
        printf 'INPUT_OK\t%s\n' "$shard"
      else
        printf 'NO_INPUT\t%s/*.[0-9][0-9][0-9][0-9].SJ.out.tab*\n' "$patient_root"
        return 1
      fi
      ;;
    spl2|spl3)
      if [ "$step" = "spl2" ]; then
        input_info="$(splicing_check_step_completion "$patient" spl1)" || true
      else
        input_info="$(splicing_check_step_completion "$patient" spl2)" || true
      fi
      input_state="${input_info%%$'\t'*}"
      input_detail="${input_info#*$'\t'}"
      if [ "$input_state" = "DONE" ]; then
        printf 'INPUT_OK\t%s\n' "$input_detail"
      else
        printf 'NO_INPUT\t%s\n' "$input_detail"
        return 1
      fi
      ;;
    spl3.5)
      input_info="$(splicing_check_step_completion "$patient" spl3)" || true
      input_state="${input_info%%$'\t'*}"
      input_detail="${input_info#*$'\t'}"
      if [ "$input_state" != "DONE" ]; then
        printf 'NO_INPUT\t%s\n' "$input_detail"
        return 1
      fi
      if [ -z "${splicing_normal_junction_ref:-}" ] || [ ! -s "$splicing_normal_junction_ref" ]; then
        printf 'NO_INPUT\t%s\n' "${splicing_normal_junction_ref:-CONFIG:splicing_normal_junction_ref missing}"
        return 1
      fi
      printf 'INPUT_OK\t%s ; %s\n' "$input_detail" "$splicing_normal_junction_ref"
      ;;
    spl4)
      spl4_input="$(splicing_check_spl4_input "$patient")" || {
        printf 'NO_INPUT\tCONFIG needs splicing_outdir or datadir\n'
        return 1
      }
      if [ -s "$spl4_input" ]; then
        printf 'INPUT_OK\t%s\n' "$spl4_input"
      else
        printf 'NO_INPUT\t%s\n' "$spl4_input"
        return 1
      fi
      ;;
    *)
      printf 'NO_INPUT\tunknown step: %s\n' "$step"
      return 1
      ;;
  esac
}

splicing_check_step_logdir() {
  local step="$1"
  local prefix root
  prefix="$(splicing_check_step_prefix "$step")" || return 1
  if [ -n "${splicing_logroot:-}" ]; then
    root="$splicing_logroot"
  elif [ "$step" = "spl1" ]; then
    root="$(splicing_check_star_root)" || return 1
    root="${root%/}/${prefix}.logs_and_reports"
  else
    root="$(splicing_check_output_root)" || return 1
    root="${root%/}/${prefix}.logs_and_reports"
  fi
  printf '%s/logs\n' "${root%/}"
}

splicing_check_find_marker() {
  local patient="$1"
  local sample_id="$2"
  local step="$3"
  local prefix logdir tag marker fallback=""
  prefix="$(splicing_check_step_prefix "$step")" || return 1
  logdir="$(splicing_check_step_logdir "$step")" || return 1

  for tag in "$patient" "$sample_id"; do
    tag="${tag//[^A-Za-z0-9_.-]/_}"
    marker="${logdir}/submitted.${prefix}.${tag}.jobid"
    [ -n "$fallback" ] || fallback="$marker"
    if [ -f "$marker" ]; then
      printf '%s\n' "$marker"
      return 0
    fi
  done
  printf '%s\n' "$fallback"
}
