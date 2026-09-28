#!/usr/bin/env bash
#-----------------------------------------------------------------------
# check_tables.sh -- does every table in the paper match numbers that are
# actually in results/?
#
#     bash check_tables.sh
#
# For each table it prints the rows as the manuscript has them, then the same
# rows recomputed from every candidate directory under results/. You compare
# by eye. That is the point: this is a crude check, and a crude check that
# runs is worth more than a precise one that does not.
#
# Why it exists
# -------------
# Table 2 was published with figures that could not be traced to anything in
# results/. The run had been overwritten, nothing recorded which model or
# sampler settings produced it, and the aggregated summaries did not match
# either. Nobody noticed until a referee-style pass asked the question. This
# script asks it on demand, before submission rather than after.
#
# A table with more than one candidate directory is not a pass. It means the
# repository cannot say which run is the published one, which is the state
# that lost Table 2.
#-----------------------------------------------------------------------

set -uo pipefail
cd "$(dirname "$0")"

# The manuscript source if it is beside this package, and otherwise the
# extract of its four table environments that ships with it.  The replication
# repository does not carry the paper, so without the fallback every "paper"
# block printed `awk: can t open file` and the check recomputed rows with
# nothing to compare them against -- which looks like a pass.
MS="../reviewed_manuscript/jcgs.tex"
if [ ! -f "$MS" ]; then MS="paper_tables.tex"; fi
if [ ! -f "$MS" ]; then
  echo "Neither ../reviewed_manuscript/jcgs.tex nor paper_tables.tex found;" >&2
  echo "there is nothing to compare the recomputed rows against." >&2
  exit 2
fi
echo "Comparing against: $MS"
R_OK=1
command -v Rscript >/dev/null 2>&1 || R_OK=0
[ "$R_OK" = 0 ] && echo "Rscript not found: only the manuscript side will be shown." && echo

hr() { printf '\n%s\n' "======================================================================"; }
paper_rows() {  # $1 = label
  local lab="$1"
  echo "  paper (jcgs.tex, table $lab):"
  # Walk out to the enclosing table environment and print the body rows inside
  # it.  Looking a fixed number of lines BACK from the label, as a first version
  # did, found nothing for Table 4, whose caption and label precede its tabular:
  # the check reported no rows for a table that was there all along.
  awk -v lab="$lab" '
    { buf[NR]=$0 }
    END {
      for (i=1;i<=NR;i++) if (buf[i] ~ ("label\\{" lab "\\}")) { at=i; break }
      if (!at) { print "    (label not found)"; exit }
      for (b=at; b>1  && buf[b] !~ /\\begin\{table/; b--) ;
      for (e=at; e<NR && buf[e] !~ /\\end\{table/;   e++) ;
      for (i=b;i<=e;i++)
        if (buf[i] ~ /\\\\ *$/ && buf[i] !~ /hline|multicolumn|backslash/)
          print "    " buf[i]
    }' "$MS" | sed 's/[[:space:]]\+/ /g'
}

agg() {  # $1 = script, $2 = in dir, $3 = grep filter
  local out; out=$(mktemp -d)
  Rscript "scripts/$1" --in="$2" --out="$out" 2>/dev/null | grep -E "$3" | head -12
  rm -rf "$out"
}

#-----------------------------------------------------------------------
hr; echo "TABLE 1  table:AEC  dimension sweep"
paper_rows 'table:AEC'
if [ "$R_OK" = 1 ]; then
  for d in results/dimension_sweep*; do
    [ -d "$d" ] || continue
    echo "  results/${d##*/}:"
    agg 20_aggregate_dimension_sweep.R "$d" '^[0-9]+ +&|backslash' | sed 's/^/    /'
  done
fi

#-----------------------------------------------------------------------
hr; echo "TABLE 2  tab:support  support recovery"
paper_rows 'tab:support'
# The directory the supplement names as the provenance of this table.  The
# comparison is against this one; the others are listed as alternatives, which
# is what they are.  Note there is no --post-draws here: the MANIFEST in the
# directory carries the draw count, so the screen cannot silently revert to
# demanding zero divergences, as it did when this number lived on the command
# line and this script did not pass it.
DECLARED_T2=results/support_recovery_ad99_slope
if [ "$R_OK" = 1 ]; then
  echo "  recomputed from ${DECLARED_T2#results/} (the run the supplement names):"
  Rscript scripts/02_aggregate_support_recovery.R --in="$DECLARED_T2" \
          --out=$(mktemp -d) 2>&1 \
    | awk '/falling back/{print} /Table 2 body/{b=1} /Fit diagnostics/{b=0} b&&/ & /{print}' \
    | sed 's/^/    /'
  echo "  other candidate directories, for the record:"
  for d in results/support_recovery*; do
    [ -d "$d" ] || continue
    [ "$d" = "$DECLARED_T2" ] && continue
    printf '    %-44s ' "${d##*/}"
    Rscript scripts/02_aggregate_support_recovery.R --in="$d" --out=$(mktemp -d) 2>&1 \
      | awk '/Stan model/{sub(/.*: /,"");m=$0} /falling back/{w=" [no draw count]"}
             END{printf "%s%s\n", (m?m:"model unrecorded"), w}'
  done
fi

#-----------------------------------------------------------------------
hr; echo "TABLE 3  tab:sbc  simulation-based calibration"
paper_rows 'tab:sbc'
if [ "$R_OK" = 1 ] && [ -d results/sbc ]; then
  # No --in: this aggregator reads results/summaries/ and ignored the flag, so
  # passing one suggested a directory was being checked that was not.  And its
  # stderr is kept, because discarding it hid a crash that left this section
  # blank -- a table with nothing under it looked like a table with nothing to
  # check.
  echo "  recomputed from results/summaries (arms of results/sbc):"
  Rscript scripts/12_aggregate_sbc.R --out=$(mktemp -d) 2>&1 \
    | grep -E 'gi_hd|Error|error|^ +[0-9]+\.[0-9]+' | head -20 | sed 's/^/    /'
fi

#-----------------------------------------------------------------------
hr; echo "TABLE 4  tab:brfss  case study"
paper_rows 'tab:brfss'
DECLARED_T4=results/case_brfss_diab_gi_hd_slope_iter4000
if [ "$R_OK" = 1 ]; then
  # One array task per held-out environment, so a single run leaves one
  # loeo_folds.csv per fold.  Counting files therefore reported forty-eight
  # ambiguous candidates for one unambiguous run; what matters is how many RUNS
  # produced folds, which is the directory above `folds/`.
  runs=$(find results -name 'loeo_folds.csv' 2>/dev/null \
         | sed -n 's|/folds/fold_[0-9]*/loeo_folds.csv$||p' | sort -u)
  [ -z "$runs" ] && echo "  no loeo_folds.csv anywhere under results/"
  for r in $runs; do
    mark="  "
    [ "$r" = "$DECLARED_T4" ] && mark=" *"
    printf '%s%-52s %s folds\n' "$mark" "${r#results/}" \
      "$(find "$r/folds" -name 'loeo_folds.csv' | wc -l | tr -d ' ')"
  done
  echo "  (* is the run the supplement names; recomputed fractions:)"
  # The selection frequencies, which is what Table 4 reports.  An earlier
  # version grepped for any `name number` line and printed the predictive
  # coverage summary instead -- plausible-looking rows that the table does not
  # contain, which is worse than printing nothing.
  Rscript scripts/07_aggregate_case_study.R --in="$DECLARED_T4/folds" \
          --out=$(mktemp -d) 2>&1 \
    | awk '/selection_frequency/{f=1;next} f&&/^ *x[0-9]+ /{print}' \
    | sed 's/^/    /'
  echo "    (covariates are x1..x13 in the order the caption lists them)"
fi

#-----------------------------------------------------------------------
# Does the prose quote numbers the table does not contain?
#
# A table can be correct and the paragraph around it still wrong.  Table 2's
# prose quotes thirteen figures -- a recall range, a false discovery
# proportion, four coverages, three interval scores -- and when the table is
# re-run every one of them is a chance to contradict the table on the same
# page.  A referee reads both.  So: pull the numbers out of the paragraph, pull
# them out of the table body, and print the ones that appear in the prose and
# nowhere in the table.
#
# Not every flagged number is an error.  `0.95` is the nominal level and `2`
# comes from an equation reference; neither belongs to a cell.  A figure quoted
# from a deliberate comparison run is also flagged, because this script cannot
# know the difference.  The test to apply to each one is whether it is traceable
# to a run named in the provenance map of the supplement: if it is, it is
# accounted for; if it is not, it is the Table 2 failure again.  The output is a
# list to look at, not a verdict -- same contract as the rest of this script.
#-----------------------------------------------------------------------
prose_check() {  # $1 = table label, $2 = human name
  local lab="$1" name="$2"
  echo "  numbers quoted in the prose but absent from the table body:"
  awk -v lab="$lab" '
    # The table body: the tabular of the environment carrying this label.
    # The prose: from the end of that environment to the next sectioning
    # command, which is where the discussion of this table stops.
    { line[NR] = $0 }
    END {
      for (i = 1; i <= NR; i++) if (line[i] ~ ("label\\{" lab "\\}")) { at = i; break }
      if (!at) { print "    (label not found)"; exit }
      for (i = at; i > 1 && line[i] !~ /begin\{tabular\}/; i--) ;
      tab_start = i
      for (i = tab_start; i <= NR && line[i] !~ /end\{tabular\}/; i++) body = body " " line[i]
      for (i = at; i <= NR && line[i] !~ /end\{table\}/; i++) ;
      for (i = i + 1; i <= NR && line[i] !~ /\\(sub)*section/; i++) prose = prose " " line[i]
      n = split(prose, w, /[^0-9.]+/)
      for (j = 1; j <= n; j++) {
        v = w[j]
        sub(/\.$/, "", v)
        if (v !~ /^[0-9]*\.[0-9]+$/) continue          # only decimals
        bare = v; sub(/^0/, "", bare)                   # .93 and 0.93 are one number
        if (index(body, v) || index(body, bare)) continue
        if (!seen[v]++) out = out "    " v "\n"
      }
      printf "%s", (out ? out : "    (none)\n")
    }' "$MS"
}

# Only against the manuscript.  paper_tables.tex holds the table environments
# and nothing between them, so "the prose after this table" is the next table,
# and every figure in it gets flagged.  A check that cries wolf on its own
# fallback teaches the reader to ignore it.
if [ "$MS" = "paper_tables.tex" ]; then
  hr; echo "PROSE AGAINST TABLES"
  echo "  Skipped: needs the manuscript source, and this package ships only the"
  echo "  table environments. Run it from a checkout that has jcgs.tex beside."
else
  hr; echo "PROSE AGAINST TABLES"
  prose_check 'tab:support' 'Table 2'
  prose_check 'table:AEC'   'Table 1'
fi

hr
cat <<'EOF'
A table passes only if exactly one directory reproduces it. More than one
candidate, or none, means the published figures cannot be traced, which is the
state Table 2 was found in.
EOF
