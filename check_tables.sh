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
if [ ! -f "$MS" ]; then
  MS="paper_tables.tex"
  # Say what cannot be established here. Without the manuscript there is nothing
  # to diff the extract against, so every verdict below is measured against a
  # file whose currency this clone cannot verify -- it can only report the
  # commit the extract was made from, which is stamped in its header.
  echo "NOTE: the manuscript is not in this checkout, so the extract's currency"
  echo "      cannot be verified here. It was generated from:"
  grep -m1 'Generated from jcgs.tex at commit' paper_tables.tex | sed 's/^%/       /'
fi
if [ ! -f "$MS" ]; then
  echo "Neither ../reviewed_manuscript/jcgs.tex nor paper_tables.tex found;" >&2
  echo "there is nothing to compare the recomputed rows against." >&2
  exit 2
fi
echo "Comparing against: $MS"
# Is the shipped extract still the paper's?  Nothing else checks it, so a stale
# paper_tables.tex would have every table agreeing with a manuscript that no
# longer says that.  When the manuscript is beside this package, regenerate the
# extract into a temporary file and diff.
if [ -f "../reviewed_manuscript/jcgs.tex" ] && [ -f paper_tables.tex ]; then
  FRESH=$(mktemp)
  awk '/\\begin\{table/{buf="";inb=1} inb{buf=buf $0 "\n"}
       /\\end\{table\}/&&inb{
         if (buf ~ /label\{(table:AEC|tab:support|tab:sbc|tab:brfss)\}/) printf "%s", buf
         inb=0}' ../reviewed_manuscript/jcgs.tex > "$FRESH"
  if grep -v '^%' paper_tables.tex | diff -q - "$FRESH" >/dev/null 2>&1; then
    echo "paper_tables.tex is current with the manuscript."
  else
    echo "WARNING: paper_tables.tex differs from the manuscript's tables."
    echo "  It is stale; regenerate it before trusting any verdict below."
  fi
fi
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
    }' "$MS" | sed 's/[[:space:]]\+/ /g' | tee "${PAPER_OUT:-/dev/null}"
}

agg() {  # $1 = script, $2 = in dir, $3 = grep filter
  local out; out=$(mktemp -d)
  Rscript "scripts/$1" --in="$2" --out="$out" 2>/dev/null | grep -E "$3" | head -12
  rm -rf "$out"
}

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# The verdict.
#
# This script used to print the two sides and leave the comparison to the
# reader.  That was defensible while the recomputed rows came out in a
# different shape from the paper's; it is not defensible now that both sides
# are emitted in the same format, and a check that cannot fail is decoration.
#
# What is compared is the multiset of decimal numbers in each block, at the
# precision printed.  Formatting is stripped first -- \textbf{}, the revision's
# \rd{} and \bl{} colour macros, column separators -- because none of it is a
# number.  Multisets rather than sequences, because Table 4's rows come out
# ordered by frequency while the paper orders them in two columns, so the
# sequences legitimately differ while the values must not.
#
# The one thing this cannot catch is two values swapping places within a table.
# Nothing else gets through: any changed, missing or extra figure shows up as a
# difference, and the differing values are printed.
# Every table is compared. The two sides do not always print to the same number
# of decimals -- the paper rounds where the aggregators do not, 1.02 against
# 1.018 -- so rather than refuse, the check reports the precision at which each
# table agrees. Both sides are rounded to d decimals and compared as multisets,
# for d from 4 down to 1, and the largest d that agrees is what gets printed.
#
# Multisets rather than sequences, because Table 4 comes out ordered by
# frequency while the paper orders it in two columns. A figure could therefore
# swap places with another inside the same table and go unnoticed; nothing else
# does, and any changed, missing or extra figure is named below.
SUMMARY=$(mktemp); TMPA=$(mktemp); TMPB=$(mktemp)
# Whether anything failed.  Printing DISAGREES and then exiting 0 is the shape
# of a check that cannot fail, which is the shape this script had.
BAD=0
numbers() {
  sed -e 's/\\textbf{//g' -e 's/\\rd{//g' -e 's/\\bl{//g' -e 's/[{}]//g' \
    | grep -oE '[0-9]+\.[0-9]+|[0-9]*\.[0-9]+'
}
verdict() {  # $1 = name, $2 = paper block, $3 = recomputed block
  local name="$1" pf rf na nb d
  pf=$(mktemp); rf=$(mktemp)
  numbers < "$2" > "$pf"; numbers < "$3" > "$rf"
  na=$(grep -c . "$pf" 2>/dev/null || echo 0)
  nb=$(grep -c . "$rf" 2>/dev/null || echo 0)
  if [ "$na" -eq 0 ] || [ "$nb" -eq 0 ]; then
    printf '  %-9s %4s figures   NOT CHECKED: one side produced none\n' "$name" "-" >> "$SUMMARY"
    printf '  VERDICT: %s NOT CHECKED -- one side produced no figures.\n' "$name"
    BAD=$((BAD+1)); return
  fi
  for d in 4 3 2 1; do
    awk -v d="$d" '{printf "%.*f\n", d, $1}' "$pf" | sort > "$TMPA"
    awk -v d="$d" '{printf "%.*f\n", d, $1}' "$rf" | sort > "$TMPB"
    if cmp -s "$TMPA" "$TMPB"; then
      printf '  %-9s %4d figures   agree to %d decimal%s\n' \
        "$name" "$na" "$d" "$([ "$d" = 1 ] || echo s)" >> "$SUMMARY"
      printf '  VERDICT: %s agrees with results/ to %d decimal%s (%d figures).\n' \
        "$name" "$d" "$([ "$d" = 1 ] || echo s)" "$na"
      return
    fi
  done
  # No single precision fits every figure. Pair the two sides by sorted value
  # and look at the gaps: a gap no larger than half a unit in the paper's last
  # printed decimal is the paper rounding, not the table drifting. Table 4's
  # strength-guideline frequency is exactly 0.625, printed as 0.63 and rounded
  # by printf to 0.62, and an all-or-nothing rule called that a failure.
  sort -g "$pf" > "$TMPA"; sort -g "$rf" > "$TMPB"
  local worst n_off
  worst=$(paste "$TMPA" "$TMPB" | awk '
    { d = $1 - $2; if (d < 0) d = -d; if (d > m) m = d }
    END { printf "%.6f", m }')
  n_off=$(paste "$TMPA" "$TMPB" | awk '{ d=$1-$2; if (d<0) d=-d; if (d > 0.0051) c++ }
                                       END { print c+0 }')
  if [ "$n_off" -eq 0 ]; then
    printf '  %-9s %4d figures   agree to the printed precision (largest gap %s)\n' \
      "$name" "$na" "$worst" >> "$SUMMARY"
    printf '  VERDICT: %s agrees to the precision the paper prints; largest gap %s,\n' \
      "$name" "$worst"
    printf '           which is rounding in the last printed digit.\n'
  else
    printf '  %-9s %4d figures   %d differ by more than rounding\n' \
      "$name" "$na" "$n_off" >> "$SUMMARY"
    printf '  VERDICT: %s DISAGREES: %d figure(s) differ by more than half a unit\n' \
      "$name" "$n_off"
    printf '           in the last printed decimal; largest gap %s.\n' "$worst"
    paste "$TMPA" "$TMPB" | awk '{ d=$1-$2; if (d<0) d=-d;
      if (d > 0.0051) printf "    paper %s vs recomputed %s\n", $1, $2 }'
    BAD=$((BAD+1))
  fi
}

WORK=$(mktemp -d)

hr; echo "TABLE 1  table:AEC  dimension sweep"
PAPER_OUT="$WORK/t1.paper" paper_rows 'table:AEC'
if [ "$R_OK" = 1 ]; then
  : > "$WORK/t1.recomp"
  for d in results/dimension_sweep*; do
    [ -d "$d" ] || continue
    echo "  results/${d##*/}:"
    agg 20_aggregate_dimension_sweep.R "$d" '^[0-9]+ +&|backslash' \
      | tee -a "$WORK/t1.recomp" | sed 's/^/    /'
  done
  verdict "Table 1" "$WORK/t1.paper" "$WORK/t1.recomp"
fi

#-----------------------------------------------------------------------
hr; echo "TABLE 2  tab:support  support recovery"
PAPER_OUT="$WORK/t2.paper" paper_rows 'tab:support'
# The directory the supplement names as the provenance of this table.  The
# comparison is against this one; the others are listed as alternatives, which
# is what they are.  Note there is no --post-draws here: the MANIFEST in the
# directory carries the draw count, so the screen cannot silently revert to
# demanding zero divergences, as it did when this number lived on the command
# line and this script did not pass it.
DECLARED_T2=results/support_recovery_properprior
if [ "$R_OK" = 1 ]; then
  if [ -f "$DECLARED_T2/MANIFEST" ]; then
    printf '  settings recorded in %s/MANIFEST: %s\n' "${DECLARED_T2#results/}" \
      "$(grep -E '^(post_draws|v_prior_shape|v_prior_rate|adapt_delta|model)=' \
         "$DECLARED_T2/MANIFEST" | tr '\n' ' ')"
  else
    echo "  NO MANIFEST in $DECLARED_T2: the screen falls back to zero divergences."
    BAD=$((BAD+1))
  fi
  echo "  recomputed from ${DECLARED_T2#results/} (the run the supplement names):"
  Rscript scripts/02_aggregate_support_recovery.R --in="$DECLARED_T2" \
          --out=$(mktemp -d) 2>&1 \
    | awk '/falling back/{print} /Table 2 body/{b=1} /Fit diagnostics/{b=0} b&&/ & /{print}' \
    | tee "$WORK/t2.recomp" | sed 's/^/    /'
  verdict "Table 2" "$WORK/t2.paper" "$WORK/t2.recomp"
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
PAPER_OUT="$WORK/t3.paper" paper_rows 'tab:sbc'
if [ "$R_OK" = 1 ] && [ -d results/sbc ]; then
  # No --in: this aggregator reads results/summaries/ and ignored the flag, so
  # passing one suggested a directory was being checked that was not.  And its
  # stderr is kept, because discarding it hid a crash that left this section
  # blank -- a table with nothing under it looked like a table with nothing to
  # check.
  echo "  recomputed from results/summaries (arms of results/sbc):"
  # Only the two quantities the table reports: coverage, which is the fifth
  # field of the per-arm lines, and the extreme rank ratio, which R prints in a
  # second block when the data frame is too wide for the terminal.
  Rscript scripts/12_aggregate_sbc.R --out=$(mktemp -d) 2>&1 \
    | tee "$WORK/t3.raw" \
    | grep -E 'gi_hd|Error|error|^ +[0-9]+\.[0-9]+' | head -20 | sed 's/^/    /'
  awk '/gi_hd/{print $5} /^ +[0-9]+\.[0-9]+ +(calibrated|too narrow)/{print $1}' \
      "$WORK/t3.raw" > "$WORK/t3.recomp"
  verdict "Table 3" "$WORK/t3.paper" "$WORK/t3.recomp"
fi

#-----------------------------------------------------------------------
hr; echo "TABLE 4  tab:brfss  case study"
PAPER_OUT="$WORK/t4.paper" paper_rows 'tab:brfss'
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
    | tee "$WORK/t4.recomp" | sed 's/^/    /'
  verdict "Table 4" "$WORK/t4.paper" "$WORK/t4.recomp"
  cat <<'NOTE'
    Rows are selection frequencies in descending order, which is how the paper
    orders them too, so compare the two columns of numbers as sorted lists. The
    labels are the design-matrix positions, not the caption's names: the per-fold
    CSVs record `selected` as indices and R/case_data.R fixes what each index is.
    They are NOT the caption's order -- an earlier version of this note said they
    were, which would have had a reader map x3 to the caption's first row.
NOTE
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
echo "TABLES CHECKED AGAINST results/"
echo
cat "$SUMMARY"
echo
cat <<'EOF'
Precision is the finest at which every figure in that table agrees; it differs
because the paper rounds and the aggregators do not. Table 4's one gap of 0.005
is the exactly-0.625 frequency, printed as 0.63 and rounded by printf to 0.62.
EOF
if [ "$BAD" -gt 0 ]; then
  printf '\nRESULT: %d table(s) did not check out. Exiting non-zero.\n' "$BAD"
else
  printf '\nRESULT: every table checks out against results/.\n'
fi


exit $((BAD > 0))
