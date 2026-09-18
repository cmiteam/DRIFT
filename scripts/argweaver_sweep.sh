#!/usr/bin/env bash
#
# ARGweaver Ne-prior sweep driver — TMR4A.md W9 / Study Design E1.
#
# Runs on the REMOTE box, inside WSL Ubuntu, on files carried across the RDP
# session. DRIFT never invokes this; see "DRIFT to ARGweaver workflow.docx" §5-§7
# for how the files get here.
#
# WHAT IT IS FOR. The constant-N coalescent predicts TMR4A ~ N generations before
# any data is consulted, so a single assumed --popsize mostly reproduces its own
# assumption. The deliverable is inferred depth as a FUNCTION of assumed N, and
# the ratio inferred/N is the result: flat near 1 means the statistic is largely
# restating its prior.
#
# THE TWO RATES ARE NEVER RETYPED. --mutrate and --recombrate are lifted out of
# the DRIFT-generated <prefix>_argweaver_cmd.sh, because they are the easiest
# thing in this pipeline to get wrong and nothing in the .sites file records them.
# If the cmd file is missing, this script stops rather than guessing.
#
# --maxtime IS FIXED ACROSS THE SWEEP, and this matters. arg-sample's own
# default is 200,000 generations over 20 log-spaced points; a DRIFT history is
# ~20-180 generations, so at the default every coalescence falls below the first
# grid point and the reported depth is discretisation plus prior rather than data.
# But tying maxtime to 4N PER RUN is also wrong: at N the prior's E[TMRCA] is
# 3.9N for 40 haplotypes, so a 4N ceiling truncates the prior at its own mean,
# and each sweep point would then carry a DIFFERENT truncation. The sweep would
# measure a mix of prior and grid ceiling instead of the prior alone. So maxtime
# is set ONCE, from 4x the LARGEST N in the sweep, leaving assumed-N as the only
# thing that varies. Override with -m.
#
# CAUTION ON THE STATISTIC. arg-summarize --tmrca is the time to the FULL MRCA,
# i.e. TMR-1-A. It is NOT the TMR-4-A this study is about, and nothing in
# arg-summarize computes TMR-K-A for arbitrary K. This script reports --tmrca as
# a first-order signal and, under -T, additionally saves the per-sample newick
# from --tree, which is what W7's parser needs to walk back to K lineages.
#
# Usage:
#   ./argweaver_sweep.sh [options] <prefix>_chr<N>.sites
#
# Options:
#   -N "LIST"  assumed effective sizes to sweep   (default: 500 1000 2000 5000 10000 20000)
#   -r RANGE   restrict to a sub-region, e.g. 1-2000000  (default: whole .sites region)
#   -i ITERS   arg-sample --iters                 (default: 1000, arg-sample's own)
#   -p STEP    arg-sample --sample-step           (default: 10, arg-sample's own)
#   -t NTIMES  arg-sample --ntimes                (default: 30)
#   -S SEED    arg-sample --randseed              (default: 4242)
#   -m MAXTIME arg-sample --maxtime, fixed for the whole sweep
#              (default: 4 x the largest value in -N)
#   -o DIR     output directory                   (default: ./sweep)
#   -T         also emit per-sample newick (--tree) for W7's parser
#   -c         preflight checks only, run nothing
#
# Exit codes: 0 all runs completed; 1 preflight failed; 2 one or more runs failed.

set -u -o pipefail

POPSIZES="500 1000 2000 5000 10000 20000"
REGION=""
ITERS=1000
STEP=10
NTIMES=30
SEED=4242
MAXTIME=""
OUTDIR="sweep"
EMIT_TREES=0
CHECK_ONLY=0

while getopts "N:r:i:p:t:S:m:o:Tch" opt; do
  case "$opt" in
    N) POPSIZES="$OPTARG" ;;
    r) REGION="$OPTARG" ;;
    i) ITERS="$OPTARG" ;;
    p) STEP="$OPTARG" ;;
    t) NTIMES="$OPTARG" ;;
    S) SEED="$OPTARG" ;;
    m) MAXTIME="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    T) EMIT_TREES=1 ;;
    c) CHECK_ONLY=1 ;;
    h) sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "run '$0 -h' for usage" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [ $# -ne 1 ]; then
  echo "ERROR: expected exactly one .sites file. Run '$0 -h' for usage." >&2
  exit 1
fi

SITES="$1"
SITES_DIR=$(cd "$(dirname "$SITES")" && pwd)
SITES_BASE=$(basename "$SITES")
SITES="$SITES_DIR/$SITES_BASE"

fail=0
note() { printf '  %s\n' "$*"; }
bad()  { printf '  FAIL  %s\n' "$*"; fail=1; }
warn() { printf '  WARN  %s\n' "$*"; }

echo "=== Preflight ==="

# ---------------------------------------------------------------- tools
# smc2bed-all shells out to three tools that are NOT part of ARGweaver:
# bgzip and tabix from the tabix package, sort-bed from bedops. Missing them
# fails deep inside the converter with an error that looks like data corruption.
for t in arg-sample arg-summarize smc2bed smc2bed-all bgzip tabix sort-bed; do
  if command -v "$t" >/dev/null 2>&1; then
    note "found    $t"
  else
    case "$t" in
      bgzip|tabix)  bad "$t missing — sudo apt install -y tabix" ;;
      sort-bed)     bad "sort-bed missing — sudo apt install -y bedops" ;;
      *)            bad "$t missing — is ~/argweaver/bin on PATH?" ;;
    esac
  fi
done

# The mdrasmus fork has no smc2bed-all and breaks the whole BED step. Do not use
# the version string to tell the forks apart: arg-summarize reports 0.3 in both.
if [ -d "$HOME/argweaver/.git" ]; then
  remote=$(git -C "$HOME/argweaver" remote get-url origin 2>/dev/null || echo "?")
  case "$remote" in
    *CshlSiepelLab*) note "argweaver repo is CshlSiepelLab (correct)" ;;
    *mdrasmus*)      bad "argweaver repo is the mdrasmus fork — no smc2bed-all; re-clone CshlSiepelLab" ;;
    *)               warn "argweaver repo origin is '$remote' — expected CshlSiepelLab" ;;
  esac
fi

# ---------------------------------------------------------------- filesystem
# Cross-filesystem I/O over /mnt/c is markedly slower and this is a long job.
case "$SITES_DIR" in
  /mnt/*) warn "running under $SITES_DIR — copy into \$HOME first; /mnt/c is slow for a long MCMC" ;;
esac

# ---------------------------------------------------------------- input files
[ -r "$SITES" ] || bad "cannot read $SITES"

# The .sites file crossed from Windows. CRLF misbehaves here in ways that look
# like data corruption, so detect it rather than letting it surface downstream.
if [ -r "$SITES" ] && head -c 4096 "$SITES" | grep -q $'\r'; then
  bad "$SITES_BASE has CRLF line endings — run: dos2unix *.sh *.sites *.bed *.txt"
fi

# Derive the DRIFT-generated command file: <run>_chr<N>.sites -> <run>_argweaver_cmd.sh
RUN_PREFIX=$(printf '%s' "$SITES_BASE" | sed -E 's/_chr[^_.]+\.sites$//')
CMDFILE="$SITES_DIR/${RUN_PREFIX}_argweaver_cmd.sh"
if [ ! -r "$CMDFILE" ]; then
  alt=$(ls "$SITES_DIR"/*_argweaver_cmd.sh 2>/dev/null | head -1)
  [ -n "$alt" ] && CMDFILE="$alt"
fi
[ -r "$CMDFILE" ] || bad "no *_argweaver_cmd.sh beside the .sites file — the rates live there and are not being guessed"

# The scaling rationale is the record of what scaling was claimed. Study Design
# §5.3 RECORD: confirm it ends in NO WARNINGS.
SCALEFILE="$SITES_DIR/${RUN_PREFIX}_argweaver_scaling.txt"
if [ -r "$SCALEFILE" ]; then
  if grep -q "NO WARNINGS" "$SCALEFILE"; then
    ratio=$(grep -E "^\s*ratio" "$SCALEFILE" | head -1 | tr -s ' ')
    note "scaling file: NO WARNINGS —${ratio#*=}"
  else
    warn "scaling file contains WARNING lines — an off-anchor run is legitimate but must not be accidental:"
    grep -n "WARNING" "$SCALEFILE" | sed 's/^/        /'
  fi
else
  warn "no ${RUN_PREFIX}_argweaver_scaling.txt — without it these numbers are uninterpretable later"
fi

# ---------------------------------------------------------------- .sites shape
if [ -r "$SITES" ]; then
  nhap=$(head -1 "$SITES" | awk '{print NF-1}')
  hdr1=$(head -1 "$SITES" | cut -f1)
  hdr2=$(sed -n '2p' "$SITES" | cut -f1)
  nlines=$(wc -l < "$SITES")
  nsites=$((nlines - 2))
  [ "$hdr1" = "NAMES" ]  || bad "line 1 is '$hdr1', expected NAMES"
  [ "$hdr2" = "REGION" ] || bad "line 2 is '$hdr2', expected REGION"
  regline=$(sed -n '2p' "$SITES")
  note "haplotypes  $nhap   (20 sampled diploids = 40)"
  note "sites       $nsites"
  note "region      $(printf '%s' "$regline" | cut -f2-4 | tr '\t' ' ')"

  # Positions must be strictly ascending; arg-sample does not check and the
  # failure downstream is silent rather than loud.
  desc=$(awk 'NR>2 { if ($1 <= prev) { print NR; exit } prev=$1 }' "$SITES")
  if [ -n "$desc" ]; then
    bad "positions not strictly ascending — first offender at line $desc"
  else
    note "positions   strictly ascending"
  fi

  # Site density is NOT covered by the scaling check (Study Design §8): the
  # rationale confirms the RATES, while realized diversity is set separately by
  # seeding. Check it by hand, every run.
  rstart=$(printf '%s' "$regline" | cut -f3)
  rend=$(printf '%s' "$regline" | cut -f4)
  if [ -n "$rstart" ] && [ -n "$rend" ] && [ "$nsites" -gt 0 ]; then
    note "density     1 segregating site per $(( (rend - rstart + 1) / nsites )) bp"
  fi
fi

# ---------------------------------------------------------------- rates
MUTRATE=""; RECOMBRATE=""
if [ -r "$CMDFILE" ]; then
  # Unfold backslash continuations, then take the arg-sample line for THIS
  # .sites file — a multi-region cmd file carries one --recombrate per chromosome.
  line=$(sed -e :a -e '/\\$/{N;s/\\\n//;ta' -e '}' "$CMDFILE" | grep -F -- "--sites $SITES_BASE" | head -1)
  if [ -z "$line" ]; then
    bad "no arg-sample line in $(basename "$CMDFILE") mentions $SITES_BASE"
  else
    MUTRATE=$(printf '%s\n' "$line"    | awk '{for(i=1;i<NF;i++) if($i=="--mutrate")    print $(i+1)}')
    RECOMBRATE=$(printf '%s\n' "$line" | awk '{for(i=1;i<NF;i++) if($i=="--recombrate") print $(i+1)}')
    [ -n "$MUTRATE" ]    || bad "could not read --mutrate from $(basename "$CMDFILE")"
    [ -n "$RECOMBRATE" ] || bad "could not read --recombrate from $(basename "$CMDFILE")"
    note "mutrate     $MUTRATE   (from $(basename "$CMDFILE"), not retyped)"
    note "recombrate  $RECOMBRATE"
  fi
fi

echo
if [ "$fail" -ne 0 ]; then
  echo "Preflight FAILED. Nothing was run."
  exit 1
fi
echo "Preflight passed."
if [ "$CHECK_ONLY" -eq 1 ]; then
  exit 0
fi

# ================================================================== the sweep
mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/sweep_summary.tsv"
if [ ! -s "$SUMMARY" ]; then
  printf 'prefix\tregion\tassumed_N\tmaxtime\tntimes\titers\tsites\thaps\twall_s\tintervals\ttmrca_mean_gen\ttmrca_stdev_gen\tratio_to_N\traw\n' > "$SUMMARY"
fi

# One ceiling for every point in the sweep, so assumed-N is the only variable.
if [ -z "$MAXTIME" ]; then
  biggest=0
  for N in $POPSIZES; do [ "$N" -gt "$biggest" ] && biggest=$N; done
  MAXTIME=$((4 * biggest))
fi
echo "--maxtime $MAXTIME, FIXED for every point in this sweep (largest N x 4)."
echo "Tying it to each run's own N would truncate each prior differently and the"
echo "sweep would measure grid ceiling as well as prior. See the header notes."

REGION_ARG=""
REGION_TAG="whole"
if [ -n "$REGION" ]; then
  REGION_ARG="--region $REGION"
  REGION_TAG="$REGION"
fi

TREE_ARG=""
[ "$EMIT_TREES" -eq 1 ] && TREE_ARG="--tree"

rc=0
for N in $POPSIZES; do
  tag="N${N}"
  [ -n "$REGION" ] && tag="${tag}_r${REGION}"
  out="$OUTDIR/${RUN_PREFIX}_${tag}"

  echo
  echo "=== assumed N=$N   --maxtime $MAXTIME   --ntimes $NTIMES   region $REGION_TAG ==="

  t0=$(date +%s)
  # shellcheck disable=SC2086
  arg-sample --sites "$SITES" \
    --popsize "$N" --mutrate "$MUTRATE" --recombrate "$RECOMBRATE" \
    --maxtime "$MAXTIME" --ntimes "$NTIMES" \
    --iters "$ITERS" --sample-step "$STEP" --randseed "$SEED" \
    $REGION_ARG --output "$out" > "$out.stdout" 2>&1
  if [ $? -ne 0 ]; then
    echo "  arg-sample FAILED — see $out.stdout"; rc=2; continue
  fi
  t1=$(date +%s)
  wall=$((t1 - t0))
  echo "  arg-sample done in ${wall}s — $(ls "$out".*.smc.gz 2>/dev/null | wc -l) sampled ARGs"

  # arg-summarize cannot read .smc; convert to an indexed BED first.
  # smc2bed-all takes the OUTPUT BASE NAME positionally — there is no -o flag —
  # and writes <base>.bed.gz itself.
  if ! smc2bed-all "$out" > "$out.smc2bed.log" 2>&1; then
    echo "  smc2bed-all FAILED — see $out.smc2bed.log"; rc=2; continue
  fi
  [ -r "$out.bed.gz" ] || { echo "  smc2bed-all produced no $out.bed.gz"; rc=2; continue; }

  # Statistics are their own flags in arg-summarize 0.8.1. -s is --subset, not a
  # statistic selector, and passing it as one fails with "need to specify a tree
  # statistic". Without --mean/--stdev/-Q you get one row per MCMC sample.
  raw="$out.tmrca.txt"
  if ! arg-summarize -a "$out.bed.gz" --tmrca --mean --stdev -Q 0.025,0.5,0.975 > "$raw" 2>"$out.summarize.log"; then
    echo "  arg-summarize FAILED — see $out.summarize.log"; rc=2; continue
  fi

  if [ "$EMIT_TREES" -eq 1 ]; then
    # Per-MCMC-sample newick per region. This is the input W7's parser needs to
    # walk back to K lineages — --tmrca alone is TMR-1-A and cannot give TMR-4-A.
    arg-summarize -a "$out.bed.gz" --tree 2>"$out.tree.log" | gzip > "$out.trees.gz" \
      && echo "  newick saved: $(basename "$out").trees.gz"
  fi

  # Span-weight the per-interval depths. §1A: TMR4A is evaluated per local tree
  # and aggregated weighted by span, so a plain mean over intervals would
  # over-weight short ones.
  read -r nint mean sd <<EOF
$(awk -F'\t' '
  # arg-summarize labels its columns; find them by name rather than by position,
  # because which columns exist depends on which statistic flags were passed.
  (NR==1 || /^#/) && /tmrca/ && !seen {
    for (i=1; i<=NF; i++) { h=$i; sub(/^#/,"",h); gsub(/^ +| +$/,"",h); col[h]=i }
    seen=1; next
  }
  /^#/ { next }
  NF < 4 { next }
  {
    mc = ("tmrca_mean" in col) ? col["tmrca_mean"] : 4
    sc = ("tmrca_stdev" in col) ? col["tmrca_stdev"] : 0
    if (NF < mc) next
    s = $2; e = $3; w = e - s; if (w <= 0) w = 1
    tw += w; sm += w * $mc
    if (sc && NF >= sc) ss += w * $sc
    n++
  }
  END { if (tw > 0) printf "%d %.1f %.1f", n, sm/tw, ss/tw; else printf "0 NA NA" }
' "$raw")
EOF
  [ -n "${nint:-}" ] || { nint=0; mean=NA; sd=NA; }

  ratio="NA"
  if [ "$mean" != "NA" ]; then
    ratio=$(awk -v m="$mean" -v n="$N" 'BEGIN{ printf "%.3f", m/n }')
  fi
  echo "  span-weighted TMR-1-A: $mean generations over $nint intervals   (inferred/N = $ratio)"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$RUN_PREFIX" "$REGION_TAG" "$N" "$MAXTIME" "$NTIMES" "$ITERS" \
    "${nsites:-NA}" "${nhap:-NA}" "$wall" "$nint" "$mean" "$sd" "$ratio" \
    "$(basename "$raw")" >> "$SUMMARY"
done

echo
echo "=== Sweep complete ==="
echo "Summary: $SUMMARY"
echo
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo
echo "The RATIO column is the E1 result. Flat near 1 across the sweep means the"
echo "inferred depth is tracking the assumed N rather than the data."
echo
echo "NOTE: tmrca_mean is TMR-1-A (full MRCA), not TMR-4-A. Re-run with -T and"
echo "feed the .trees.gz to W7's parser for the TMR-K-A the study is about."
exit $rc
