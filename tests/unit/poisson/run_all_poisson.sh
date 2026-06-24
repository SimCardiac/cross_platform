#!/bin/bash
set -e
BINDIR="${1:-.}"

SOLVERS=("poisson_vertex_2D" "poisson_DMStag_2D" "poisson_staggered_2D")
LABELS=("Vertex (DMDA)" "Cell (DMStag)" "Staggered (Mixed)")

# MMS cases: name, expected rate with compatible BC
MMS_CASES=("sinpi:2.0" "poly2:2.0" "cospi:2.0")

# BC: label left right bottom top (same as solver naming)
BCS=("D-D-D-D" "N-N-N-N" "D-D-N-N" "N-N-D-D")
BC_LABELS=("All Dirichlet" "All Neumann" "D-LR / N-TB" "N-LR / D-TB")

GRID=(16 32 64 128)

rate(){ python3 -c "import math;print(f'{math.log(float($1)/float($2))/math.log(2):.3f}')" 2>/dev/null||echo " --- "; }

echo "============================================================"
echo " Poisson — MMS × BC × Grid Convergence Matrix"
echo " $(date)"
echo "============================================================"
echo ""

for si in "${!SOLVERS[@]}"; do
  s="${SOLVERS[$si]}"; lb="${LABELS[$si]}"
  [ ! -f "$BINDIR/$s" ] && { echo "=== $lb: SKIP (no $s) ==="; echo; continue; }
  echo "======================================================================"
  echo "  $lb  ($s)"
  echo "======================================================================"

  for mm in "${MMS_CASES[@]}"; do
    mn="${mm%%:*}"; mr="${mm##*:}"
    echo "  MMS: $mn (expected rate ~$mr)"
    echo ""

    for bi in "${!BCS[@]}"; do
      bc="${BCS[$bi]}"; bl="${BC_LABELS[$bi]}"
      printf "  BC: %-20s" "$bl ($bc)"
      printf " | %-10s %-16s %s\n" "Grid" "Error" "Rate"
      echo "  $(printf '%.0s-' {1..55})"
      pe="---"
      for n in "${GRID[@]}"; do
        pn="${n}x${n}"
        out=$("$BINDIR/$s" -nx "$n" -ny "$n" -mms "$mn" -poisson_check_error -convergence_test 2>&1)
        cl=$(echo "$out"|grep "CONVERGENCE:"||true)
        if [ -n "$cl" ]; then
          e=$(echo "$cl"|awk '{print $4}')
          r=$(rate "$pe" "$e")
          printf "  %-22s | %-10s %-16s %s\n" "" "$pn" "$e" "$r"
          pe="$e"
        else
          printf "  %-22s | %-10s FAILED\n" "" "$pn"
        fi
      done
      echo ""
    done
    echo "  --------------------------------------------------"
    echo ""
  done
done
echo "Done."
