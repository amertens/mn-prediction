#!/usr/bin/env bash
# Build the illustrated concept slides for a deck: quantities + icons (R), then
# the slides, previews and placement table (python-pptx).
#
#   bash scripts/concept_slides/build.sh docs/slides/MN-proxy-Ghana-concept-slides-2026-09.yaml
#
# The spec names its deck qmd; the deck's setup chunk is evaluated so every number
# on the slides is the one the deck prints. Outputs sit next to the spec:
# <spec>.pptx, <spec>-PLACEMENT.md, concept_previews/slideNN.png.
set -euo pipefail
SPEC="$1"
DECK="$(dirname "$SPEC")/$(grep -m1 '^deck:' "$SPEC" | sed 's/^deck:[[:space:]]*//')"
RSCRIPT="${RSCRIPT:-/c/Program Files/R/R-4.4.2/bin/Rscript.exe}"
"$RSCRIPT" scripts/concept_slides/build_quantities.R "$DECK" "$SPEC"
python scripts/concept_slides/build_concept_slides.py "$SPEC"
