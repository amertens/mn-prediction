#!/usr/bin/env bash
# Render a .qmd to .pptx and repair the two defects pandoc leaves behind with
# our reference document. Always use this rather than calling quarto directly:
# a raw `quarto render ... --to pptx` produces a file PowerPoint will not open.
#
#   bash scripts/render_deck.sh docs/slides/MN-proxy-policy-deck-2026-09.qmd
set -euo pipefail
QMD="$1"
PPTX="${QMD%.qmd}.pptx"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin:$PATH"
quarto render "$QMD" --to pptx
python scripts/fix_pptx_dangling_rels.py "$PPTX"
echo "rendered and repaired: $PPTX"
