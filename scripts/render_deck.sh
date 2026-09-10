#!/usr/bin/env bash
# Render a .qmd to .pptx and repair the defects pandoc leaves behind with our
# reference documents. Always use this rather than calling quarto directly:
# a raw `quarto render ... --to pptx` produces a file PowerPoint will not open.
#
#   bash scripts/render_deck.sh docs/slides/MN-proxy-policy-deck-2026-09.qmd
#   bash scripts/render_deck.sh docs/slides/MN-proxy-policy-deck-2026-09.qmd --forum
#
# --forum: the deck uses the Micronutrient Forum template; also enlarge the
# figures to the slide and add the unpublished-data label
# (scripts/pptx_forum_postprocess.py).
set -euo pipefail
QMD="$1"
PPTX="${QMD%.qmd}.pptx"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin:$PATH"
quarto render "$QMD" --to pptx
python scripts/fix_pptx_dangling_rels.py "$PPTX"
if [ "${2:-}" = "--forum" ]; then
  python scripts/pptx_forum_postprocess.py "$PPTX"
fi
echo "rendered and repaired: $PPTX"
