#!/usr/bin/env bash
# Render a .qmd to .pptx and repair the defects pandoc leaves behind with our
# reference documents. Always use this rather than calling quarto directly:
# a raw `quarto render ... --to pptx` produces a file PowerPoint will not open.
#
#   bash scripts/render_deck.sh docs/slides/MN-proxy-policy-deck-2026-09.qmd
#   bash scripts/render_deck.sh docs/slides/MN-proxy-policy-deck-2026-09.qmd --forum
#   bash scripts/render_deck.sh docs/slides/MN-proxy-data-sources-2026-09.qmd --text 18,11
#
# --forum: the deck uses the Micronutrient Forum template; also enlarge the
# figures to the slide and add the unpublished-data label
# (scripts/pptx_forum_postprocess.py).
# --text BODY[,TABLE]: stamp explicit font sizes on body text and tables
# (scripts/pptx_set_text_size.py); the Ghana template's master style is 28 pt,
# too large for text-heavy slides. e.g. --text 18,11
set -euo pipefail
QMD="$1"; shift
PPTX="${QMD%.qmd}.pptx"
FORUM=0; TEXT=""
while [ $# -gt 0 ]; do
  case "$1" in
    --forum) FORUM=1 ;;
    --text) TEXT="$2"; shift ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin:$PATH"
quarto render "$QMD" --to pptx
python scripts/fix_pptx_dangling_rels.py "$PPTX"
if [ "$FORUM" = 1 ]; then
  python scripts/pptx_forum_postprocess.py "$PPTX"
fi
if [ -n "$TEXT" ]; then
  python scripts/pptx_set_text_size.py "$PPTX" --body "${TEXT%%,*}" --table "${TEXT#*,}"
fi
echo "rendered and repaired: $PPTX"
