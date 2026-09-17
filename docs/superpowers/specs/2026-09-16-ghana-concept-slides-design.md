# Ghana deck: illustrated concept slides

Date: 2026-09-16. Status: approved in conversation.

## Goal

Give the concept slides of `docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd`
the icon-illustrated look of the DHS-style decks (icon + heading + caption rows,
cards, a rounded title box), on the existing Ghana template. The figure and
result-table slides are unchanged. The user assembles the final deck by hand:
the deliverable is a standalone pptx of replacement slides plus a placement
table, not a change to the Quarto render.

## Deliverables

1. `docs/slides/MN-proxy-Ghana-concept-slides-2026-09.pptx`: ten slides on the
   Ghana reference template (13.33 x 7.5 in). Each slide's speaker notes name
   the deck slide it replaces and carry that slide's existing notes.
2. `docs/slides/concept_previews/*.png`: one preview per slide (LibreOffice
   PDF export, then rasterised).
3. `docs/slides/MN-proxy-Ghana-concept-slides-PLACEMENT.md`: table of
   new slide -> deck slide number and title -> what to copy across.

## Slides

| # | Replaces (deck slide) | Layout |
|---|---|---|
| 1 | 3 Motivation | icon rows (3), photo right |
| 2 | 4 Aims | icon rows (3) |
| 3 | 7 Outcomes | five nutrient cards: icon, population, biomarker, cut-off |
| 4 | 9 Data assembly | six-step pipeline of icon cards |
| 5 | 11 Predictor domains | domain grid: icon, domain, example layers, grouped |
| 6 | 13 Estimators | two panels, January vs current, icons |
| 7 | 14 Estimands | three cards: in-fill, region, country; fold and baseline |
| 8 | 43 Conclusions | icon rows (4), headline number large |
| 9 | 44 What the models can and cannot do yet | checklist tiles (yes / partly / no) |
| 10 | 45 Next steps | two panels, planned vs open questions |

## Visual language

- Title: template title placeholder, unchanged (so the deck stays consistent
  when the slides are pasted in); the rounded title box from the reference
  deck is NOT used, the template's own title style is.
- Icons: Font Awesome 6 Free (solid), rendered by `fontawesome::fa_png()` in the
  template blue `#1f4e79` at 256 px, committed as PNGs under
  `docs/slides/img/icons/`. Status tiles use green / amber / grey.
- Text: Calibri (template default); heading 20 pt bold in `#1f4e79`, caption
  14-16 pt in dark grey; card fill `#eef3f8`, no outline; 0.15 in corner radius.
- Numbers: every figure in the text is read from the result tables through the
  qmd's own setup chunk (see Mechanics), never typed.

## Mechanics

- `scripts/concept_slides/build_quantities.R`: `knitr::purl()` the deck qmd,
  evaluate its `setup` chunk with the working directory set to `docs/slides/`,
  write the `Q` list to `docs/slides/MN-proxy-Ghana-presentation-2026-09.quantities.json`,
  and render the icon PNGs the spec names (idempotent; skips existing files).
- `docs/slides/MN-proxy-Ghana-concept-slides-2026-09.yaml`: one entry per slide
  with `replaces`, `title`, `layout` (`icon_rows`, `cards`, `pipeline`, `grid`,
  `two_panel`, `checklist`) and `items`; text may contain `{name}` placeholders
  resolved from the quantities JSON (`{name:.2f}`, `{name:pc}` formats).
- `scripts/concept_slides/build_concept_slides.py`: python-pptx; opens the
  Ghana reference deck, removes its slides, adds one "Title Only" slide per
  spec entry, draws the layout, writes the notes, saves the pptx; then exports
  previews through LibreOffice and writes the placement table.
- One command: `bash scripts/concept_slides/build.sh`.

## Out of scope

Changing the Quarto render, the other decks, or the figure slides. If the
look is adopted for good, a later step can fold the YAML + Python into
`render_deck.sh` (populate-in-place by title), which this design leaves open.

## First milestone

Three slides first (Data assembly, Estimands, What the models can and cannot
do), previews sent for a look before the remaining seven are written.

## Addendum (2026-09-16, same day): the full 30-45 minute talk

Built on the same machinery, now wired into the render:

- `docs/slides/MN-proxy-full-talk-2026-09.qmd`: the data-sources deck and the
  Ghana analysis deck combined (63 slides + section headers; appendix holds
  what the speaker cuts). Setup chunk = Ghana chunk + data-sources additions,
  one `Q` list. New slides: big-data/small-data framing, the four surveys,
  how the zero-tuning index works, Malawi observed-vs-predicted maps, the
  surrogate-marker maps (women's B12 vs the IHME anaemia surface vs fish
  consumption, one percentile scale), survey augmentation (what is supported,
  in design, not shown), the dashboard (live screenshot).
- `docs/slides/MN-proxy-full-talk-2026-09.concept.yaml`: 19 concept slides
  drawn IN PLACE by `render_deck.sh --concept` (builder `--into` mode: match
  by title, skip section headers, keep title/footer placeholders, drop the
  empty body); `fill_figures: true` enlarges lone figures to the content area.
- Numbers: qmd prose, notes and YAML all read `Q`; one command refreshes all:
  `bash scripts/render_deck.sh docs/slides/MN-proxy-full-talk-2026-09.qmd --text 18,10 --concept docs/slides/MN-proxy-full-talk-2026-09.concept.yaml`
- `scripts/concept_slides/dashboard_screenshots.py` (playwright) refreshes
  `docs/slides/img/dashboard_*.png`.
