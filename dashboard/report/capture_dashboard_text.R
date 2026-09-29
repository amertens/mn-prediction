# =============================================================================
# dashboard/report/capture_dashboard_text.R
#
# Visits every page and sub-tab of a running dashboard with headless Chrome,
# expands folded sections, takes a full-page screenshot, and records the
# visible text in reading order (headings, paragraphs, lists, notes, controls,
# tables; charts and maps are marked, not transcribed). Writes one folder with
# PNGs and manifest.json, which build_text_review_docx.js turns into a Word
# document for reviewing the writing.
#
#   Rscript dashboard/report/capture_dashboard_text.R <out_dir> [url]
#   (start the app first: shiny::runApp("dashboard", port = 7791))
# =============================================================================
suppressPackageStartupMessages({library(chromote); library(jsonlite)})
args <- commandArgs(TRUE)
OUT <- if (length(args)) args[1] else "dashboard/report/out/text_capture"
URL <- if (length(args) > 1) args[2] else "http://127.0.0.1:7791"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

b <- ChromoteSession$new(width = 1400, height = 900)
on.exit(try(b$close(), silent = TRUE))
for (attempt in 1:4) {   # the first request after runApp can be slow
  ok <- tryCatch({ b$Page$navigate(URL, timeout_ = 60); TRUE }, error = function(e) FALSE)
  if (ok) break
  Sys.sleep(10)
}
Sys.sleep(15)
js <- function(code) b$Runtime$evaluate(code, returnByValue = TRUE, timeout_ = 120)$result$value

# The walker: visible text in DOM order, with a type per block.
WALK <- "
function(root){
  const out = [];
  const vis = el => {
    const cs = getComputedStyle(el);
    if (cs.display === 'none' || cs.visibility === 'hidden') return false;
    if (cs.display === 'contents') return true;        // Shiny output containers inside bslib layouts
    return !!(el.offsetWidth || el.offsetHeight || el.getClientRects().length);
  };
  const txt = el => (el.innerText || '').replace(/\\s+/g, ' ').trim();
  const has = (el, c) => el.classList && el.classList.contains(c);
  function walk(el){
    if (!vis(el)) return;
    const tag = el.tagName.toLowerCase();
    if (['script','style','svg','noscript','canvas'].includes(tag)) return;
    if (has(el,'leaflet-container')) { out.push({t:'widget', s:'[Interactive map]'}); return; }
    if (has(el,'plotly') || (has(el,'html-widget') && el.querySelector('.plotly'))) { out.push({t:'widget', s:'[Interactive chart]'}); return; }
    if (has(el,'reactable')) {
      const head = Array.from(el.querySelectorAll('.rt-thead .rt-tr:last-child .rt-th')).map(txt);
      const rows = Array.from(el.querySelectorAll('.rt-tbody .rt-tr')).slice(0, 12).map(r => Array.from(r.querySelectorAll('.rt-td')).map(txt));
      out.push({t:'rtable', head: head, rows: rows}); return; }
    if (has(el,'selectize-dropdown') || has(el,'irs')) return;
    if (has(el,'card-header') && el.querySelector('.nav')) {   // a row of sub-tabs, not a heading
      const tabs = Array.from(el.querySelectorAll('.nav-link')).map(txt).filter(Boolean);
      const act = el.querySelector('.nav-link.active');
      out.push({t:'ctrl', s:'Tabs on this page', opts: tabs, sel: act ? txt(act) : ''}); return; }
    if (/^h[1-6]$/.test(tag) || has(el,'card-header') || tag === 'summary' || has(el,'accordion-button')) { const s = txt(el); if (s) out.push({t:'h', s: s}); return; }
    if (has(el,'methods-note')) { out.push({t:'note', s: txt(el)}); return; }
    if (has(el,'alert')) { out.push({t:'note', s: txt(el)}); return; }
    if (tag === 'table') {
      const rows = Array.from(el.querySelectorAll('tr')).map(r => Array.from(r.querySelectorAll('th,td')).map(txt));
      out.push({t:'table', head: rows[0] || [], rows: rows.slice(1)}); return; }
    if (tag === 'li') { out.push({t:'li', s: txt(el)}); return; }
    if (tag === 'dt') { out.push({t:'dt', s: txt(el)}); return; }
    if (tag === 'p' || tag === 'dd' || tag === 'pre') { const s = txt(el); if (s) out.push({t:'p', s: s}); return; }
    if (has(el,'shiny-input-container')) {
      const lab = el.querySelector('.control-label');
      const opts = Array.from(el.querySelectorAll('.radio label, .checkbox label')).map(txt).filter(Boolean);
      const sel = el.querySelector('.selectize-input .item');
      const btn = el.querySelector('.btn');
      out.push({t:'ctrl', s: lab ? txt(lab) : (btn ? txt(btn) : ''), opts: opts, sel: sel ? txt(sel) : ''}); return; }
    if (tag === 'button' || (tag === 'a' && has(el,'btn'))) { const s = txt(el); if (s) out.push({t:'btn', s: s}); return; }
    let direct = '';
    el.childNodes.forEach(n => { if (n.nodeType === 3) direct += n.textContent; });
    if (direct.replace(/\\s+/g, '').length) { const s = txt(el); if (s) out.push({t:'p', s: s}); return; }
    Array.from(el.children).forEach(walk);
  }
  walk(root);
  return JSON.stringify(out);
}"

extract_pane <- function() {
  s <- js(sprintf("(function(){ const W = %s;
    const panes = Array.from(document.querySelectorAll('.tab-content > .tab-pane.active'));
    const top = panes.find(p => !p.parentElement.closest('.tab-pane'));
    return top ? W(top) : '[]'; })()", WALK))
  fromJSON(s, simplifyVector = FALSE)
}
extract_sel <- function(selector) {
  s <- js(sprintf("(function(){ const W = %s; const el = document.querySelector(%s); return el ? W(el) : '[]'; })()",
                  WALK, toJSON(selector, auto_unbox = TRUE)))
  fromJSON(s, simplifyVector = FALSE)
}
click_nav <- function(label, wait = 8) {
  r <- js(sprintf("(function(){
    const l = Array.from(document.querySelectorAll('a.nav-link, a.dropdown-item')).find(a => a.textContent.trim().startsWith(%s));
    if (!l) return 'missing';
    const d = l.closest('.nav-item.dropdown'); if (d) d.querySelector('a.dropdown-toggle').click();
    l.click(); return 'ok'; })()", toJSON(label, auto_unbox = TRUE)))
  if (!identical(r, "ok")) warning("nav not found: ", label)
  Sys.sleep(wait)
  js("document.body.click(); 1")   # close any open dropdown
}
open_all <- function(wait = 3) {
  js("document.querySelectorAll('.tab-pane.active .accordion-button.collapsed').forEach(b => b.click());
      document.querySelectorAll('.tab-pane.active details').forEach(d => d.open = true); 1")
  Sys.sleep(wait)
}
shot <- function(file) {
  h <- js("Math.max(document.documentElement.scrollHeight, document.body.scrollHeight)")
  res <- b$Page$captureScreenshot(format = "png", captureBeyondViewport = TRUE,
                                  clip = list(x = 0, y = 0, width = 1400, height = h, scale = 1), timeout_ = 120)
  writeBin(base64_dec(res$data), file.path(OUT, file))
  file
}
set_select <- function(id, value, wait = 5) {
  js(sprintf("(function(){ const el = $('#%s')[0]; if (el && el.selectize) { el.selectize.setValue(%s); return 1; } return 0; })()",
             id, toJSON(value, auto_unbox = TRUE)))
  Sys.sleep(wait)
}

manifest <- list(banner = js("(document.querySelector('.alert-info') || {innerText:''}).innerText.replace(/\\s+/g,' ').trim()"),
                 footer = js("(Array.from(document.querySelectorAll('div')).reverse().find(d => /data built/.test(d.textContent) && d.children.length === 0) || {innerText:''}).innerText.trim()"),
                 captured = format(Sys.time(), "%Y-%m-%d %H:%M"), url = URL, pages = list())
add_page <- function(id, title, blocks, shots, extra = NULL) {
  manifest$pages[[length(manifest$pages) + 1]] <<- list(id = id, title = title, shots = as.list(shots), blocks = blocks, extra = extra)
  cat(sprintf("  captured %-28s %3d blocks\n", id, length(blocks)))
}

# ── pages ───────────────────────────────────────────────────────────────────
PAGES <- getOption("capture.pages", list(
  list(id = "start", title = "Start here", nav = "Start here", open = TRUE),
  list(id = "map", title = "Map explorer", nav = "Map explorer", map = TRUE),
  list(id = "district", title = "District profiles", nav = "District profiles"),
  list(id = "civ", title = "Cote d'Ivoire", nav = "Cote d'Ivoire", open = TRUE),
  list(id = "imp_top", title = "What drives the estimate: leading data layers", nav = "What drives the estimate", sub = "Leading data layers"),
  list(id = "imp_search", title = "What drives the estimate: search all layers", sub = "Search all layers"),
  list(id = "imp_groups", title = "What drives the estimate: which data groups matter", sub = "Which data groups matter"),
  list(id = "imp_twenty", title = "What drives the estimate: twenty public layers", sub = "Twenty public layers"),
  list(id = "imp_pattern", title = "What drives the estimate: a shared pattern", sub = "A shared pattern"),
  list(id = "nut_by", title = "The data behind it: what tracks which nutrient", nav = "The data behind it", sub = "What tracks which nutrient"),
  list(id = "nut_fam", title = "The data behind it: indicator families", sub = "Indicator families"),
  list(id = "nut_ind", title = "The data behind it: individual indicators", sub = "Individual indicators"),
  list(id = "catalogue", title = "The data behind it: data layer catalogue", sub = "Data layer catalogue", catalogue = TRUE),
  list(id = "trust_tests", title = "How well it works: three tests", nav = "How well it works", sub = "Three tests"),
  list(id = "trust_ceiling", title = "How well it works: best achievable score", sub = "Best achievable score"),
  list(id = "trust_curve", title = "How well it works: each survey added", sub = "Each survey added"),
  list(id = "trust_geo", title = "How well it works: compared with a geostatistical model", sub = "Compared with a geostatistical model"),
  list(id = "trust_tried", title = "How well it works: what else was tried", sub = "What else was tried"),
  list(id = "external", title = "Tested in six more countries", nav = "Tested in six more countries"),
  list(id = "target_burden", title = "What the ranking buys: burden reached", nav = "What the ranking buys", sub = "Burden reached"),
  list(id = "target_signal", title = "What the ranking buys: how reliable is the worst-fifth signal?", sub = "How reliable is the worst-fifth signal?"),
  list(id = "target_who", title = "What the ranking buys: WHO severity bands", sub = "WHO severity bands"),
  list(id = "roadmap", title = "What more data buys", nav = "What more data buys"),
  list(id = "methods", title = "Methods", nav = "Methods", open = TRUE),
  list(id = "plan_where", title = "Plan a survey: where to sample next", nav = "Plan a survey", sub = "Where to sample next"),
  list(id = "plan_size", title = "Plan a survey: how small can a survey be", sub = "How small can a survey be"),
  list(id = "plan_choose", title = "Plan a survey: does choosing districts well matter?", sub = "Does choosing districts well matter?")
))
if (!identical(Sys.getenv("CAPTURE_VARIANT"), "full")) {
  # the current dashboard (concise text since 2026-09-27): no "What else was tried" sub-tab, and
  # Technical notes replaces Methods. CAPTURE_VARIANT=full captures the archived full-text layout.
  PAGES <- Filter(function(p) !p$id %in% c("trust_tried", "methods"), PAGES)
  PAGES <- c(PAGES, list(list(id = "technical", title = "Technical notes (appendix page)", nav = "Technical notes", open = TRUE)))
}

for (pg in PAGES) {
  if (!is.null(pg$nav)) click_nav(pg$nav, 10)
  if (!is.null(pg$sub)) click_nav(pg$sub, 8)
  if (isTRUE(pg$open)) open_all()
  extra <- NULL
  if (isTRUE(pg$catalogue)) {
    js("(function(){ const r = document.querySelector('#catalogue-table .rt-tbody .rt-tr'); if (r) r.click(); return 1; })()")
    Sys.sleep(6)
  }
  if (isTRUE(pg$map)) {
    # a clicked district shows the click-through panel
    js("Shiny.setInputValue('map-map_shape_click', {id: 'Tamale', '.nonce': Math.random()}, {priority: 'event'}); 1")
    Sys.sleep(5)
  }
  f <- shot(sprintf("%02d_%s.png", length(manifest$pages) + 1, pg$id))
  blocks <- extract_pane()
  shots <- f
  if (isTRUE(pg$map)) {
    # every layer's caption, then the cross-hatched rank-range view
    caps <- list()
    for (ly in c("priority", "p_worst_fifth", "rank_width", "prev_anchored", "p_modplus_cal", "survey_prev", "who_class", "people_affected")) {
      set_select("map-layer", ly, 4)
      caps[[length(caps) + 1]] <- list(layer = js("(function(){ const el = $('#map-layer')[0]; return el && el.selectize ? el.selectize.getItem(el.selectize.getValue()).text() : ''; })()"),
                                       caption = js("(document.getElementById('map-caption') || {innerText:''}).innerText.trim()"))
    }
    set_select("map-layer", "rank_width", 3)
    js("(function(){ const c = document.getElementById('map-hatch_unstable'); if (c && !c.checked) { c.click(); } return 1; })()")
    Sys.sleep(6)
    shots <- c(shots, shot(sprintf("%02d_%s_hatch.png", length(manifest$pages) + 1, pg$id)))
    js("(function(){ const c = document.getElementById('map-hatch_unstable'); if (c && c.checked) { c.click(); } return 1; })()")
    set_select("map-layer", "priority", 3)
    extra <- list(captions = caps)
  }
  add_page(pg$id, pg$title, blocks, shots, extra)
}

# the two popovers
for (pp in list(list(id = "glossary", title = "Glossary (pop-up)", btn = "Glossary"),
                list(id = "about", title = "About (pop-up)", btn = "About"))) {
  js(sprintf("(function(){ const b = Array.from(document.querySelectorAll('button')).find(x => x.textContent.trim() === %s); if (b) b.click(); return 1; })()",
             toJSON(pp$btn, auto_unbox = TRUE)))
  Sys.sleep(3)
  blocks <- extract_sel(".popover .popover-body")
  f <- shot(sprintf("%02d_%s.png", length(manifest$pages) + 1, pp$id))
  add_page(pp$id, pp$title, blocks, f)
  js(sprintf("(function(){ const b = Array.from(document.querySelectorAll('button')).find(x => x.textContent.trim() === %s); if (b) b.click(); return 1; })()",
             toJSON(pp$btn, auto_unbox = TRUE)))
  Sys.sleep(2)
}

writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(OUT, "manifest.json"), useBytes = TRUE)
cat(sprintf("wrote %s: %d pages\n", file.path(OUT, "manifest.json"), length(manifest$pages)))
