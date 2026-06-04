# =============================================================================
# scripts/fix_quarto_chunks.R
#
# Convert mixed old-knitr / Quarto chunk-option syntax in docs/manuscript_mcn.qmd
# to pure Quarto YAML chunk options. Examples:
#
#   ```{r tbl-cohort, tbl-cap="Sample sizes ..."}
# → ```{r tbl-cohort}
#   #| tbl-cap: "Sample sizes ..."
#
#   ```{r fig-all-methods-iron, fig.cap="...", fig.width=11, fig.height=6}
# → ```{r fig-all-methods-iron}
#   #| fig-cap: "..."
#   #| fig-width: 11
#   #| fig-height: 6
# =============================================================================
qmd <- here::here("docs", "manuscript_mcn.qmd")
src <- readLines(qmd, warn = FALSE)

# Find chunk headers with multiple options.
hdr_re <- "^```\\{r ([a-zA-Z0-9_.-]+), (.+)\\}$"
matches <- grep(hdr_re, src)
cat(sprintf("matched %d chunk headers\n", length(matches)))

new_src <- src
shifted <- 0L
for (i in matches) {
  line_now <- new_src[i + shifted]
  m <- regmatches(line_now, regexec(hdr_re, line_now))[[1]]
  if (length(m) != 3) next
  label <- m[2]
  opts  <- m[3]
  # Split options carefully respecting quoted strings.
  parts <- character()
  cur <- ""; in_quote <- FALSE
  for (ch in strsplit(opts, "", fixed = TRUE)[[1]]) {
    if (ch == "\"") in_quote <- !in_quote
    if (ch == "," && !in_quote) {
      parts <- c(parts, trimws(cur)); cur <- ""
    } else cur <- paste0(cur, ch)
  }
  parts <- c(parts, trimws(cur))
  # Convert each option to YAML.
  yaml_lines <- character()
  for (p in parts) {
    eq <- regexpr("=", p, fixed = TRUE)
    if (eq == -1) next
    key <- trimws(substr(p, 1, eq - 1))
    val <- trimws(substr(p, eq + 1, nchar(p)))
    # Translate knitr keys to Quarto keys.
    key_q <- switch(key,
                     "fig.cap"    = "fig-cap",
                     "fig.width"  = "fig-width",
                     "fig.height" = "fig-height",
                     "fig.dpi"    = "fig-dpi",
                     "tbl.cap"    = "tbl-cap",
                     "out.width"  = "out-width",
                     "out.height" = "out-height",
                     key)
    yaml_lines <- c(yaml_lines, sprintf("#| %s: %s", key_q, val))
  }
  new_hdr <- sprintf("```{r %s}", label)
  new_block <- c(new_hdr, yaml_lines)
  # Replace line in new_src at position i + shifted with new_block.
  pos <- i + shifted
  new_src <- c(new_src[seq_len(pos - 1)], new_block,
                if (pos < length(new_src)) new_src[(pos + 1):length(new_src)]
                else character())
  shifted <- shifted + (length(new_block) - 1)
}

writeLines(new_src, qmd)
cat(sprintf("rewrote %s (%d lines -> %d)\n", qmd, length(src), length(new_src)))
