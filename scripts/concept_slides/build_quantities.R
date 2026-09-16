# Dump the quantities a deck's setup chunk computes (its `Q` list) to JSON, and
# render the Font Awesome icons a concept-slide spec names, so that
# scripts/concept_slides/build_concept_slides.py can build the illustrated
# slides with the same numbers the deck prints.
#
#   Rscript scripts/concept_slides/build_quantities.R docs/slides/MN-proxy-Ghana-presentation-2026-09.qmd \
#           docs/slides/MN-proxy-Ghana-concept-slides-2026-09.yaml
#
# Writes <deck>.quantities.json next to the qmd and PNGs to docs/slides/img/icons/
# (skipped when present). The setup chunk is evaluated as the deck would: with
# the working directory at docs/slides/ and `root` resolved by here::here().
suppressPackageStartupMessages({library(jsonlite); library(fontawesome); library(yaml)})
args <- commandArgs(trailingOnly = TRUE)
qmd <- normalizePath(args[1]); spec <- normalizePath(args[2])
root <- normalizePath(file.path(dirname(qmd), "..", ".."))
ICON_DIR <- file.path(root, "docs", "slides", "img", "icons"); dir.create(ICON_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. evaluate the setup chunk -------------------------------------------------------
lines <- readLines(qmd, warn = FALSE)
open <- grep("^```[{]r setup", lines); close <- grep("^```[[:space:]]*$", lines); close <- close[close > open[1]][1]
code <- lines[(open[1] + 1):(close - 1)]
env <- new.env(); assign("root", root, envir = env)
owd <- setwd(dirname(qmd)); on.exit(setwd(owd), add = TRUE)
suppressPackageStartupMessages(eval(parse(text = code), envir = env))
Q <- get("Q", envir = env)
# vectors named d/b/n (paired comparisons) become objects; data frames become row lists
Q <- lapply(Q, function(x) if (is.data.frame(x)) x else if (!is.null(names(x))) as.list(x) else x)

# ---- 2. extra quantities the concept slides use ---------------------------------------
MD <- get("MD", envir = env); SY <- get("SY", envir = env); TG <- get("TG", envir = env)
if (!is.null(MD)) {
  dc <- as.data.frame(table(MD$domain), stringsAsFactors = FALSE); Q$domain_counts <- setNames(as.list(dc$Freq), dc$Var1)
  src <- split(MD$source, MD$domain); Q$domain_sources <- lapply(src, function(s) sort(unique(sub(" [(].*$", "", s))))
}
if (!is.null(SY) && !is.null(TG)) {
  sy <- SY[SY$in_protocol %in% c(TRUE, "TRUE"), ]
  dist <- aggregate(Admin2 ~ country, unique(TG[, c("country", "Admin1", "Admin2")]), length)
  Q$surveys <- lapply(seq_len(nrow(sy)), function(i) list(country = sy$country[i], survey = sy$survey[i], year = sy$survey_year[i],
    districts = dist$Admin2[match(sy$country[i], dist$country)]))
}
out <- sub("[.]qmd$", ".quantities.json", qmd)
write_json(Q, out, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
cat("wrote", out, "with", length(Q), "quantities\n")

# ---- 3. icons ---------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a
sp <- yaml::read_yaml(spec)
want <- unique(unlist(lapply(sp$slides, function(s) {
  it <- c(s$items, unlist(lapply(s$panels, `[[`, "items"), recursive = FALSE), unlist(lapply(s$groups, `[[`, "items"), recursive = FALSE))
  st <- c(yes = "check@green", partly = "triangle-exclamation@amber", no = "xmark@grey")
  c(sapply(it, function(x) x$icon %||% NULL), sapply(it, function(x) if (is.null(x$status)) NULL else st[[x$status]]),
    sapply(s$panels, function(p) if (is.null(p$icon)) NULL else if (grepl("@", p$icon)) p$icon else paste0(p$icon, "@white")), sapply(s$groups, function(g) g$icon %||% NULL))
})))
cols <- c(blue = "#1f4e79", light = "#6baed6", green = "#2e7d32", amber = "#e0a100", grey = "#8a8a8a", white = "#ffffff")
n <- 0
for (w in want) {   # "name" or "name@colour"
  parts <- strsplit(w, "@", fixed = TRUE)[[1]]; nm <- parts[1]; col <- if (length(parts) > 1) parts[2] else "blue"
  f <- file.path(ICON_DIR, paste0(nm, "_", col, ".png"))
  if (file.exists(f)) next
  ok <- tryCatch({ fontawesome::fa_png(nm, file = f, fill = cols[[col]], height = 256); TRUE }, error = function(e) { message("no icon: ", nm, " (", conditionMessage(e), ")"); FALSE })
  if (ok) n <- n + 1
}
cat("rendered", n, "new icon(s) into", ICON_DIR, "\n")
