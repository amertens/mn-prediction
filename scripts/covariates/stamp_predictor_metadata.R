# =============================================================================
# scripts/covariates/stamp_predictor_metadata.R   [TP-01 / TA-01, 2026-09-15]
#
# STAMP THE SHARED PREDICTOR METADATA WITH THE THREE THINGS THE POLICY NEEDS
#
# The block builders write `column, domain, source, n_countries, countries,
# completeness, subnational`. This step, run after the last block is appended
# and the exclusion policy applied (rebuild order: builder -> 07 -> 08 -> 59 ->
# [60] -> 62 -> 53 -> STAMP -> audit), adds and reconciles:
#
#   subnational        recomputed from the live data (varies within at least
#                      one country: >= 4 finite values and sd > 0). A block
#                      that declared TRUE for a column that is constant
#                      everywhere is corrected and reported.
#   tier               open / survey_public / survey_dhs, from the column's
#                      source by metadata/covariates/predictor_tiers.csv
#                      (assign_tier_v2 in R/protocol_v2.R)
#   domain             regrouped by metadata/covariates/domain_overrides.csv
#                      (the block's own label is kept in domain_block)
#   modelled_surface   TRUE for IHME / MIMI model outputs (V2_DROP_MODELLED)
#   alignment_rule     the row of metadata/covariates/temporal_alignment.csv
#   year_used          the data year behind the column, per country where it
#                      differs ("Gambia=2019;Ghana=2014;..."), one number
#                      where it does not, "static" for time-invariant layers
#   year_offset_max_abs  max |year_used - survey_year| over the countries in
#                      which the column has data (NA for static layers)
#
# Every column must match an alignment rule; an unmatched column stops the
# script, so a new block cannot enter the set without declaring its year.
#
#   Rscript -e "source('scripts/covariates/stamp_predictor_metadata.R')"
# -> data/covariates/harmonized/predictors_admin2_shared_metadata.csv (in place)
# -> results/tables/predictor_alignment_<date>.csv   (column x country years)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R"); source("R/survey_years.R")
HDIR <- "data/covariates/harmonized"; STAMP <- format(Sys.Date(), "%Y-%m-%d")
SURVEY_YEAR <- survey_years()

S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
M <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
cols <- setdiff(names(S), c("country", "Admin1", "Admin2"))
stopifnot(setequal(cols, M$column), !anyDuplicated(M$column))
M <- M[match(cols, M$column), ]
COUNTRIES <- names(SURVEY_YEAR); stopifnot(all(COUNTRIES %in% S$country))

# ── 1. subnational, from the data ────────────────────────────────────────────
has_data <- sapply(COUNTRIES, function(cc) vapply(cols, function(v) sum(is.finite(S[[v]][S$country == cc])) >= 4, NA))
varies   <- sapply(COUNTRIES, function(cc) vapply(cols, function(v) { z <- S[[v]][S$country == cc]; z <- z[is.finite(z)]; length(z) >= 4 && stats::sd(z) > 0 }, NA))
sub_new <- rowSums(varies) > 0
sub_old <- .v2_as_logical(M$subnational)
chg <- which(is.na(sub_old) | sub_old != sub_new)
if (length(chg)) cat(sprintf("[stamp] subnational corrected for %d column(s): %s\n", length(chg), paste(sprintf("%s (%s -> %s)", cols[chg], sub_old[chg], sub_new[chg]), collapse = ", ")))
M$subnational <- sub_new

# ── 2. tier, domain overrides, modelled-surface flag ─────────────────────────
M$tier <- assign_tier_v2(M)
# Domains: the block builders assign prefix-driven labels (86 columns landed in
# one "Infant and child morbidity/mortality" catch-all). metadata/covariates/
# domain_overrides.csv regroups by column regex, first match wins; the block's
# own label is kept in domain_block so nothing is lost.
if (!"domain_block" %in% names(M)) M$domain_block <- M$domain
M$domain_block <- ifelse(is.na(M$domain_block) | !nzchar(M$domain_block), M$domain, M$domain_block)   # rows appended since the last stamp carry only `domain`
DO <- read.csv("metadata/covariates/domain_overrides.csv", stringsAsFactors = FALSE)
dom <- M$domain_block; hit_any <- rep(FALSE, nrow(M))
for (i in seq_len(nrow(DO))) { h <- !hit_any & grepl(DO$column_regex[i], M$column, perl = TRUE); dom[h] <- DO$domain[i]; hit_any <- hit_any | h }
chg <- sum(dom != M$domain); M$domain <- dom
cat(sprintf("[stamp] domain overrides: %d rule(s), %d column(s) regrouped; %d domains\n", nrow(DO), chg, length(unique(M$domain))))
# Modelled surfaces (someone else's model output, near the outcome): flagged by
# SOURCE, so the V2_DROP_MODELLED sensitivity no longer depends on a domain label.
M$modelled_surface <- grepl("IHME|Tang et al|MODELLED", M$source, ignore.case = TRUE) | grepl("MODELLED SURFACE", M$domain_block, fixed = TRUE)

# ── 3. temporal alignment ────────────────────────────────────────────────────
R <- read.csv("metadata/covariates/temporal_alignment.csv", stringsAsFactors = FALSE)
stopifnot(all(c("rule_id", "column_regex", "year_rule") %in% names(R)), !anyDuplicated(R$rule_id))
rule_of <- rep(NA_integer_, length(cols))
for (i in seq_len(nrow(R))) { hit <- is.na(rule_of) & grepl(R$column_regex[i], cols, perl = TRUE); rule_of[hit] <- i }
if (anyNA(rule_of)) stop("[stamp] no alignment rule for: ", paste(cols[is.na(rule_of)], collapse = ", "), " - add a row to metadata/covariates/temporal_alignment.csv")

year_for <- function(v, rule, country) {
  sy <- SURVEY_YEAR[[country]]; kind <- sub(":.*$", "", rule); arg <- sub("^[^:]*:?", "", rule)
  switch(kind,
    survey_year = sy,
    static = NA_integer_,
    fixed = as.integer(arg),
    from_name = as.integer(regmatches(v, regexpr("(19|20)[0-9]{2}", v))),
    per_country = { kv <- strsplit(strsplit(arg, ";")[[1]], "="); y <- stats::setNames(as.integer(vapply(kv, `[`, "", 2)), vapply(kv, `[`, "", 1))
                    if (is.na(y[country])) stop("per_country rule lacks ", country, ": ", rule) else unname(y[country]) },
    nearest_in_range = { r <- as.integer(strsplit(arg, "-")[[1]]); min(max(sy, r[1]), r[2]) },
    epoch5_nearest = { r <- as.integer(strsplit(arg, "-")[[1]]); e <- seq(r[1], r[2], by = 5); e[which.min(abs(e - sy))] },
    stop("unknown year rule: ", rule))
}
long <- bind_rows(lapply(seq_along(cols), function(i) {
  rule <- R$year_rule[rule_of[i]]
  data.frame(column = cols[i], rule_id = R$rule_id[rule_of[i]], country = COUNTRIES,
             survey_year = unname(unlist(SURVEY_YEAR[COUNTRIES])),
             year_used = vapply(COUNTRIES, function(cc) as.integer(year_for(cols[i], rule, cc)), NA_integer_),
             has_data = unname(has_data[i, COUNTRIES]), stringsAsFactors = FALSE)
}))
long$year_offset <- long$year_used - long$survey_year
write.csv(long, sprintf("results/tables/predictor_alignment_%s.csv", STAMP), row.names = FALSE)

M$alignment_rule <- R$rule_id[rule_of]
M$year_used <- vapply(cols, function(v) {
  d <- long[long$column == v, ]
  if (all(is.na(d$year_used))) return("static")
  if (length(unique(d$year_used)) == 1L) return(as.character(d$year_used[1]))
  paste(sprintf("%s=%d", d$country, d$year_used), collapse = ";") }, "")
M$year_offset_max_abs <- vapply(cols, function(v) { d <- long[long$column == v & long$has_data, ]; if (!nrow(d) || all(is.na(d$year_offset))) NA_real_ else max(abs(d$year_offset), na.rm = TRUE) }, 0)

write.csv(M, file.path(HDIR, "predictors_admin2_shared_metadata.csv"), row.names = FALSE)

# ── report ───────────────────────────────────────────────────────────────────
cat(sprintf("\n[stamp] %d predictors | tiers: %s | national constants: %d | modelled surfaces: %d\n", nrow(M),
            paste(sprintf("%s %d", names(table(M$tier)), table(M$tier)), collapse = ", "), sum(!M$subnational), sum(M$modelled_surface)))
dd <- sort(table(M$domain), decreasing = TRUE); cat("domains:", paste(sprintf("%s %d", names(dd), dd), collapse = " | "), "\n")
al <- long |> filter(has_data, !is.na(year_offset)) |> group_by(rule_id) |>
  summarise(columns = n_distinct(column), mean_abs_offset = round(mean(abs(year_offset)), 2), max_abs_offset = max(abs(year_offset)),
            worst_country = country[which.max(abs(year_offset))], .groups = "drop") |> arrange(desc(max_abs_offset), desc(columns))
cat("\nalignment by rule (columns with data; static layers excluded):\n"); print(as.data.frame(al), row.names = FALSE)
cat(sprintf("\ncolumns with |offset| >= 3 years in some country with data: %d of %d dated columns\n",
            sum(M$year_offset_max_abs >= 3, na.rm = TRUE), sum(!is.na(M$year_offset_max_abs))))
cat("DONE\n")
