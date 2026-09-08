# =============================================================================
# scripts/protocol_v2/53_apply_exclusion_policy.R   [LK-01]
#
# APPLY metadata/covariates/exclusions.csv TO THE LIVE SHARED SET
#
# The shared set is rebuilt by scripts/covariates/build_shared_predictor_set.R
# only rarely (two rebuild traps, see the protocol-v2 memory note), and its
# extra derivers used to bypass the exclusion file (AU-01 finding 5). This
# applies every rule in exclusions.csv to the live set in place, so a policy
# change reaches the vocabulary without a rebuild: matched columns are dropped
# from predictors_admin2_shared.csv and its metadata (backup *.pre_policy the
# first time), and the drops are listed with the rule that caught them.
#
#   Rscript scripts/protocol_v2/53_apply_exclusion_policy.R          # apply
#   V2_POLICY_DRY=1 Rscript scripts/protocol_v2/53_apply_exclusion_policy.R   # list only
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"; DRY <- identical(Sys.getenv("V2_POLICY_DRY", "0"), "1")
EX <- read.csv("metadata/covariates/exclusions.csv", stringsAsFactors = FALSE)
if (!"policy" %in% names(EX)) EX$policy <- "data_defect"
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
M <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
cols <- setdiff(names(S), c("country", "Admin1", "Admin2"))
cat(sprintf("live shared set: %d rows x %d predictors; %d exclusion rules (%d leakage)\n", nrow(S), length(cols), nrow(EX), sum(EX$policy == "leakage")))
hits <- list()
for (i in seq_len(nrow(EX))) {
  h <- cols[grepl(EX$canonical_regex[i], cols, perl = TRUE)]
  if (length(h)) hits[[length(hits) + 1L]] <- data.frame(column = h, policy = EX$policy[i], rule = EX$canonical_regex[i], stringsAsFactors = FALSE)
}
H <- if (length(hits)) bind_rows(hits) |> distinct(column, .keep_all = TRUE) else data.frame(column = character(0), policy = character(0), rule = character(0))
if (!nrow(H)) { cat("no live column matches any rule; nothing to do\n"); quit(save = "no") }
cat(sprintf("\n%d column(s) match:\n", nrow(H)))
for (i in seq_len(nrow(H))) cat(sprintf("  %-40s %-12s %s\n", H$column[i], H$policy[i], substr(H$rule[i], 1, 60)))
if (DRY) { cat("\n(dry run; nothing written)\n"); quit(save = "no") }
for (f in c("predictors_admin2_shared.csv", "predictors_admin2_shared_metadata.csv"))
  if (!file.exists(file.path(HDIR, paste0(f, ".pre_policy")))) file.copy(file.path(HDIR, f), file.path(HDIR, paste0(f, ".pre_policy")))
S2 <- S[, !names(S) %in% H$column]; M2 <- M[!M$column %in% H$column, ]
write.csv(S2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(M2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"), row.names = FALSE)
write.csv(H, "results/tables/protocol_v2/exclusion_policy_applied.csv", row.names = FALSE)
cat(sprintf("\nshared set: %d -> %d predictors; metadata %d -> %d rows\n-> results/tables/protocol_v2/exclusion_policy_applied.csv\nDONE\n", length(cols), ncol(S2) - 3, nrow(M), nrow(M2)))
