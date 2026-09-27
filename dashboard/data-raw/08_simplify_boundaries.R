# =============================================================================
# dashboard/data-raw/08_simplify_boundaries.R
#
# The boundary bundles ship full GADM resolution (admin2 ~11 MB of the app's
# ~17 MB), which is most of the first-load time on shinyapps. This simplifies
# them once with mapshaper (topology-preserving, so neighbouring districts keep
# shared borders and no slivers open up), after backing up the originals as
# *_full.rds. Validation before writing: every country keeps its feature count
# and every Admin1|Admin2 key, all geometries valid, no empty geometries.
#
#   Rscript dashboard/data-raw/08_simplify_boundaries.R
# =============================================================================
suppressPackageStartupMessages({library(sf); library(rmapshaper); library(here)})
setwd(here::here())
D <- "dashboard/data"
KEEP <- list(admin2_boundaries = 0.10, admin1_boundaries = 0.20)

for (nm in names(KEEP)) {
  f <- file.path(D, paste0(nm, ".rds"))
  bak <- file.path(D, paste0(nm, "_full.rds"))
  b <- readRDS(f)
  if (!file.exists(bak)) { file.copy(f, bak); cat("backed up ->", bak, "\n") }
  out <- list()
  for (ck in names(b)) {
    x <- sf::st_make_valid(b[[ck]])
    s <- rmapshaper::ms_simplify(x, keep = KEEP[[nm]], keep_shapes = TRUE, sys = FALSE)
    s <- sf::st_make_valid(s)
    keycols <- intersect(c("Admin1", "Admin2"), names(x))
    k0 <- do.call(paste, c(sf::st_drop_geometry(x)[keycols], sep = "|"))
    k1 <- do.call(paste, c(sf::st_drop_geometry(s)[keycols], sep = "|"))
    stopifnot(nrow(s) == nrow(x), setequal(k0, k1), !any(sf::st_is_empty(s)), all(sf::st_is_valid(s)))
    out[[ck]] <- s
    cat(sprintf("  %-14s %-12s %4d features, vertices %7d -> %6d\n", nm, ck, nrow(s),
                sum(vapply(sf::st_geometry(x), function(g) nrow(sf::st_coordinates(g)), 1L)),
                sum(vapply(sf::st_geometry(s), function(g) nrow(sf::st_coordinates(g)), 1L))))
  }
  saveRDS(out, f)
  cat(sprintf("%s: %.1f MB -> %.1f MB\n", nm, file.size(bak) / 1e6, file.size(f) / 1e6))
}
cat("done; rerun smoke_test.R (joins) and 07_build_country_briefs.R (map images) before deploying\n")
