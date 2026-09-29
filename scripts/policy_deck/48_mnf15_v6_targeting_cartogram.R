# =============================================================================
# scripts/policy_deck/48_mnf15_v6_targeting_cartogram.R
#
# One figure for the v6 MNF15 talk: why "rank by rate" and "go to the most
# populous districts" pick different places (TC-01, script 74).
#
#   Left:  all 260 Ghana districts, coloured by thirds of the model's predicted
#          rate of iron deficiency in children 6-59 months (worst third darkest).
#   Right: the same districts as non-overlapping circles at their centroids,
#          circle area proportional to children 6-59 months (Dorling-style
#          cartogram, laid out here in base R: packcircles / cartogram are not
#          installed), same colours, Ghana's outline faint behind.
#
# The model is the deployed index: fitted on Ghana's 75 surveyed districts
# (targets_v2.csv, child_iron, y_prev on the logit scale via .v2_logit) and
# applied to every district, with predictors built for all 260 districts as in
# the "Model, all 260 districts" panel of 26_mnf15_v2_figures.R (rank-normalised
# over all 260 districts, domain PCs with rotations learned on the surveyed
# rows). The rate is the calibrated index (domain_index_cal, prevalence target,
# as script 74's in-country rate), expit of its logit. Its order equals the
# ranking index (domain_index) exactly; the script checks this and reports the
# Spearman against script 26's all-district panel (level target), recomputed
# here because that panel's values are not exported.
#
# Population: dashboard/data/admin2_population.rds, pop_child, on the Admin1 +
# Admin2 pair key (admin2_population_v2 / join_admin2_v2).
#
#   Rscript scripts/policy_deck/48_mnf15_v6_targeting_cartogram.R
# -> results/figures/mnf15_v6/a6_targeting_cartogram.png
#    results/tables/policy_deck/ghana_child_iron_priority_all_districts.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(sf); library(patchwork)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")          # the benchmark's headline set
Sys.unsetenv(c("V2_INDEX_SHRINK", "V2_PREP_SCALE", "V2_DOMAIN_REP"))   # protocol defaults
source("R/protocol_v2.R")

OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
P2 <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"
INK <- "#1A1A1A"; GREY <- "grey40"
FILL3 <- c(worst = "#0F7B8A", middle = "#7FB8C0", best = "#E4E7E8")
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)
key <- function(d) paste(d$Admin1, d$Admin2, sep = "||")
N_ALL <- 260L; TOPFRAC <- 0.20

# ---------------------------------------------------------------- data
TG  <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
B_gh <- readRDS("dashboard/data/admin2_boundaries.rds")[["ghana"]]      # read only
POP  <- readRDS("dashboard/data/admin2_population.rds")
chk(nrow(B_gh) == N_ALL && !anyDuplicated(key(B_gh)), "260 unique Ghana districts in the boundary file")

# ---------------------------------------------------------------- model, all 260 districts
# (26_mnf15_v2_figures.R, "every district: fit on the 75 surveyed rows, apply to all 260")
all5 <- S[S$country == "Ghana", c("Admin1", "Admin2", PREDS)]
chk(nrow(all5) == N_ALL && setequal(key(all5), key(B_gh)), "predictor table covers the 260 boundary districts")
Xa5 <- prep_predictors_v2(as.matrix(all5[, PREDS]))

# prevalence target (the figure's rate)
tp <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_prev), is.finite(n_eff))
tr_p <- match(key(tp), key(all5)); chk(!anyNA(tr_p) && length(tr_p) == 75, "75 surveyed rows in the all-district table")
Yp <- rep(NA_real_, N_ALL); Yp[tr_p] <- tp$y_prev
Da_p <- domain_representation_v2(Xa5, domain_of, sign_rows = tr_p)
z_cal <- ARMS_V2[["domain_index_cal"]](tr_p, seq_len(N_ALL), .v2_logit(Yp), Xa5, Da_p, list(y_nat = Yp, target = "prev"))
z_idx <- ARMS_V2[["domain_index"]](tr_p, seq_len(N_ALL), .v2_logit(Yp), Xa5, Da_p, list(y_nat = Yp))
chk(isTRUE(all.equal(rank(z_cal), rank(z_idx))), "calibrated rate and ranking index give the same order")
all5$pred_rate <- .v2_expit(z_cal)
cat(sprintf("fit: %d surveyed of %d; %d predictors -> %d domain axes; mean predicted rate %.1f%% (surveyed mean %.1f%%)\n",
            length(tr_p), N_ALL, ncol(Xa5), ncol(Da_p), 100 * mean(all5$pred_rate), 100 * mean(tp$y_prev)))

# sanity check: script 26's all-district panel (level target), recomputed exactly as there
tl <- TG |> filter(country == "Ghana", outcome == "child_iron", is.finite(y_level), is.finite(n_eff_cont))
tr_l <- match(key(tl), key(all5)); chk(!anyNA(tr_l), "level-target surveyed rows in the all-district table")
Yl <- rep(NA_real_, N_ALL); Yl[tr_l] <- tl$y_level
Da_l <- domain_representation_v2(Xa5, domain_of, sign_rows = tr_l)
pred26 <- ARMS_V2[["domain_index"]](tr_l, seq_len(N_ALL), Yl, Xa5, Da_l, list(y_nat = Yl))
rho26 <- cor(all5$pred_rate, pred26, method = "spearman")
third_of <- function(v) dplyr::ntile(-v, 3)            # 1 = highest predicted rate (87, 87, 86 districts)
agree26 <- mean(third_of(all5$pred_rate) == third_of(pred26))
cat(sprintf("sanity: Spearman with script 26's all-district panel (level target) %.3f; same third for %.0f%% of districts\n",
            rho26, 100 * agree26))
chk(rho26 > 0.8, "ordering close to script 26's all-district panel")
# for reference only: the dashboard's deployment ranking (built 27 Sep, its own tier set and fit)
IX <- readRDS("dashboard/data/admin2_index.rds")$districts
ix <- IX[IX$country == "Ghana" & IX$outcome == "child_iron", ]
if (nrow(ix) == N_ALL) cat(sprintf("reference: Spearman with the dashboard's deployment score %.3f\n",
                                   cor(all5$pred_rate, ix$score_logit[match(key(all5), key(ix))], method = "spearman")))

# ---------------------------------------------------------------- population and thirds
pp <- admin2_population_v2(POP, "Ghana", "pop_child")
D <- join_admin2_v2(all5[, c("Admin1", "Admin2", "pred_rate")], pp, how = "left", what = "pop Ghana child")
chk(all(is.finite(D$pop) & D$pop > 0), "every district has a child population")
D$children <- D$pop
D$rank_worst <- rank(-D$pred_rate, ties.method = "first")
D$third <- factor(c("worst", "middle", "best")[third_of(D$pred_rate)], levels = c("worst", "middle", "best"))
D$surveyed <- key(D) %in% key(tp)
tot <- sum(D$children)
n_worst <- sum(D$third == "worst"); sh_d_worst <- n_worst / N_ALL
sh_c3 <- tapply(D$children, D$third, sum) / tot
k20 <- round(TOPFRAC * N_ALL)                                          # 52 districts, script 12's rule
top_pop <- order(D$children, decreasing = TRUE)[seq_len(k20)]
sh_c_pop20 <- sum(D$children[top_pop]) / tot
sh_c_rate20 <- sum(D$children[D$rank_worst <= k20]) / tot
ov_pop_rate <- sum(D$rank_worst[top_pop] <= k20)
exp_cases <- D$pred_rate * D$children                                 # the model's expected deficient children
sh_x_worst <- sum(exp_cases[D$third == "worst"]) / sum(exp_cases)
sh_x_pop20 <- sum(exp_cases[top_pop]) / sum(exp_cases)
sh_x_rate20 <- sum(exp_cases[D$rank_worst <= k20]) / sum(exp_cases)
cat(sprintf("\nworst third: %d of %d districts (%.1f%%), %.1f%% of children 6-59 months (middle %.1f%%, best %.1f%%)\n",
            n_worst, N_ALL, 100 * sh_d_worst, 100 * sh_c3["worst"], 100 * sh_c3["middle"], 100 * sh_c3["best"]))
cat(sprintf("most populous fifth (%d districts): %.1f%% of children; highest-rate fifth: %.1f%% of children; %d districts in both\n",
            k20, 100 * sh_c_pop20, 100 * sh_c_rate20, ov_pop_rate))
cat(sprintf("model's expected deficient children: worst third %.1f%%, highest-rate fifth %.1f%%, most populous fifth %.1f%%\n",
            100 * sh_x_worst, 100 * sh_x_rate20, 100 * sh_x_pop20))
cat(sprintf("predicted rate: worst third %.1f-%.1f%%, best third %.1f-%.1f%%\n",
            100 * min(D$pred_rate[D$third == "worst"]), 100 * max(D$pred_rate[D$third == "worst"]),
            100 * min(D$pred_rate[D$third == "best"]), 100 * max(D$pred_rate[D$third == "best"])))
TC <- read.csv(file.path(P2, "tc01_targeting_by_cases_cells.csv"))
tcg <- TC[TC$country == "Ghana" & TC$outcome == "child_iron" & TC$estimand == "infill", c("arm", "framing", "capture")]
cat("TC-01, Ghana children's iron, in-country (surveyed districts, share of deficient children reached):\n")
print(tcg |> mutate(capture = round(100 * capture, 1)), row.names = FALSE)

write.csv(D |> arrange(rank_worst) |>
            transmute(Admin1, district = Admin2, pred_rate = round(pred_rate, 4), rank_worst, third,
                      children = round(children), surveyed),
          file.path(PD, "ghana_child_iron_priority_all_districts.csv"), row.names = FALSE)

# ---------------------------------------------------------------- Dorling-style layout
Bp <- sf::st_transform(B_gh, 32630) |> left_join(D, by = c("Admin1", "Admin2"))
chk(nrow(Bp) == N_ALL && all(is.finite(Bp$children)), "boundary join")
xy0 <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(Bp)))
outline <- sf::st_union(sf::st_make_valid(sf::st_geometry(Bp)))
# circle area proportional to children; all circles together cover CIRC_SHARE of Ghana's land area
CIRC_SHARE <- 0.45
k_r <- sqrt(CIRC_SHARE * as.numeric(sf::st_area(outline)) / (pi * tot))
r <- k_r * sqrt(Bp$children)

#' Pairwise repulsion removes overlaps (the smaller circle of a pair moves
#' more), a weak pull back to each centroid keeps circles near home. The pull
#' decays geometrically, so the last iterations are pure repulsion and the
#' layout ends free of overlaps. Stops once the largest overlap is below `tol`
#' of the smaller radius of the pair.
dorling <- function(x0, y0, r, pull = 0.05, decay = 0.995, step = 0.5, tol = 1e-3, maxit = 20000) {
  x <- x0; y <- y0
  R <- outer(r, r, "+"); rmin <- outer(r, r, pmin)
  share <- outer(r^2, r^2, function(a, b) b / (a + b))    # share of the overlap that circle i takes
  for (it in seq_len(maxit)) {
    dx <- outer(x, x, "-"); dy <- outer(y, y, "-")
    d <- sqrt(dx^2 + dy^2); diag(d) <- Inf
    d[d < 1e-6] <- 1e-6
    ov <- R - d; ov[ov < 0] <- 0
    worst <- max(ov / rmin)
    a <- pull * decay^it
    if (worst < tol && a < 1e-4) break
    f <- ov * share / d
    x <- x + step * rowSums(f * dx) + a * (x0 - x)
    y <- y + step * rowSums(f * dy) + a * (y0 - y)
  }
  list(x = x, y = y, iter = it, max_overlap = worst,
       max_overlap_m = max(ov), moved = sqrt((x - x0)^2 + (y - y0)^2))
}
lay <- dorling(xy0[, 1], xy0[, 2], r)
chk(lay$max_overlap < 1e-3, "circles do not overlap")
cat(sprintf("\nlayout: %d iterations; largest remaining overlap %.2f m (%.3f%% of the smaller radius); radius %.1f-%.1f km; moved median %.1f km, 90th pct %.1f km, max %.1f km\n",
            lay$iter, lay$max_overlap_m, 100 * lay$max_overlap, min(r) / 1e3, max(r) / 1e3,
            median(lay$moved) / 1e3, quantile(lay$moved, 0.9) / 1e3, max(lay$moved) / 1e3))
circ <- sf::st_sf(sf::st_drop_geometry(Bp)[, c("Admin1", "Admin2", "third", "children")],
                  geometry = sf::st_buffer(sf::st_sfc(lapply(seq_len(N_ALL), function(i) sf::st_point(c(lay$x[i], lay$y[i]))),
                                                      crs = 32630), r, nQuadSegs = 20))
circ <- circ[order(-circ$children), ]                                  # small circles drawn on top

# land area by third, for "looks like much of the country"
area_km2 <- as.numeric(sf::st_area(Bp)) / 1e6
sh_a3 <- tapply(area_km2, Bp$third, sum) / sum(area_km2)
cat(sprintf("land area: highest-rate third %.1f%%, middle %.1f%%, lowest %.1f%%\n",
            100 * sh_a3["worst"], 100 * sh_a3["middle"], 100 * sh_a3["best"]))

# ---------------------------------------------------------------- city labels
# Labels sit at the same fixed points (UTM 30N, km) in both panels, outside the
# country or at sea; the leader runs to the district centroid (left) or to the
# edge of the district's circle (right).
CITY <- data.frame(Admin1 = c("Greater Accra", "Ashanti", "Northern", "Western"),
                   Admin2 = c("Accra", "Kumasi", "Tamale", "Sekondi Takoradi"),
                   lab = c("Accra", "Kumasi", "Tamale", "Sekondi-Takoradi"),
                   Lx = 1e3 * c(905, 435, 990, 575), Ly = 1e3 * c(505, 760, 1040, 470),
                   hj = c(0, 1, 0, 1))
ci <- match(key(CITY), key(Bp)); chk(!anyNA(ci), "city districts found in the boundary file")
CITY$x0 <- xy0[ci, 1]; CITY$y0 <- xy0[ci, 2]
CITY$cx <- lay$x[ci]; CITY$cy <- lay$y[ci]; CITY$cr <- r[ci]
u <- sqrt((CITY$Lx - CITY$cx)^2 + (CITY$Ly - CITY$cy)^2)
CITY$ex <- CITY$cx + CITY$cr * (CITY$Lx - CITY$cx) / u; CITY$ey <- CITY$cy + CITY$cr * (CITY$Ly - CITY$cy) / u
print(CITY[, c("lab", "cr")] |> mutate(cr_km = round(cr / 1e3, 1), children = round(D$children[match(key(CITY), key(D))])) |> select(-cr))

# ---------------------------------------------------------------- figure
# One extent for both panels (same scale), widened to about the shape of each
# half of the figure so patchwork does not squeeze the two panels together.
bb <- sf::st_bbox(c(sf::st_geometry(circ), outline))
pad <- 15e3; ASP <- 1.30
LY <- c(min(bb["ymin"], CITY$Ly) - pad, max(bb["ymax"], CITY$Ly) + pad)
LX <- mean(c(bb["xmin"], bb["xmax"])) + c(-0.5, 0.5) * ASP * diff(LY)
chk(all(CITY$Lx > LX[1] + 200e3 | CITY$hj == 0) && all(CITY$Lx < LX[2] - 150e3 | CITY$hj == 1) &&
      LX[1] < bb["xmin"] && LX[2] > bb["xmax"], "labels and circles inside the panel extent")
LEG <- c(worst = "Highest third", middle = "Middle third", best = "Lowest third")
fill_sc <- scale_fill_manual(values = FILL3, labels = LEG, name = "Predicted rate of iron deficiency, children 6-59 months:", drop = FALSE)
th <- theme_void(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 15, colour = INK, hjust = 0.5, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 12.5, colour = GREY, hjust = 0.5, lineheight = 1.05, margin = margin(b = 2)),
        legend.position = "bottom", legend.title = element_text(size = 13, colour = GREY),
        legend.text = element_text(size = 13, colour = INK), legend.key.size = unit(0.5, "cm"))
lab_layer <- function(x, y) geom_label(data = CITY, aes(x = {{ x }}, y = {{ y }}, label = lab, hjust = hj),
                                       size = 4.4, colour = INK, fontface = "bold", linewidth = 0,
                                       fill = scales::alpha("white", 0.85), label.padding = unit(0.08, "lines"))
pL <- ggplot() +
  geom_sf(data = Bp, aes(fill = third), colour = "white", linewidth = 0.1) +
  geom_segment(data = CITY, aes(x = x0, y = y0, xend = Lx, yend = Ly), colour = INK, linewidth = 0.35) +
  geom_point(data = CITY, aes(x0, y0), shape = 21, size = 2.2, stroke = 0.8, colour = INK, fill = "white") +
  lab_layer(Lx, Ly) +
  fill_sc + coord_sf(crs = 32630, datum = NA, xlim = LX, ylim = LY, expand = FALSE) +
  labs(title = "Districts, by predicted priority",
       subtitle = sprintf("Highest-rate third: %.0f%% of districts\nand %.0f%% of the land",
                          100 * sh_d_worst, 100 * sh_a3["worst"])) + th
pR <- ggplot() +
  geom_sf(data = outline, fill = "grey96", colour = "grey78", linewidth = 0.3) +
  geom_sf(data = circ, aes(fill = third), colour = "grey55", linewidth = 0.12, show.legend = FALSE) +
  geom_segment(data = CITY, aes(x = ex, y = ey, xend = Lx, yend = Ly), colour = INK, linewidth = 0.35) +
  lab_layer(Lx, Ly) +
  fill_sc + coord_sf(crs = 32630, datum = NA, xlim = LX, ylim = LY, expand = FALSE) +
  labs(title = "The same districts, sized by number of young children",
       subtitle = sprintf("Highest-rate third: %.0f%% of young children\nMost populous fifth of districts: %.0f%%",
                          100 * sh_c3["worst"], 100 * sh_c_pop20)) + th
fig <- (pL | pR) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT, "a6_targeting_cartogram.png"), fig, width = 12, height = 5.2, dpi = 220, bg = "white")
cat("wrote", file.path(OUT, "a6_targeting_cartogram.png"), "and", file.path(PD, "ghana_child_iron_priority_all_districts.csv"), "\n")
