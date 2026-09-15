# =============================================================================
# Same-survey overlap between a micronutrient survey and the DHS round that
# supplies its proxy predictors.
#
# Malawi's MNS 2015-16 was a subsample of the MDHS 2015-16 (105 of 850
# clusters; the README says to merge on MCLUSTER / MNUMBER / M01, and
# docs/findings/WSC4_MALAWI_BLOCK.md verified 3,097 of 3,099 rows link to a DHS
# person). Every DHS-derived Admin-2 predictor for Malawi must therefore be
# built from the clusters the MNS did NOT sample, or the outcome individuals sit
# inside their own district's predictors. The other three surveys were
# independent samples, so they have no overlap and nothing is filtered.
#
# metadata/mns_dhs_overlap_clusters.csv is the single source of truth, written
# by scripts/covariates/build_mns_overlap_clusters.R.
# =============================================================================

test_that("the overlap helper and its metadata exist", {
  expect_true(file.exists(here::here("R", "mns_dhs_overlap.R")))
  expect_true(file.exists(here::here("metadata", "mns_dhs_overlap_clusters.csv")))
})

if (file.exists(here::here("R", "mns_dhs_overlap.R")))
  source(here::here("R", "mns_dhs_overlap.R"))

test_that("Malawi's MNS subsample is 105 distinct DHS clusters", {
  skip_if_not(exists("mns_overlap_clusters"))
  cl <- mns_overlap_clusters("Malawi")
  expect_type(cl, "integer")
  expect_length(cl, 105L)
  expect_equal(anyDuplicated(cl), 0L)
  expect_true(all(cl >= 1L & cl <= 850L))
})

test_that("the overlap list reproduces the clusters in the Malawi merged dataset", {
  skip_if_not(exists("mns_overlap_clusters"))
  p <- here::here("data", "IPD", "Malawi", "Malawi_merged_dataset.rds")
  skip_if_not(file.exists(p), "Malawi merged dataset absent")
  mns <- readRDS(p)
  expect_setequal(mns_overlap_clusters("Malawi"),
                  sort(unique(as.integer(mns$gw_cnum))))
})

test_that("surveys that were independent of the DHS have no DHS overlap", {
  skip_if_not(exists("mns_overlap_clusters"))
  for (cn in c("Gambia", "Ghana", "SierraLeone", "Sierra Leone", "Tanzania"))
    expect_length(mns_overlap_clusters(cn), 0L)
})

test_that("The Gambia's overlap is with MICS 2018, kept apart from its DHS cluster numbers", {
  skip_if_not(exists("mns_overlap_clusters"))
  p <- here::here("metadata", "mns_dhs_overlap_clusters.csv")
  m <- read.csv(p, stringsAsFactors = FALSE)
  skip_if_not("programme" %in% names(m) && any(m$programme == "MICS"), "MICS overlap rows not yet built")
  cl <- mns_overlap_clusters("Gambia", programme = "MICS")
  expect_equal(length(cl), 70L)
  expect_length(mns_overlap_clusters("Gambia"), 0L)            # DHS 2019-20 is a separate sample
  expect_length(mns_overlap_clusters("Malawi", programme = "MICS"), 0L)
  df <- data.frame(HH1 = c(cl[1], 9999L), x = 1:2)
  expect_equal(nrow(drop_mns_overlap(df, "Gambia", "HH1", quiet = TRUE, programme = "MICS")), 1L)
  expect_equal(nrow(drop_mns_overlap(df, "Gambia", "HH1", quiet = TRUE)), 2L)
})

test_that("drop_mns_overlap removes exactly the overlapping clusters, and only for Malawi", {
  skip_if_not(exists("drop_mns_overlap"))
  cl <- mns_overlap_clusters("Malawi")
  df <- data.frame(v001 = c(cl[1:3], 9001L, 9002L), x = 1:5)
  out <- drop_mns_overlap(df, "Malawi", "v001")
  expect_equal(out$v001, c(9001L, 9002L))
  expect_equal(out$x, 4:5)
  expect_identical(drop_mns_overlap(df, "Ghana", "v001"), df)
  # a column the frame does not carry is a no-op, not an error
  expect_identical(drop_mns_overlap(df, "Malawi", "hv001"), df)
})

test_that("drop_mns_overlap handles haven-labelled cluster ids", {
  skip_if_not(exists("drop_mns_overlap"))
  skip_if_not_installed("haven")
  cl <- mns_overlap_clusters("Malawi")
  v <- haven::labelled(c(cl[1], 9001L), labels = c(one = 1L))
  df <- data.frame(x = 1:2); df$hv001 <- v
  out <- drop_mns_overlap(df, "Malawi", "hv001")
  expect_equal(nrow(out), 1L)
  expect_equal(out$x, 2L)
})
