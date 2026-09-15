# =============================================================================
# Checked Admin-2 joins (JK-01, 2026-09-15): the spine with GADM codes, the
# fan-refusing pair-key join, the spine-validated population loader and the
# reviewed name matcher.
# =============================================================================

source(here::here("R", "admin2_key_hygiene.R"))
source(here::here("R", "admin2_keys.R"))

test_that("the spine has one row per polygon, unique names and unique GADM codes, and matches the shared set", {
  p <- here::here("metadata", "admin2_spine.csv"); skip_if_not(file.exists(p), "spine absent")
  sp <- admin2_spine()
  expect_equal(nrow(sp), 554L)
  expect_false(anyDuplicated(sp$gid_2) > 0)
  expect_true(all(grepl("^(GMB|GHA|MWI|SLE)", sp$gid_2)))   # GADM 4.1 codes Ghana as GHA1.1_2, the others as GMB.1.1_1
  s <- here::here("data", "covariates", "harmonized", "predictors_admin2_shared.csv")
  skip_if_not(file.exists(s), "shared set absent")
  S <- read.csv(s, check.names = FALSE, stringsAsFactors = FALSE)[, c("country", "Admin1", "Admin2")]
  expect_setequal(paste(S$country, S$Admin1, S$Admin2), paste(sp$country, sp$Admin1, sp$Admin2))
})

test_that("join_admin2_v2 joins on the pair key, refuses to fan and reports unmatched rows", {
  x <- data.frame(Admin1 = c("A", "A", "B"), Admin2 = c("Lundu", "Mbwana", "Lundu"), y = 1:3, stringsAsFactors = FALSE)
  y <- data.frame(Admin1 = c("A", "B"), Admin2 = c("Lundu", "Lundu"), pop = c(10, 20), stringsAsFactors = FALSE)
  out <- join_admin2_v2(x, y, what = "test", quiet = TRUE)
  expect_equal(nrow(out), 3L); expect_equal(out$pop, c(10, NA, 20))
  # the same-name pair in two regions is exactly the case a name-only join fans
  expect_equal(nrow(join_admin2_v2(x, y, how = "inner", quiet = TRUE)), 2L)
  y_dup <- rbind(y, data.frame(Admin1 = "A", Admin2 = "Lundu", pop = 99))
  expect_error(join_admin2_v2(x, y_dup, quiet = TRUE), "would fan rows")
  expect_error(join_admin2_v2(x, y[, c("Admin2", "pop")], quiet = TRUE), "need Admin1 and Admin2")
})

test_that("the population table loads for every country on the pair key and validates against the spine", {
  p <- here::here("dashboard", "data", "admin2_population.rds"); skip_if_not(file.exists(p), "population table absent")
  skip_if_not(file.exists(here::here("metadata", "admin2_spine.csv")), "spine absent")
  POP <- readRDS(p)
  for (cn in c("Gambia", "Ghana", "Malawi", "SierraLeone")) {
    pp <- admin2_population_v2(POP, cn, "pop_child")
    expect_true(nrow(pp) > 0, info = cn)
    expect_equal(names(pp), c("Admin1", "Admin2", "pop"))
    expect_false(anyDuplicated(paste(pp$Admin1, pp$Admin2)) > 0, info = cn)
  }
  # the file's own "Sierra Leone" spelling must not drop the country
  expect_gt(nrow(admin2_population_v2(POP, "SierraLeone", "pop_child")), 10)
  # JK-02: every spine unit has a population row (the name-only dedupe that
  # dropped one TA of each Malawi same-name pair is fixed)
  sp <- admin2_spine()
  for (cn in c("Gambia", "Ghana", "Malawi", "SierraLeone")) {
    pp <- admin2_population_v2(POP, cn, "population"); s1 <- sp[sp$country == cn, ]
    expect_true(all(paste(s1$Admin1, s1$Admin2) %in% paste(pp$Admin1, pp$Admin2)), info = paste(cn, "units without population"))
  }
})

test_that("admin2_match_v2 matches exact, then alias, then fuzzy, and writes a review file", {
  tmp <- withr::local_tempdir()
  al <- file.path(tmp, "aliases.csv"); writeLines(c("source,target", "Sabach Sanjal,Upper Baddibu"), al)
  rev <- file.path(tmp, "review.csv")
  out <- admin2_match_v2(c("upper baddibu", "Sabach Sanjal", "Janjabureh", "Nowhere"), c("Upper Baddibu", "Janjanbureh", "Kuntaur"),
                         aliases_csv = al, review_csv = rev, label = "test")
  expect_equal(out, c("Upper Baddibu", "Upper Baddibu", "Janjanbureh", NA))
  r <- read.csv(rev, stringsAsFactors = FALSE)
  expect_setequal(r$method[!is.na(r$target)], c("exact", "alias", "fuzzy"))
  expect_true(is.na(r$target[r$source == "Nowhere"]))
  # an alias pointing outside the target vocabulary is an error, not a silent NA
  writeLines(c("source,target", "Sabach Sanjal,Upper Badibu"), al)
  expect_error(admin2_match_v2("Sabach Sanjal", c("Upper Baddibu"), aliases_csv = al), "not in the target vocabulary")
})

test_that("admin2_join_by reports the first degradation to a name-only key", {
  .a2_degrade_env$reported <- NULL
  x <- data.frame(Admin2 = "a"); y <- data.frame(Admin1 = "r", Admin2 = "a")
  expect_message(admin2_join_by(x, y), "name-only")
  expect_silent(admin2_join_by(x, y))   # reported once per session
  expect_equal(admin2_join_by(y, y), c("Admin1", "Admin2"))
})
