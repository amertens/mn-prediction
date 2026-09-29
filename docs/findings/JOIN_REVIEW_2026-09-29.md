# Critical review: Admin-2 name cleaning, merging, and whether finer units would help

*2026-09-29. Scope: the district-name harmonisation and join machinery across
`R/`, `scripts/`, `dashboard/`, and the question of whether smaller admin units
— Sierra Leone chiefdoms in particular — would increase usable sample size.*

---

## Part 1 — is the cleaning and merging implemented correctly?

### The short answer

**Where it matters, yes.** The end-to-end merge is exact and complete:

| check | result |
|---|---|
| spine units vs `predictors_admin2_shared.csv` | 554 vs 554, **0** mismatched either direction |
| outcome units (`targets_v2.csv`) with no predictor row | **0** of 206 |
| fuzzy matches ever made across all crosswalk builders | **3**, all three correct |

There are no silent drops and no orphaned units. The machinery behind this is
genuinely well built: `R/admin2_keys.R` gives a pair-keyed join that refuses to
fan rows, validates against the spine, and logs unmatched keys;
`R/lint_admin2_joins.R` is a ratchet lint with a recorded baseline and a
testthat test, written after the name-only-join defect "recurred ten times".

Four defects remain. They are ranked by what they actually cost.

---

### 1. `scripts/protocol_v2/32_admin1_aggregation_weights.R` is broken, and re-running it destroys its own result

Line 47:

```r
sc <- S[S$country == cn, c("Admin1", "Admin2", PREDS)] |> left_join(pp, by = "Admin2")
```

**`pp` is never defined** — not in the script, not in `R/protocol_v2.R`, not in
anything it sources. Verified by sourcing `R/protocol_v2.R` in a clean session:
`exists("pp")` is `FALSE`.

The failure is silent by construction. `build_a1()` is called as

```r
z <- tryCatch(build_a1(cn, on, target, scheme), error = function(e) NULL)
...
if (length(cl) < 3) next
```

so the "object not found" error is swallowed for every country, every
combination is skipped, and the script ends with

```r
R <- bind_rows(rows); write.csv(R, file.path(OUTDIR, "admin1_aggregation_weights.csv"), ...)
```

`bind_rows(list())` is a 0-row frame, so **the write overwrites the existing
132-row AG-01 result with an empty file** and prints no error.

The current `admin1_aggregation_weights.csv` (132 rows, written 9 Sep) was
therefore produced in a session where a leftover `pp` existed in the global
environment. Line 47 is unchanged since the file was created on 4 Sep
(`ca0c5d4`), so that result is **not reproducible from the script as committed**.

This is the one live name-only join in the whole protocol-v2 path, which is
otherwise clean: 70 `by = c("Admin1","Admin2")`, 13 `by = c("country","Admin1",
"Admin2")`, 12 `join_admin2_v2()`.

It also matters *which* table it joins. `targets_v2.csv` (87 Malawi outcome
units) has **no** duplicated Admin2 names, so name-only joins against outcomes
are harmless. `predictors_admin2_shared.csv` (243 Malawi rows) has **four**
— `TA Lundu`, `TA Ngabu`, `TA Pemba`, `TA Malemia`, each in two Admin1 regions.
Line 47 joins `S`, the predictor table. A name-only self-join there takes
243 rows to 251.

**Fix:** define `pp` (evidently the population table — the line above builds it
correctly via `admin2_population_v2()`), join on `c("Admin1","Admin2")`, and
replace the blanket `tryCatch` with one that reports rather than returns `NULL`.

---

### 2. The linter has a blind spot that hid exactly that site

```r
is_pair <- grepl("[\"']Admin1[\"']", code)
hit <- which(grepl(p$rx, code, perl = TRUE) & !is_pair)
```

Any line containing the literal `"Admin1"` **anywhere** is treated as a
pair-key join and skipped. The intent (documented in the file) was to avoid
matching inside `by = c("Admin1","Admin2")`, but the test is line-level, not
match-level. Selecting the pair columns and *then* joining on the name escapes:

```r
S[S$country == cn, c("Admin1", "Admin2", PREDS)] |> left_join(pp, by = "Admin2")
#                    ^^^^^^^^ satisfies is_pair            ^^^^^^^^^^^^^^ name-only
```

Re-scanning with the Admin1 test applied **to the `by=` clause only** finds
exactly one such site across `R/`, `scripts/` and `dashboard/`: script 32
line 47. The blind spot is narrow, but it hid the one defect that mattered.

**Fix:** decide `is_pair` from the matched `by=` argument, not the whole line.

---

### 3. All 59 grandfathered join sites are marked `"not assessed"`

The ratchet stopped *new* name-only joins but nobody triaged the existing ones:
every row of `tests/testthat/admin2_join_baseline.csv` carries
`audit = "not assessed"`. The linter's own header states the decision rule —
a name-only join is safe only when neither side can carry a duplicate name —
so the `audit` column is the whole point and it is empty.

The rule is now cheap to apply, because the duplicate names are known exactly:
**a site is safe unless one side is the 243-row Malawi predictor set** (or any
table spanning multiple Malawi Admin1 regions). Sites worth checking first are
the four in `scripts/malawi_admin2_oos.R` and the `R/corrected/*` merges, since
`R/corrected/` is the audited layer the main pipeline is reconciled against.

---

### 4. The fuzzy threshold is looser than the class separation

`admin2_match_v2()` accepts any Jaro-Winkler match with `max_jw = 0.15`. The
distance between *genuinely different* districts is far smaller than that:

| country | closest distinct pair | JW | distinct pairs within 0.15 |
|---|---|---|---|
| Ghana | `ahafoanosoutheast` ~ `ahafoanosouthwest` | 0.024 | **147** |
| Malawi | `tampando` ~ `tamponda` | 0.025 | **322** |
| Gambia | `fulladueast` ~ `fulladuwest` | 0.036 | 15 |
| Sierra Leone | `westernrural` ~ `westernurban` | 0.087 | 1 |

So the matcher decides between "Ahafo Ano South East" and "South West" on a
margin of 0.024 while tolerating error up to 0.15 — a threshold six times the
separation it must resolve. The East/West and North/South pairs that dominate
Ghanaian and Gambian district names are precisely the ones it cannot safely
discriminate.

**The realised risk today is zero**: only three fuzzy matches were ever made
(`Janjanbureh→Janjabureh`, `WESTERN AREA RURAL/URBAN→Western Rural/Urban`) and
all are correct. The exposure is latent, and it is bounded by the review CSVs,
which are written for every call. But the matcher is only ever used on small
vocabularies (8–28 units); nothing has yet pointed it at Ghana's 260 or
Malawi's 243, where it would be unsafe.

**Fix:** reject a match unless the best candidate beats the second-best by a
margin (e.g. `d2 - d1 > d1`), and refuse to match Admin2 names without an
Admin1 restriction when the target vocabulary has near-duplicates.

---

### Smaller notes

- **`CLAUDE.md` is wrong about testing.** It says "There is no test suite (no
  `testthat`)". There are 15 test files under `tests/testthat/`, including
  `test-admin2-join-lint.R`, `test-admin2-keys.R` and `test-leakage-guard.R`.
  Anyone trusting the doc will not run them.
- **Three different name normalisers are in use.** `admin2_kk()` and several
  scripts use `[^a-z0-9]` (keeps digits); at least nine scripts use `[^a-z]`
  (strips digits); `ws1c_simulation.R` strips spaces only. No current admin
  name contains a digit, so nothing is broken today — but two of the three
  would silently disagree on a name like "Bo 2".

---

## Part 2 — would smaller admin units increase sample size?

### The premise is inverted for Sierra Leone

Sierra Leone has the **fewest areas but by far the best-measured ones**:

| country | areas | median effective n per area | min |
|---|---|---|---|
| Gambia | 30 | 13.2 | 4.9 |
| Ghana | 75 | 4.8 | 1.0 |
| Malawi | 87 | 4.6 | 1.0 |
| **Sierra Leone** | **14** | **22.4** | 6.3 |

Ghana and Malawi already *are* the fine-grained case, and their per-area
effective n has collapsed below 5. Sierra Leone is the one country whose areas
carry enough data to estimate a prevalence at all. It is therefore the country
where subdividing costs the most, not the least.

### What chiefdoms would actually give (measured, not argued)

Sierra Leone's survey is **1,477 individuals in 60 clusters across 14
districts** — 2 to 8 clusters per district, median 4. Assigning those 60
clusters to the OCHA `sle_admin3` chiefdom polygons by point-in-polygon
(0 fell outside):

| | |
|---|---|
| chiefdoms in the country | 167 |
| chiefdoms containing ≥1 survey cluster | **42** (25%) |
| of those, chiefdoms with **exactly one** cluster | **29** (69%) |
| chiefdoms with ≥3 clusters | 3 |

A chiefdom holding one cluster **is** that cluster. So a chiefdom panel is the
cluster panel for 29 of its 42 populated rows, and 125 of 167 chiefdoms would
carry no data at all. Median effective n falls from 22.4 (district) to 9.0
(cluster).

Subdividing does not create information. It redistributes the same 60 clusters
into more, noisier buckets.

### And the extreme case has already been tested

`scripts/cluster_level/` fitted the full 323-cluster panel across all four
countries under the same three estimands. Cluster-fitted models **do not beat
district-fitted on any estimand**, and transport is clearly worse — climate +
soil index at Admin-2, **0.244 vs 0.368** on the biomarker level. Script `04`
added Kish-n weighting and empirical-Bayes shrinkage of cluster outcomes toward
the district, which "recovers a third to a half of the gap, never all of it".

Chiefdom-level sits between district and cluster and, for Sierra Leone,
overwhelmingly at the cluster end. There is no reason to expect it to land
above a result that the cluster panel already failed to reach.

### The honest counter-argument, and why it does not rescue this

More rows do help *fitting*: the index is currently estimated on 14 Sierra
Leone rows, and 42 would give more degrees of freedom even with a noisier
outcome. That is a real errors-in-variables trade rather than a free lunch, and
it is exactly what script `04` tested with shrinkage — the honest instrument
for it. It recovered part of the gap and never closed it.

### Where the sample-size constraint actually binds

Not Sierra Leone. Ghana (n_eff 4.8) and Malawi (n_eff 4.6) are where areas are
too thinly measured, and for them the promising direction is the **opposite**:
aggregate up, trading spatial resolution for measurement precision. That is
precisely what AG-01 was built to test — and AG-01 is script 32, the one that
is currently broken. Fixing defect 1 is therefore the prerequisite for
answering the resolution question properly.

---

## Recommended order of work

1. Fix `pp` in script 32 and re-run AG-01 (defect 1). Until then its result is
   unreproducible and a re-run silently empties it.
2. Narrow the linter's `is_pair` test to the `by=` clause (defect 2) and
   re-baseline.
3. Triage the 59 grandfathered sites against the now-explicit rule, starting
   with anything touching the 243-row Malawi predictor table (defect 3).
4. Add a margin rule to `admin2_match_v2()` before it is ever pointed at a
   large vocabulary (defect 4).
5. Correct the testing claim in `CLAUDE.md`.

Sierra Leone chiefdoms: **do not pursue.** The measurement above is the reason,
and it is cheap to re-check if the survey is ever extended.
