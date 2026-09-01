# WS-G. Côte d'Ivoire: an anchored export, and why not on Ghana

`results/deliverables/civ_anchored_regions.csv`,
`docs/requests/civ_labels_request.md`.

## The proposal, and the number that settled it

The proposal was to apply this project's relative-ranking machinery to Côte
d'Ivoire, anchoring on Ghana's national mean. The anchoring half is right. The
Ghana half is not, and VMNIS holds both countries so it can be checked directly:

| Non-pregnant women | Côte d'Ivoire (2007) | Ghana (2017) | gap |
|:---|---:|---:|---:|
| Folate deficiency | **86.1%** | 53.8% | **32 pp** |
| Vitamin B12 deficiency | **18.1%** | 6.9% | **11 pp** |

Borrowing Ghana's level would inject a 32-point error, which is precisely the
failure mode WS-I exists to route around. **Côte d'Ivoire's own VMNIS entry is
used instead**, and it is available at nine sub-national units, so the anchor
sits at region rather than national resolution. WS-5's budget curve puts Admin-1
error well below Admin-2 at every survey fraction, so that is the better product
regardless.

## What shipped

24 rows: 6 national, 18 eco-region, covering **folate and vitamin B12 only**.

The B12 gradient is real targeting signal -- a 48-point spread, coherent
north to south:

| Vitamin B12, worst first | prevalence | vs national |
|:---|---:|---:|
| North | 48.6% | +30.5 |
| North East | 40.4% | +22.3 |
| North West | 36.4% | +18.3 |
| West | 29.8% | +11.7 |
| Central | 20.0% | +1.9 |
| Central West | 11.5% | −6.6 |
| South East | 7.0% | −11.1 |
| South | 6.2% | −11.9 |
| Abidjan | 0.0% | −18.1 |

**Folate is not usable for targeting.** It runs 74 to 96 percent against a
national 86.1 percent. At that level nearly everyone is deficient and the
ranking carries little decision value, whatever its internal ordering.

## No covariate tilt, for two independent reasons

**It is not feasible.** The pooled model matches covariates by exact name and
treats an absent name as zero rather than as an error. The harmonised set is 294
columns dominated by blocks Côte d'Ivoire has none of -- `dhs` 97, `aef` 64,
`soil` 38, `tclim` 13, `spam` 12, `precip` 10. Its 83 cached GEE rasters would
populate roughly 20 to 30 names and silently zero the rest, producing a
confident number with nothing behind it.

**It would not help.** WS-I measured the flat anchor at 8.94 pp against 17.49 pp
for the un-anchored transported model, with the tilt worth nothing pooled.

## Every caveat travels with the rows

| Column | Value |
|:---|:---|
| `scored` | **FALSE** on all 24 rows — no CIV biomarker labels exist in this project |
| `cutoff_verified` | **FALSE** on all 24 rows |
| `years_stale` | **19** — the survey is 2007 |
| `unit_is_administrative` | **FALSE** — the nine units are survey eco-regions, not GADM |
| `covariate_tilt_applied` | FALSE, with the reason |

`Deficiencycutoff` is NA for Côte d'Ivoire *and* Ghana, so part of that 32-point
folate gap may be definitional rather than biological. It cannot be checked from
the data in hand, and the export says so rather than implying the gap is
entirely real.

## What would unlock more

`docs/requests/civ_labels_request.md` lists it: individual biomarker records
with cluster IDs (for any scored estimate at all), iron and vitamin A national
estimates (absent from VMNIS for CIV), the 2007 deficiency cutoffs, a crosswalk
from eco-regions to GADM, and a survey more recent than 2007.
