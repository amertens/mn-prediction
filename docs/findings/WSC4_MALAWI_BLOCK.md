# WS-C4. Admitting Malawi to the individual-level arms

`data/IPD/Malawi/Malawi_dhs_questionnaire_block.rds`,
`results/tables/wsc4_malawi_block_manifest.csv`.

## The blockage, and that it was not the one recorded

WS-C1 reported Malawi's four cells as `not_computed` for both questionnaire
arms, correctly: the merged dataset carries **one** `gw_` column against 430 to
599 in the other three countries, so a questionnaire arm there would have been
the proxy arm under another label. The previous session recorded this as
requiring a codebook request.

It did not. The micronutrient survey was nested in the 2015-16 DHS. The README
says to merge on MCLUSTER, MNUMBER and M01; MCLUSTER does not exist under that
name, but `gw_cnum` is it. Verified before anything was built:

```
(gw_cnum, mnumber, m01) -> (hv001, hv002, hvidx)     3097 of 3099 rows
```

105 of the DHS's 850 clusters, which is the MNS subsample.

## The guard came first, and it found a twelfth leak class

Every leakage guard in this project matches **analyte names** -- `Hb`,
`ferritin`, `anemia`, `Inflamm`. DHS recodes are named by questionnaire
position, so not one of them fires:

| Column | What it is | Classified by the old guard as |
|:---|:---|:---|
| `gw_hw53` | child haemoglobin, g/dL | **questionnaire** |
| `gw_v453` | women's haemoglobin | **questionnaire** |
| `gw_ha53`, `gw_hc53`, `gw_hb53` | the same in the person recode | **questionnaire** |
| `gw_hml35` | malaria rapid-test result | **questionnaire** |

All of these are present in the Malawi cache. Every one would have entered the
arm whose entire purpose is to exclude blood draws.

`dhs_measurement_class()` closes it, scoped to the survey prefix so DHS-derived
Admin-2 aggregates from other rounds -- which are proxy covariates, available to
a country that never ran a survey -- are untouched. **This is the first instance
of the leak class caught before the data it would have leaked existed.**

Anthropometry is deliberately included as `hb_field` rather than
`questionnaire`: height, weight and their z-scores are not blood draws, but they
are measurements taken from the same respondent at the same visit, which is the
distinction the questionnaire arm rests on.

Six tests pin it, including one that `gw_v459` -- a bednet question sitting
inside the numeric band -- must survive, because an earlier band guard
over-blocked it.

## What was built

| Recode | Population | Matched | Candidates | Kept |
|:---|:---|---:|---:|---:|
| MWHR7ADT household | all | 3097 / 3099 (100%) | 3617 | 137 |
| MWPR7ADT person | all | 3097 / 3099 (100%) | 305 | 96 |
| MWIR7ADT women | women | 809 / 1066 (76%) | 4834 | 240 |
| MWKR7ADT children | children | 1107 / 2033 (54%) | 1027 | 225 |

**450 questionnaire-class columns** after removing 248 duplicated across
recodes, against 430 to 599 in the other three countries. Every column passes
`biomarker_column_class()` as `questionnaire`, asserted on the built object
rather than on a pattern.

## A silent rejection, and the rule that produced it

The first build dropped **every IR and KR column** and said nothing. Completeness
was measured across all in-population rows, but IR reaches only 76 percent of
the MNS women, so an 80 percent threshold exceeded the attainable maximum and
rejected the entire block.

Those are the recodes carrying dietary recall, vitamin A capsule receipt,
deworming and iron in pregnancy -- every nutrition-proximal item the household
block does not have.

Two things were wrong and both are fixed. **Completeness now means item
non-response**, measured among matched rows, and roster coverage is reported
separately as `match_rate`. And **a block that yields nothing now says so**,
naming the threshold it could not reach, because silence read as "the recode was
not there" rather than "the rule was unsatisfiable".

## What this does and does not license

Malawi can now enter the questionnaire and questionnaire-plus-haemoglobin arms,
raising WS-C1 from twelve computable cells to sixteen.

**The instrument is not the same one.** The other three countries' `gw_` blocks
are their own micronutrient survey questionnaires. Malawi's is the DHS
questionnaire, administered to the same households. That is the same *construct*
-- a household questionnaire put to the same people -- but not the same
instrument, and a questionnaire gain computed for Malawi is not strictly
comparable to the other three. The sidecar carries a `wsc4_provenance` attribute
saying so, and the merged dataset is untouched.

**KR matches only 54 percent** of MNS children, because the birth recode covers
births in the last five years to interviewed mothers. That is above the floor
but low enough that Malawi's child arms will rest on about eleven hundred
records.
