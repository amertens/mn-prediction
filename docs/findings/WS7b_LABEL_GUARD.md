# WS7b. The label-derived assay guard

Section 7.5 states the residual risk and the correct fix: the merged data carries
no variable labels, so `is_biomarker_column()` is a regex over opaque names whose
failure mode is silent, and the recommended action is to obtain the original
Stata files and re-derive the guard from labels.

That is what this workstream does. WS7a had already shown the failure mode is
real rather than hypothetical.

## Coverage

| Country | Columns in the `.dta` files | Labelled | Carried by the pipeline |
|:---|---:|---:|---:|
| Gambia | 812 | 710 | 575 |
| Ghana | 586 | 345 | 367 |
| Sierra Leone | 752 | 189 | 446 |
| **Malawi** | **none** | **none** | not coverable |

**Malawi cannot be covered and this is not a small gap.** No `.dta` exists. Its
source is `data/IPD/Malawi/Malawi_merged_dataset.rds`, whose 242 questionnaire
columns are coded `m01`, `m115a`, `m220h` with **zero** labels, and the only
local documentation is a 349-character README pointing to `immpact@cdc.gov`. The
open request is for the Malawi MNS questionnaire codebook, or for the DHS merge
the README describes, which would supply labelled equivalents.

## Result

Of **865 labelled columns the pipeline carries**, **828 agree** with the
name-based guard (source: `results/tables/assay_guard_label_comparison.csv`,
column `agreement`).

| Agreement | Count |
|:---|---:|
| agree | 828 |
| LABEL BLOCKS, NAME ALLOWS | 22 |
| NAME BLOCKS, LABEL ALLOWS | 12 |
| different class | 3 |

### The 22 the label blocks and the name allowed

**Thirteen are genuine and are now blocked.** Each was adjudicated against its
label, and the pattern added was verified against all 4,906 columns to match
exactly these and no others:

| Column | Label | Why it is a draw |
|:---|:---|:---|
| `gw_cInflammMarker`, `gw_wInflammMarker` | "CRP and AGP elevation" | computed from the assay panel |
| `gw_cInflammYN`, `gw_wInflammYN` | "Elevated CRP or AGP" | same |
| `gw_cAnemMalarInflamm`, `gw_wAnemMalarInflamm` | "Status of anemia, malaria, and inflamm" | same |
| six Gambia hyperglycaemia columns | "RECODE of bs3_hba1c (Glycated hemoglobin (%))" | derived from HbA1c |
| `gw_cpcn` | "16. Phlebotomist's ID number" | the draw's metadata |

**None of them names an analyte**, so no pattern over analyte names could have
caught any of them. This is the same lesson Section 7.2 records and it has now
recurred twice more.

`gw_wAnemMalarInflamm` was the **top-ranked eligible column for Sierra Leone
women's iron in the WS7a leakage report at |r| 0.358**, and WS7a flagged it there
for label verification without asserting what it was. The label confirms it. The
correlation report and the label guard found the same column independently,
which is the intended behaviour of having both.

**Nine are false positives of the label predicate and are correctly left
unblocked:** the Gambia `gw_wAnemCauses_*` and `gw_wAnemPrev_*` knowledge blocks
("Anemia prevented by eating iron rich food??"), `gw_nkaq2` and `gw_nkaq3`
("What are the causes of anemia?"), and `gw_cREOil_day` ("Retinol Equiv (RE)
consumer per day from vegetable oil"), which is dietary intake rather than a
measurement of the respondent.

**One is a confirmation worth recording.** Section 7.5 lists `gw_wFFAnemia` as an
assumption: it is *assumed* to be a fortified-food belief item on the strength of
its `gw_wFF*` siblings, and Section 7.5 says that if it is actually a measured
status it should be blocked. Its label reads **"Woman says Fort Flour reduces
anemia"**. The assumption is correct and can stop being an assumption.

### The 12 the name blocks and the label allowed

Eight are correctly blocked despite uninformative labels: the Gambia
`gw_bs3`-`gw_bs6m` block ("Tube Filling", "Malaria RDT result", "Time of blood
taking") is draw metadata, and `gw_wCRPHigh`, `gw_wAGPHigh`, `gw_cCRPHigh`,
`gw_cAGPHigh` ("High is 5.0+", "High is 1.0+") are CRP and AGP elevation flags.

**Two are genuine over-blocking.** Gambia's `gw_cID` and `gw_wID` are labelled
"Child ID number constructed using cluster #, HH #, and child line #" and the
equivalent for women. The guard's `_[cw]ID($|_|s[TR])` pattern was written to
catch iron-deficiency indicators and Section 7.3 states it does exactly that.
In Gambia those two columns are identifiers. They are left blocked, because the
same name may denote iron deficiency in another survey and blocking an
identifier costs nothing as a predictor, but the case is recorded: the guard is
fragile in both directions and the label is the only thing that can tell which.

## What this changes about the guard's standing

The guard now blocks 90 of the 4,906 columns, against 77 before WS7a and WS7b.
It has been checked against ground truth for three of the four surveys, and 828
of 865 labelled columns agree. That is a materially stronger position than
Section 7.5 describes.

It is not a solved problem. The label predicate is itself a regex over English
text and it produced nine false positives on 865 columns. The 137 to 189 columns
in Sierra Leone's files that carry no label are unchecked. Malawi is entirely
unchecked and contributes 242 opaque columns that no automated method can
classify.

## Reviewer pass, statistical

**No model is fitted.** This workstream compares two classifiers of column
names against each other and against text. The statistical risk it addresses is
leakage, and the adjudication of all 22 disagreements is set out above rather
than summarised, so a reviewer can disagree with any individual call.

**Scope of the claim.** "828 of 865 agree" is over columns that are both labelled
and carried by the pipeline. It is not a statement about the 4,906 columns, most
of which are external covariates with no Stata label at all.

**The adjudication is mine, not the label's.** Nine of the 22 were judged false
positives on a reading of the label text. Another reader could take
`gw_cREOil_day` differently, though "consumer per day from vegetable oil" is
plainly a dietary quantity.

## Reviewer pass, reproducibility

**Targets graph.** WS7a's `leakage_report` target covers the re-verification
hook the workstream asks for: it depends on every `outcome_data_*` target, so a
newly ingested country cannot reach a model without its columns being ranked.
The label comparison is a script rather than a target because its inputs are the
raw `.dta` files, which are not pipeline targets.

**Stamp targets.** Reads six `.dta` files under `data/IPD/`. They are tracked
inputs. The script names each and reports missing ones rather than failing.

**Determinism.** `haven::read_dta(n_max = 1)` reads only the header, so the run
takes seconds and depends on no stochastic call.

**Guard verification.** After each pattern addition the guard was re-checked
against all 4,906 columns for the count of newly blocked columns, and every one
was listed and matched to its label before the change was kept. The household
block and every exposure and identifier in Section 7.4 were re-confirmed
unblocked.
