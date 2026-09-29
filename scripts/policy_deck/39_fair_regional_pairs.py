"""The survey's regional average, scored fairly, for the per-nutrient ranking slide.

The regional-average arm of script 28 (region_mean_jk) gives each hidden district the mean
of the OTHER surveyed districts in its region. Two surveyed districts in the same region
then come out in reversed order by construction: the higher one's figure leaves out the
high value, the lower one's keeps it. A planner's real regional figure gives every
district in a region the same value, a tie. Scored that way (region figure from the
other districts, same-region pairs and tied figures counted as half a pair, i.e. a coin
toss), the comparator is fair. The survey's published regional figure scored against
each district's own value is NOT used: it contains that district's own data.

Found by the 28 September review session (pre-fix: model 67%, regional as plotted 62%,
fair 66%); this recomputes it on the post-fix targets for the same cells as script 34.

    python scripts/policy_deck/39_fair_regional_pairs.py
<- results/tables/protocol_v2/targets_v2.csv, results/tables/policy_deck/v3_percell_pairs.csv,
   results/figures/mnf15/cell_master.csv
-> results/tables/policy_deck/v6_fair_regional_pairs.csv
"""
import itertools

import numpy as np
import pandas as pd

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction/"
TG = pd.read_csv(ROOT + "results/tables/protocol_v2/targets_v2.csv")
CM = pd.read_csv(ROOT + "results/figures/mnf15/cell_master.csv")
PP = pd.read_csv(ROOT + "results/tables/policy_deck/v3_percell_pairs.csv")
PP = PP[PP.estimand == "infill"]
keep = set(map(tuple, CM[CM.keep][["country", "outcome"]].values))

rows = []
for (cn, on), g in TG.groupby(["country", "outcome"]):
    if (cn, on) not in keep or cn == "SierraLeone":   # as script 34: no in-country test for Sierra Leone
        continue
    g = g[np.isfinite(g.y_level) & np.isfinite(g.n_eff_cont)]
    if len(g) < 12 or g.Admin1.nunique() < 3:
        continue
    y, a, n = g.y_level.values, g.Admin1.values, len(g)
    loo = np.array([y[(a == a[i]) & (np.arange(n) != i)].mean() if (a == a[i]).sum() > 1 else np.delete(y, i).mean()
                    for i in range(n)])
    agree = decided = ties = 0
    for i, j in itertools.combinations(range(n), 2):
        do = y[i] - y[j]
        if do == 0:
            continue
        dp = 0.0 if a[i] == a[j] else loo[i] - loo[j]
        if dp == 0:
            ties += 1
            continue
        agree += (do > 0) == (dp > 0)
        decided += 1
    r = PP[(PP.country == cn) & (PP.outcome == on)].set_index("arm").pairs
    rows.append(dict(country=cn, outcome=on, n=n, share_same_region=np.mean([a[i] == a[j] for i, j in itertools.combinations(range(n), 2)]),
                     index=r.get("domain_index"), regional_as_plotted=r.get("region_mean_jk"),
                     regional_fair=(agree + 0.5 * ties) / (decided + ties)))
R = pd.DataFrame(rows)
assert len(R) == 14, f"expected the 14 in-country combinations of script 34, got {len(R)}"
R.to_csv(ROOT + "results/tables/policy_deck/v6_fair_regional_pairs.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(3).to_string(index=False))
print(R[["index", "regional_as_plotted", "regional_fair"]].mean().round(3).to_string())
print("model ahead of the fair regional average in", int((R["index"] > R.regional_fair).sum()), "of", len(R))
