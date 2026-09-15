"""
scripts/covariates/source_inventory.py   [AB-01, 2026-09-15]

Source inventory of the shared Admin-2 predictor set at three granularities,
plus the sources that exist only in the per-country (individual-level) merged
datasets. Reproduces the counts quoted in the data-source review.

Granularities
  label    the `source` label in predictors_admin2_shared_metadata.csv
  product  the data product behind each column (a GEE bundle is ~19 products;
           "SoilGrids / iSDA" is two), resolved from the stage-3 data
           dictionary for the base columns and from column prefixes otherwise
  paper    the entry of the West Africa data-landscape review (Tables 1-2) the
           product falls under, or "not in paper"

    python scripts/covariates/source_inventory.py
-> results/tables/source_inventory_<date>.csv          one row per column
-> results/tables/source_inventory_<date>_products.csv one row per product
-> prints the counts
"""
import os, re, datetime, subprocess
import pandas as pd

ROOT = "C:/Users/andre/OneDrive/Documents/mn-prediction"
os.chdir(ROOT)
STAMP = datetime.date.today().isoformat()
H = "data/covariates/harmonized"

S = pd.read_csv(f"{H}/predictors_admin2_shared.csv", low_memory=False)
M = pd.read_csv(f"{H}/predictors_admin2_shared_metadata.csv")
D = pd.read_csv(f"{H}/data_dictionary.csv")
cols = [c for c in S.columns if c not in ("country", "Admin1", "Admin2")]
assert set(cols) == set(M.column), "metadata and data disagree"
COUNTRIES = sorted(S.country.unique())

# ── product resolution ──────────────────────────────────────────────────────
fam = dict(zip(D.canonical, D.family)); prov = dict(zip(D.canonical, D.provider))
FAMILY_PRODUCT = {  # stage-3 families -> (product, provider, paper entry)
    "alphaearth": ("AlphaEarth Foundations 2017 embedding", "Google DeepMind (GEE)", "not in paper (GEE catalogue)"),
    "dhs": ("DHS recodes, surveyPrev Admin-2 (BYM2)", "DHS Program", "T1 Demographic and Health Surveys"),
    "soil": ("iSDAsoil 30 m soil properties", "iSDA", "not in paper"),
    "soilgrids": ("SoilGrids v2 0-5 cm", "ISRIC", "T2 ISRIC SoilGrids"),
    "mapspam": ("MapSPAM 2010 crop allocation", "IFPRI", "not in paper"),
    "trmm": ("TRMM 3B43 precipitation", "NASA", "T2 NASA EarthData / GEE"),
    "terraclimate": ("TerraClimate", "University of Idaho", "T2 Google Earth Engine"),
    "lst_night": ("Oxford MAP night LST (MODIS)", "Oxford MAP", "T2 Malaria Atlas Project"),
    "ndvi": ("AVHRR NDVI CDR v5", "NOAA", "T2 NASA EarthData / GEE"),
    "ndvi_mean_anomaly": ("FEWS NET NDVI anomaly", "FEWS NET / USGS", "T2 FEWS NET"),
    "dailyevi": ("MODIS MOD13Q1 EVI", "NASA MODIS", "T2 NASA EarthData / GEE"),
    "lai8days": ("MODIS MOD15A2H LAI", "NASA MODIS", "T2 NASA EarthData / GEE"),
    "productivity": ("MODIS MOD17A3 NPP", "NASA MODIS", "T2 NASA EarthData / GEE"),
    "wapor": ("FAO WaPOR NPP", "FAO", "T2 Google Earth Engine"),
    "landcoverlayers": ("Copernicus Global Land Cover", "Copernicus", "T2 Copernicus / GEE"),
    "gpw_grasslands": ("Global Pasture Watch grasslands", "Global Pasture Watch", "not in paper (GEE catalogue)"),
    "popdensity": ("GPW v4.11 population density", "CIESIN", "T2 NASA EarthData / GEE"),
    "ghsbuilts": ("GHSL built-up surface", "JRC", "not in paper (GEE catalogue)"),
    "ghspop": ("GHSL population", "JRC", "not in paper (GEE catalogue)"),
    "wsf": ("World Settlement Footprint 2015", "DLR", "not in paper (GEE catalogue)"),
    "ccnl": ("CCNL night-time lights", "BNU", "not in paper (GEE catalogue)"),
    "globalhumanmodification": ("Global Human Modification", "CSP", "not in paper (GEE catalogue)"),
    "elevation": ("SRTM elevation", "USGS", "T2 NASA EarthData / GEE"),
    "aerosoloptical": ("MODIS MOD08 aerosol optical depth", "NASA MODIS", "T2 NASA EarthData / GEE"),
    "accessibility": ("Oxford MAP travel time to healthcare", "Oxford MAP", "T2 Malaria Atlas Project"),
}
PREFIX_PRODUCT = [  # (regex, product, provider, paper entry)
    (r"^dhs_", "DHS recodes, surveyPrev Admin-2 (BYM2)", "DHS Program", "T1 Demographic and Health Surveys"),
    (r"^aef_", "AlphaEarth Foundations 2017 embedding", "Google DeepMind (GEE)", "not in paper (GEE catalogue)"),
    (r"^soilgrids_|^sg_", "SoilGrids v2 0-5 cm", "ISRIC", "T2 ISRIC SoilGrids"),
    (r"^soil_", "iSDAsoil 30 m soil properties", "iSDA", "not in paper"),
    (r"^ihme_", "IHME LBD geospatial surfaces (13 families)", "IHME", "T2 IHME GHDx"),
    (r"^map_", "Malaria Atlas Project rasters at survey year", "Malaria Atlas Project", "T2 Malaria Atlas Project"),
    (r"^spam_", "MapSPAM 2010 crop allocation", "IFPRI", "not in paper"),
    (r"^fprice_", "WFP Food Price Database (via HDX HAPI)", "WFP / OCHA HDX", "T2 WFP Food Price Database"),
    (r"^fao_", "FAOSTAT Food Balance Sheets", "FAO", "T2 FAOSTAT / FBS / SUA"),
    (r"^fsec_", "HFID (IPC / CH / FEWS NET / WFP mVAM)", "FEWS NET / IPC / WFP", "T2 HFID / IPC / FEWS NET"),
    (r"^(koppen|kg_|aez|zone_)", "Koppen-Geiger 1991-2020 + HarvestChoice AEZ16", "Beck et al. / IFPRI", "not in paper"),
    (r"^glw_", "Gridded Livestock of the World 4 (2020)", "FAO", "not in paper"),
    (r"^wdist_|^coast_|^water_", "JRC Global Surface Water + LSIB coastline", "JRC / US DoS", "not in paper (GEE catalogue)"),
    (r"^espen_", "WHO ESPEN implementation-unit database", "WHO AFRO", "T2 WHO ESPEN"),
    (r"^wpop_|^smod_|^ghsl_smod", "WorldPop age-sex + GHS-SMOD", "WorldPop / JRC", "T2 WorldPop"),
    (r"^rwi_", "Meta Relative Wealth Index", "Meta Data for Good", "not in paper"),
    (r"^gpw2010_|^popdens_", "GPW v4.11 age-sex / density", "CIESIN", "T2 NASA EarthData / GEE"),
    (r"^vas_", "UNICEF vitamin A supplementation coverage (WDI mirror)", "UNICEF / World Bank", "T2 UNICEF Global Databases"),
    (r"^flunet_", "WHO FluNet", "WHO", "T1 FluNet"),
    (r"^gfdx_(anemia|zinc_def|ntd)", "GFDx nutrition-status fields (WHO 2011 anaemia, W&B zinc)", "GFDx", "T2 GFDx"),
    (r"^gfdx_", "GFDx fortification programme fields", "GFDx", "T2 GFDx"),
    (r"^who_anaemia_", "WHO Global Anaemia Estimates", "WHO / UNICEF", "T2 WHO Global Anaemia Estimates"),
    (r"^mimi_", "Tang et al. 2026 (WFP / MIMI) nutrient inadequacy", "WFP / MIMI", "T2 MIMI"),
    (r"^gdl_", "Global Data Lab subnational HDI", "GDL", "not in paper"),
    (r"^mics_heat_", "MICS by region via WHO Health Inequality Data Repository", "WHO / UNICEF", "T1 MICS via T2 WHO Health Inequality Data Repository"),
    (r"^acled_", "ACLED conflict events", "ACLED", "T1 ACLED"),
    (r"^(ndvi|evi|lai|npp|wapor|grassland|tclim|precip|lst|aod|ghs|built|wsf|ntl|human|lcover|elevation|access)", "GEE raster (family unresolved)", "GEE", "T2 Google Earth Engine"),
]

def resolve(col):
    f = fam.get(col)
    if f in FAMILY_PRODUCT:
        return FAMILY_PRODUCT[f]
    for rx, prod, pv, paper in PREFIX_PRODUCT:
        if re.search(rx, col):
            return (prod, pv, paper)
    return ("UNRESOLVED", "", "")

R = M.copy()
R[["product", "provider", "paper_entry"]] = pd.DataFrame([resolve(c) for c in R.column], index=R.index)
fin = {cc: [float((S.loc[S.country == cc, c].notna()).mean()) for c in R.column] for cc in COUNTRIES}
for cc in COUNTRIES: R[f"finite_{cc}"] = fin[cc]
R["n_countries_finite"] = sum((R[f"finite_{cc}"] > 0.5).astype(int) for cc in COUNTRIES)
R["subnational"] = R.subnational.astype(str).str.upper().eq("TRUE")
R.to_csv(f"results/tables/source_inventory_{STAMP}.csv", index=False)

P = (R.groupby(["product", "provider", "paper_entry"])
       .agg(columns=("column", "size"), all4=("n_countries_finite", lambda x: int((x == len(COUNTRIES)).sum())),
            subnational=("subnational", "sum"), national_constant=("subnational", lambda x: int((~x).sum())))
       .reset_index().sort_values("columns", ascending=False))
P.to_csv(f"results/tables/source_inventory_{STAMP}_products.csv", index=False)

print(f"shared set: {S.shape[0]} rows x {len(cols)} predictors, countries {COUNTRIES}")
print(f"source labels: {M.source.nunique()} | products: {P['product'].nunique()} | providers: {P.provider.nunique()}")
print(f"paper entries touched: {P.paper_entry[~P.paper_entry.str.startswith('not in paper')].nunique()} | products not in the paper: {int((P.paper_entry.str.startswith('not in paper')).sum())}")
print(f"unresolved columns: {int((R['product'] == 'UNRESOLVED').sum())}")
print(f"columns present in all {len(COUNTRIES)} countries: {int((R.n_countries_finite == len(COUNTRIES)).sum())}; national constants: {int((~R.subnational).sum())}")
print()
print(P.to_string(index=False))

# ── per-country (individual-level) datasets: domains not in the shared set ──
print("\n=== per-country merged datasets: domain column counts (metadata/<country>_variable_categories.rds) ===")
rs = 'for (cn in c("gambia","ghana","SL","malawi")) { x <- unclass(readRDS(file.path("metadata", paste0(cn, "_variable_categories.rds")))); cat(cn, ":", paste(sprintf("%s=%d", sub("_vars$", "", names(x)), sapply(x, length)), collapse=" "), "\\n") }'
out = subprocess.run(["C:/Program Files/R/R-4.4.2/bin/Rscript.exe", "-e", rs], capture_output=True, text=True)
print(out.stdout.strip() or out.stderr[-500:])
