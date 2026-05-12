# =============================================================================
# National Prediction Pipeline — Configuration
# =============================================================================

#' Get data paths for the national prediction pipeline
get_data_paths <- function() {
  list(
    # National-level proxy predictors (DHS, IHME, etc.)
    predictor_data = "C:/Users/andre/OneDrive/Documents/mn-proxies/data/national/combined_dataset.dta",

    # BRINDA micronutrient survey prevalence (in mn-proxies repo)
    brinda_data = "C:/Users/andre/OneDrive/Documents/mn-proxies/data/clean_brinda_data.rds",

    # WHO VMNIS indicator database (in mn-proxies repo)
    vmnis_data = "C:/Users/andre/OneDrive/Documents/mn-proxies/data/VMNIS/VMNISIndicator_long_format.dta"
  )
}


#' Get BRINDA outcome configurations
#'
#' Each entry defines an outcome variable and the SL stack to use.
#' prev_rol is commented out in the legacy script so excluded here.
get_brinda_outcomes <- function() {
  list(
    prev_ferr = list(
      outcome  = "prev_ferr",
      label    = "Iron deficiency (ferritin)",
      sl_stack = "full"
    ),
    prev_stfr = list(
      outcome  = "prev_stfr",
      label    = "Iron deficiency (sTfR)",
      sl_stack = "full"
    ),
    prev_rbp = list(
      outcome  = "prev_rbp",
      label    = "Vitamin A deficiency (RBP)",
      sl_stack = "full"
    ),
    prev_zinc = list(
      outcome  = "prev_zinc",
      label    = "Zinc deficiency",
      sl_stack = "full"
    )
  )
}


#' Get VMNIS outcome configurations
#'
#' Each entry defines: VMNIS indicator name, metric (Mean or Prevalenceofdeficiency),
#' and the SL stack. B12 and VitD use sl_simple because slfull fails on these.
get_vmnis_outcomes <- function() {
  list(
    # --- Zinc ---
    mean_zinc = list(
      indicator = "Zinc (plasma or serum)",
      metric    = "Mean",
      label     = "Zinc — Mean",
      sl_stack  = "full"
    ),
    prev_zinc = list(
      indicator = "Zinc (plasma or serum)",
      metric    = "Prevalenceofdeficiency",
      label     = "Zinc — Prevalence",
      sl_stack  = "full"
    ),

    # --- Vitamin B12 (sl_simple: slfull fails) ---
    mean_b12 = list(
      indicator = "Vitamin B12",
      metric    = "Mean",
      label     = "Vitamin B12 — Mean",
      sl_stack  = "simple"
    ),
    prev_b12 = list(
      indicator = "Vitamin B12",
      metric    = "Prevalenceofdeficiency",
      label     = "Vitamin B12 — Prevalence",
      sl_stack  = "simple"
    ),

    # --- Retinol ---
    mean_retinol = list(
      indicator = "Retinol (plasma or serum)",
      metric    = "Mean",
      label     = "Retinol — Mean",
      sl_stack  = "full"
    ),
    prev_retinol = list(
      indicator = "Retinol (plasma or serum)",
      metric    = "Prevalenceofdeficiency",
      label     = "Retinol — Prevalence",
      sl_stack  = "full"
    ),

    # --- Vitamin D (sl_simple: slfull fails) ---
    mean_vitD = list(
      indicator = "25-Hydroxyvitamin D",
      metric    = "Mean",
      label     = "Vitamin D — Mean",
      sl_stack  = "simple"
    ),
    prev_vitD = list(
      indicator = "25-Hydroxyvitamin D",
      metric    = "Prevalenceofdeficiency",
      label     = "Vitamin D — Prevalence",
      sl_stack  = "simple"
    ),

    # --- Folate ---
    mean_folate = list(
      indicator = "Folate (plasma or serum)",
      metric    = "Mean",
      label     = "Folate — Mean",
      sl_stack  = "full"
    ),
    prev_folate = list(
      indicator = "Folate (plasma or serum)",
      metric    = "Prevalenceofdeficiency",
      label     = "Folate — Prevalence",
      sl_stack  = "full"
    ),

    # --- Haemoglobin (mean only, no prevalence; sl_simple due to N=27k rows) ---
    mean_haemoglobin = list(
      indicator = "Haemoglobin",
      metric    = "Mean",
      label     = "Haemoglobin — Mean",
      sl_stack  = "simple"
    )
  )
}


#' Get SL fitting parameters
get_sl_params <- function() {
  list(
    seed  = 12345L,
    folds = 2L,
    CV    = FALSE,
    prescreen_pval = 0.2
  )
}


#' LMIC country list for VMNIS heatmap filtering
get_lmic_countries <- function() {
  c(
    # Africa
    "Malawi", "Mozambique", "Zambia", "Zimbabwe", "Angola", "Benin",
    "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
    "Chad", "Comoros", "Congo", "Côte d'Ivoire",
    "Democratic Republic of the Congo", "Djibouti", "Equatorial Guinea",
    "Eritrea", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea",
    "Guinea-Bissau", "Kenya", "Lesotho", "Liberia", "Madagascar", "Mali",
    "Mauritania", "Mauritius", "Namibia", "Niger", "Nigeria", "Rwanda",
    "Sao Tome and Principe", "Senegal", "Seychelles", "Sierra Leone",
    "Somalia", "South Africa", "South Sudan", "Sudan", "Swaziland",
    "Tanzania", "Togo", "Uganda",
    # Asia
    "Afghanistan", "Bangladesh", "Bhutan", "Cambodia", "China",
    "Democratic People's Republic of Korea", "India", "Indonesia",
    "Iran (Islamic Republic of)", "Iraq", "Jordan", "Kiribati",
    "Kyrgyzstan", "Lao People's Democratic Republic", "Lebanon",
    "Malaysia", "Maldives", "Marshall Islands",
    "Micronesia (Federated States of)", "Mongolia", "Myanmar", "Nauru",
    "Nepal", "Pakistan", "Palau", "Papua New Guinea", "Philippines",
    "Samoa", "Solomon Islands", "Sri Lanka", "Tajikistan", "Thailand",
    "Timor-Leste", "Tonga", "Turkmenistan", "Tuvalu", "Uzbekistan",
    "Vanuatu", "Viet Nam", "Yemen",
    # Europe / Middle East
    "Albania", "Armenia", "Azerbaijan", "West Bank and Gaza Strip",
    "Syrian Arab Republic",
    # Americas
    "Belize", "Bolivia (Plurinational State of)", "Brazil", "Chile",
    "Colombia", "Costa Rica", "Dominican Republic", "Ecuador",
    "El Salvador", "Guatemala", "Guyana", "Haiti", "Jamaica", "Mexico",
    "Nicaragua", "Panama", "Paraguay", "Peru",
    "Venezuela (Bolivarian Republic of)", "Argentina",
    # North Africa
    "Algeria", "Botswana", "Egypt", "Kazakhstan", "Morocco", "Tunisia",
    "United Republic of Tanzania"
  )
}
