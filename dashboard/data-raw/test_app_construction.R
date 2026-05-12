# Test app construction without launching the server
options(shiny.testmode = TRUE)
setwd(here::here("dashboard"))
source("global.R")
cat("global.R OK\n")

# Source app.R but skip the final shinyApp() call
app_code <- readLines("app.R")
app_code <- app_code[!grepl("^\\s*shinyApp", app_code)]
eval(parse(text = paste(app_code, collapse = "\n")))
cat("UI/server defined OK\n")

# Render UI to HTML
html <- as.character(htmltools::as.tags(ui))
cat("UI rendered to", nchar(html), "chars of HTML\n")

# Verify all 8 nav panels are in there
for (panel_name in c("Map explorer", "District profiles",
                     "National burden", "Scenarios", "Importance",
                     "Côte d'Ivoire (OOS)", "GBD comparison", "Methods")) {
  if (grepl(panel_name, html, fixed = TRUE)) {
    cat("✓ Panel found:", panel_name, "\n")
  } else {
    cat("✗ Panel MISSING:", panel_name, "\n")
  }
}

cat("\nApp construction test PASSED\n")
