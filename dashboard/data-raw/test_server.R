# Test that the Shiny server starts cleanly
options(shiny.host = "127.0.0.1", shiny.port = 4567L)
setwd(here::here("dashboard"))

later::later(function() {
  cat("\n=== Stopping server after grace period ===\n")
  shiny::stopApp()
}, delay = 5)

shiny::runApp(".", launch.browser = FALSE)
cat("\n=== Server stopped cleanly ===\n")
