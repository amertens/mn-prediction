# Detached runner: propagate uniform BRINDA VitA outcome (DC-H2) through targets.
Sys.setenv(PIPELINE_MODE = "full")
suppressWarnings(suppressMessages(library(targets)))
t0 <- Sys.time()
cat("=== BRINDA regen START", format(t0), "===\n")
ok <- tryCatch({ tar_make_future(workers = 4); TRUE },
               error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); FALSE })
dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
od <- tryCatch(length(tar_outdated()), error = function(e) NA)
cat(sprintf("=== BRINDA regen %s in %s min; outdated remaining=%s ===\n",
            if (ok) "COMPLETE" else "FAILED", dt, od))
