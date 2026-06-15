# Success check for the full corrected run. Exit 0 = success (the comparison
# bundle built with content), 1 = failure/incomplete. Errored slices are logged
# but do not by themselves fail the run unless the bundle is empty.
library(targets)
err <- tryCatch(tar_errored(store = "_targets_corrected"),
                error = function(e) "VERIFY_ERROR")
mc <- tryCatch(tar_read(methods_comparison, store = "_targets_corrected"),
               error = function(e) NULL)
n_slices <- tryCatch(length(mc$slices), error = function(e) 0L)
n_cv <- tryCatch(nrow(mc$cv_compare), error = function(e) 0L)
cat(sprintf("[verify] errored targets (%d): %s\n", length(err),
            paste(err, collapse = ", ")))
cat(sprintf("[verify] methods_comparison: slices=%s cv_compare_rows=%s\n",
            n_slices, n_cv))
ok <- !is.null(mc) && !identical(err, "VERIFY_ERROR") &&
      isTRUE(n_cv > 0)
cat(sprintf("[verify] SUCCESS=%s\n", ok))
quit(status = if (isTRUE(ok)) 0L else 1L)
