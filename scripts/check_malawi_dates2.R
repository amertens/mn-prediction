d <- readRDS("data/IPD/Malawi/Malawi_merged_dataset.rds")
dates <- d$date_interview
dates <- dates[!is.na(dates)]
cat("Interview date range:", format(min(dates)), "to", format(max(dates)), "\n")
cat("N interviews:", length(dates), "\n\n")
cat("By month:\n")
print(table(format(dates, "%Y-%m")))
