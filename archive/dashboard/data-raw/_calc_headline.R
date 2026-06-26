suppressWarnings({
d   <- readRDS("dashboard/data/transport_calibration.rds")$adm2_area
pop <- readRDS("dashboard/data/admin2_population.rds")
lab <- c(gambia="Gambia", ghana="Ghana", sierraleone="Sierra Leone", malawi="Malawi")
d$country_label <- lab[d$country]
cat(sprintf("%-12s %-12s %4s %6s %6s %8s %6s\n",
            "outcome","scope","n","r","MAE","reach20","lift"))
for (oc in unique(d$outcome)) {
  for (scope in c("POOLED", unique(d$country_label))) {
    sub <- if (scope == "POOLED") d[d$outcome == oc, ]
           else d[d$outcome == oc & d$country_label == scope, ]
    sub <- sub[is.finite(sub$true_prev) & is.finite(sub$pred_prev), ]
    if (nrow(sub) < 3) next
    isch <- startsWith(oc, "child_"); pc <- if (isch) "pop_child" else "pop_women"
    popk <- data.frame(country_label = pop$country, Admin2 = pop$Admin2, w = pop[[pc]])
    popk <- popk[!duplicated(popk[, c("country_label","Admin2")]), ]
    m <- merge(sub, popk, by = c("country_label","Admin2"), all.x = TRUE, sort = FALSE)
    w <- m$w; w[!is.finite(w) | w <= 0] <- NA
    if (all(is.na(w))) w <- rep(1, nrow(m)) else w[is.na(w)] <- median(w, na.rm = TRUE)
    burden <- m$true_prev * w; totB <- sum(burden); totW <- sum(w)
    om <- order(m$pred_prev, decreasing = TRUE)
    cx <- c(0, cumsum(w[om])/totW); cy <- c(0, cumsum(burden[om])/totB)
    reach <- approx(cx, cy, xout = 0.2, rule = 2)$y
    r <- cor(m$pred_prev, m$true_prev); mae <- mean(abs(m$pred_prev - m$true_prev))*100
    cat(sprintf("%-12s %-12s %4d %6.2f %6.1f %7.0f%% %5.2fx\n",
                oc, scope, nrow(m), r, mae, reach*100, reach/0.2))
  }
}
})
