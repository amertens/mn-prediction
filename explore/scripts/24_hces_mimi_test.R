# =============================================================================
# explore/scripts/24_hces_mimi_test.R   [probe HC-03]
# Does the HCES/MIMI dietary block predict the BIOMARKER, nutrient-matched?
# The store already carries MIMI's own modelled inadequate-intake percentages
# (Tang et al. 2026) plus 15 own-derived HCES consumption indicators. If the
# intake construct helps, it should help MOST when matched to its own nutrient.
# Same mismatch control as MX-01.
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
MD <- E$MD
HCES <- MD$column[grepl("HCES", MD$source) | grepl("HCES", MD$domain)]
HCES <- intersect(HCES, E$PREDS)
MIMI <- grep("^mimi_", HCES, value = TRUE)
OWN  <- setdiff(HCES, MIMI)
cat("HCES block:", length(HCES), "columns |", length(MIMI), "MIMI modelled intake,",
    length(OWN), "own-derived consumption indicators\n")
nut_of <- function(o) if (grepl("zinc",o)) "zinc" else if (grepl("iron",o)) "iron" else
  if (grepl("vitA",o)) "vita" else if (grepl("folate",o)) "folate" else
  if (grepl("b12",o)) "b12" else "other"
MIS <- c(zinc="vita", vita="zinc", iron="folate", folate="iron", b12="zinc", other="vita")
mimi_for <- function(n) grep(paste0("^mimi_", n, "_"), MIMI, value = TRUE)
ridge_on <- function(pick) function(tr,te,y,X,D,aux){
  cc <- intersect(pick, colnames(X)); if(length(cc)<1) return(rep(mean(y[tr]),length(te)))
  if(length(cc)==1){b<-stats::lm.fit(cbind(1,X[tr,cc,drop=FALSE]),y[tr])$coefficients
    b[!is.finite(b)]<-0; return(as.numeric(cbind(1,X[te,cc,drop=FALSE])%*%b))}
  .v2_enet(X[tr,cc,drop=FALSE],y[tr],X[te,cc,drop=FALSE],alpha=0)}
idx_plus <- function(pick) function(tr,te,y,X,D,aux){
  cc<-intersect(pick,colnames(X)); D2<-if(length(cc)) cbind(D,X[,cc,drop=FALSE]) else D
  arm_domain_index_v2(tr,te,y,X,D2,aux)}
rows <- list(); ix <- exp_cell_index(E)
for (i in seq_len(nrow(ix))) {
  n <- nut_of(ix$outcome[i]); mm <- mimi_for(n); mx <- mimi_for(MIS[[n]])
  arms <- c(exp_baseline_arms()[c("null_train_mean","domain_index")],
            list(mimi_matched = ridge_on(mm), mimi_mismatched = ridge_on(mx),
                 hces_own = ridge_on(OWN), hces_all = ridge_on(HCES),
                 index_plus_hces = idx_plus(HCES)))
  for (tg in c("level","prev")) {
    cell <- tryCatch(exp_cell(E, ix$country[i], ix$outcome[i], tg), error=function(e) NULL)
    if (is.null(cell)) next
    r <- exp_infill(cell, arms, reps=REPS); r$nutrient <- n; r$n_mimi <- length(mm)
    rows[[paste(i,tg)]] <- r
  }
  message("  ", ix$country[i], " ", ix$outcome[i], " (", n, ", ", length(mm), " MIMI cols)")
}
SM <- exp_summarise(dplyr::bind_rows(rows)); exp_write(SM, "24_hces_mimi")
cat("\n== in-fill, level: median Spearman ==\n")
a <- SM[SM$estimand=="infill" & SM$target=="level", ]
print(aggregate(spearman ~ arm, data=a, FUN=function(z) round(median(z,na.rm=TRUE),3)), row.names=FALSE)
cat("\n== paired against the index, blocks = country ==\n")
w <- reshape(a[,c("country","outcome","arm","spearman")], idvar=c("country","outcome"),
             timevar="arm", direction="wide"); names(w) <- sub("^spearman[.]","",names(w))
for (ar in c("index_plus_hces","hces_all","hces_own","mimi_matched")) {
  g <- w[[ar]] - w$domain_index; ok <- is.finite(g)
  b <- aggregate(list(g=g[ok]), by=list(b=w$country[ok]), FUN=mean)
  cat(sprintf("  %-16s mean %+.4f  cells %2d/%2d  countries %d/%d\n", ar, mean(g[ok]),
      sum(g[ok]>0), sum(ok), sum(b$g>0), nrow(b)))
}
cat("\n== MIMI nutrient specificity: matched vs mismatched ==\n")
g <- w$mimi_matched - w$mimi_mismatched; ok <- is.finite(g)
cat(sprintf("  matched better in %d of %d cells, mean %+.4f, sign p=%.3f\n",
    sum(g[ok]>0), sum(ok), mean(g[ok]), binom.test(sum(g[ok]>0), sum(ok), 0.5)$p.value))
