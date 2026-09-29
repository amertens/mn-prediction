# =============================================================================
# explore/scripts/27_mimi_mw_score.R   [HC-06, task 4 on the Malawi MIMI block]
# Does MIMI's Malawi district inadequacy predict the BIOMARKER, nutrient-matched?
# MX-01's matched-vs-mismatched control throughout.
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
REPS <- as.integer(Sys.getenv("EXP_REPS","10"))
E <- exp_load()
G <- read.csv(file.path(EXP_OUT,"26_mimi_malawi_admin2.csv"), stringsAsFactors=FALSE)
PCT <- grep("_inadequate_pct$", names(G), value=TRUE)
G2 <- G[, c("country","Admin1","Admin2", PCT)]
cat("MIMI Malawi block:", length(PCT), "prevalence columns,", nrow(G2), "rows\n")
nut <- function(o) if(grepl("zinc",o)) "zinc" else if(grepl("iron",o)) "iron" else
  if(grepl("vitA",o)) "vita" else if(grepl("folate",o)) "folate" else
  if(grepl("b12",o)) "b12" else if(grepl("selenium",o)) "selenium" else "other"
MIS <- c(zinc="vita", vita="zinc", iron="folate", folate="iron", b12="zinc",
         selenium="iron", other="vita")
col_for <- function(n) grep(paste0("^mimiMW_",n,"_inadequate_pct$"), PCT, value=TRUE)
ridge <- function(pick) function(tr,te,y,X,D,aux){
  cc <- intersect(pick, colnames(X)); if(!length(cc)) return(rep(mean(y[tr]),length(te)))
  if(length(cc)==1){b<-stats::lm.fit(cbind(1,X[tr,cc,drop=FALSE]),y[tr])$coefficients
    b[!is.finite(b)]<-0; return(as.numeric(cbind(1,X[te,cc,drop=FALSE])%*%b))}
  .v2_enet(X[tr,cc,drop=FALSE],y[tr],X[te,cc,drop=FALSE],alpha=0)}
iplus <- function(pick) function(tr,te,y,X,D,aux){
  cc<-intersect(pick,colnames(X)); D2<-if(length(cc)) cbind(D,X[,cc,drop=FALSE]) else D
  arm_domain_index_v2(tr,te,y,X,D2,aux)}
ix <- exp_cell_index(E); ix <- ix[ix$country=="Malawi",]
rows <- list()
for (i in seq_len(nrow(ix))) {
  n <- nut(ix$outcome[i]); mm <- col_for(n); mx <- col_for(MIS[[n]])
  if (!length(mm)) { message("  no MIMI column for ", n); next }
  arms <- c(exp_baseline_arms()[c("null_train_mean","domain_index")],
            list(mimi_matched=ridge(mm), mimi_mismatched=ridge(mx),
                 mimi_all=ridge(PCT), index_plus_mimi=iplus(PCT)))
  for (tg in c("level","prev")) {
    cell <- tryCatch(exp_cell(E,"Malawi",ix$outcome[i],tg, extra=G2,
                     extra_domain="Dietary inadequacy (MIMI Malawi)"),
                     error=function(e) NULL)
    if (is.null(cell)) next
    r <- exp_infill(cell, arms, reps=REPS); r$nutrient <- n
    rows[[paste(i,tg)]] <- r
  }
  message("  ", ix$outcome[i], " (", n, ")")
}
SM <- exp_summarise(dplyr::bind_rows(rows)); exp_write(SM,"27_mimi_mw_score")
a <- SM[SM$estimand=="infill" & SM$target=="level",]
cat("\n== Malawi in-fill, level: median Spearman ==\n")
print(aggregate(spearman ~ arm, data=a, FUN=function(z) round(median(z,na.rm=TRUE),3)), row.names=FALSE)
w <- reshape(a[,c("outcome","arm","spearman")], idvar="outcome", timevar="arm", direction="wide")
names(w) <- sub("^spearman[.]","",names(w))
w$gain <- round(w$index_plus_mimi - w$domain_index,3)
w$spec <- round(w$mimi_matched - w$mimi_mismatched,3)
cat("\n== per cell ==\n")
print(w[order(-w$gain), c("outcome","domain_index","index_plus_mimi","gain","mimi_matched","mimi_mismatched","spec")],
      row.names=FALSE, digits=3)
cat(sprintf("\nindex + MIMI vs index: mean %+.4f, better in %d of %d cells\n",
    mean(w$gain,na.rm=TRUE), sum(w$gain>0,na.rm=TRUE), sum(is.finite(w$gain))))
cat(sprintf("specificity (matched - mismatched): mean %+.4f, matched better in %d of %d\n",
    mean(w$spec,na.rm=TRUE), sum(w$spec>0,na.rm=TRUE), sum(is.finite(w$spec))))
