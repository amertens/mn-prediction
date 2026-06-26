#!/usr/bin/env Rscript
# Summary figure for the feature-engineering investigation.
suppressMessages(library(here))
png(here::here("sandbox_fe","FE_summary.png"), width=1400, height=1050, res=130)
par(mfrow=c(2,2), mar=c(7,4,3,1))

# Panel A: effective n (individuals vs distinct area predictor-rows)
cn<-c("SierraLeone","Gambia","Ghana","Malawi"); ind<-c(486,1020,1165,1102); area<-c(14,30,72,87)
bp<-barplot(rbind(ind,area), beside=TRUE, names.arg=cn, col=c("grey70","firebrick"),
            las=2, ylab="count", main="A. Apparent vs effective sample size\n(child_vitA)")
legend("topleft", c("individuals (apparent n)","distinct area predictor-rows (effective n)"),
       fill=c("grey70","firebrick"), bty="n", cex=0.8)

# Panel B: within-country PCA(~27 feat) vs full(~1000+); ranger
r9<-read.csv(here::here("sandbox_fe","results_09_reduction.csv"))
ag<-aggregate(auc~strategy,r9[r9$strategy %in% c("full","pca_dom3","clusterres"),],function(x)mean(x,na.rm=TRUE))
ag<-ag[match(c("full","pca_dom3","clusterres"),ag$strategy),]
barplot(ag$auc, names.arg=c("full\n(~1000 feat)","per-domain PCA\n(~27 feat)","cluster-res\nonly"),
        col="steelblue", ylim=c(0.5,0.6), ylab="within-country CV AUC", las=1,
        main="B. ~27 components match ~1000 features\n(mean over outcomes/countries)")
abline(h=0.5,lty=2)

# Panel C: transform rank vs zscore (lasso, within)
r4<-read.csv(here::here("sandbox_fe","results_04_within.csv"))
sub<-r4[r4$model=="glmnet_lasso" & r4$variant %in% c("zscore+corr","rank+corr"),]
ag<-aggregate(auc~variant+outcome,sub,function(x)mean(x,na.rm=TRUE))
m<-tapply(ag$auc, list(ag$outcome,ag$variant), I)
barplot(t(m[,c("zscore+corr","rank+corr")]), beside=TRUE, col=c("orange","darkgreen"),
        ylim=c(0.5,0.58), las=2, ylab="within-country CV AUC",
        main="C. Rank-normalize > z-score (linear model)")
legend("topright", c("z-score (current)","rank/quantile"), fill=c("orange","darkgreen"), bty="n", cex=0.8)

# Panel D: domain-specific transferable signal (area-level LOCO)
r12<-read.csv(here::here("sandbox_fe","results_12_domain.csv"))
r12<-r12[r12$domain %in% c("GEE","IHME","MAP","FSEC","ALL"),]
m<-tapply(r12$mean_spearman, list(r12$domain, r12$outcome), I)
m<-m[c("GEE","IHME","MAP","FSEC","ALL"),]
barplot(m, beside=TRUE, col=c("forestgreen","firebrick","purple","tan","grey40"),
        las=2, ylab="area-level transfer (Spearman)",
        main="D. Domain-outcome match beats pooling\n(LOCO)")
legend("topright", rownames(m), fill=c("forestgreen","firebrick","purple","tan","grey40"), bty="n", cex=0.75)
abline(h=0,lty=2)
dev.off()
cat("wrote FE_summary.png\n")
