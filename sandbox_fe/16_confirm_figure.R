#!/usr/bin/env Rscript
# Confirmation figure: production-SL recipe variants + outcome-specific bundles.
suppressMessages(library(here))
f14<-read.csv(here::here("sandbox_fe","results_14_sl_recipe_fast.csv"))
f14b<-read.csv(here::here("sandbox_fe","results_14b_fullstack.csv"))
bw<-read.csv(here::here("sandbox_fe","results_15_bundles_within.csv"))
bl<-read.csv(here::here("sandbox_fe","results_15_bundles_loco.csv"))

png(here::here("sandbox_fe","FE_confirmation.png"), width=1400, height=1050, res=130)
par(mfrow=c(2,2), mar=c(8,4,3,1))

# A: recipe variants on real SL (fast stack), per cell
f14$cell<-paste(substr(f14$country,1,3),f14$outcome)
m<-tapply(f14$auc, list(f14$cell, f14$variant), I)[, c("current","rank","pca","rankpca")]
barplot(t(m), beside=TRUE, las=2, col=c("orange","darkgreen","steelblue","purple"),
        ylim=c(0.45,0.72), ylab="CV AUC (real SL, fast)", main="A. Recipe variants via production SuperLearner")
legend("topleft", c("current(z)","rank","per-domain PCA","rank+PCA"),
       fill=c("orange","darkgreen","steelblue","purple"), bty="n", cex=0.7)

# B: full-stack confirmation
f14b$cell<-paste(substr(f14b$country,1,3),f14b$outcome)
m2<-tapply(f14b$auc, list(f14b$cell,f14b$variant), I)[,c("current","rank","pca")]
bp<-barplot(t(m2), beside=TRUE, las=2, col=c("orange","darkgreen","steelblue"),
        ylim=c(0.5,0.62), ylab="CV AUC (real SL, FULL stack)", main="B. Full-stack: rank>=current; PCA wins, ~5x faster")
legend("topleft", c("current(z)","rank","PCA"), fill=c("orange","darkgreen","steelblue"), bty="n", cex=0.75)

# C: bundles within-country (matched vs all vs mismatched)
bw$cell<-paste(substr(bw$country,1,3),bw$outcome)
cols<-ifelse(bw$role=="MATCHED","forestgreen",ifelse(bw$role=="all","grey50","firebrick"))
# order within cell: all, MATCHED, mismatched
ord<-order(bw$cell, factor(bw$role,levels=c("all","MATCHED","mismatched")))
bwo<-bw[ord,]
barplot(bwo$auc, names.arg=paste(bwo$cell,substr(bwo$role,1,4)), las=2, col=cols[ord],
        ylim=c(0.45,0.75), ylab="CV AUC", main="C. Bundles within-country\n(green=matched grey=all red=mismatched)", cex.names=0.6)

# D: bundles area-level transfer
bl$lab<-paste(bl$outcome, ifelse(bl$bundle=="all","ALL","bundle"))
cols<-ifelse(bl$bundle=="all","grey50","forestgreen")
barplot(bl$mean_spearman, names.arg=bl$lab, las=2, col=cols, ylim=c(-0.02,0.35),
        ylab="area-level transfer (Spearman)", main="D. Bundles improve/match transfer\nwith 1/3 the features", cex.names=0.7)
abline(h=0,lty=2)
dev.off()
cat("wrote FE_confirmation.png\n")
