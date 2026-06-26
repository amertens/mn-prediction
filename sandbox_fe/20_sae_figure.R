#!/usr/bin/env Rscript
suppressMessages(library(here))
S <- read.csv(here::here("sandbox_fe","results_19_sae_summary.csv"))
res <- readRDS(here::here("sandbox_fe","results_19_sae.rds"))
S$cell <- paste(substr(S$country,1,3), sub("child_","c.",sub("women_","w.",S$outcome)))

png(here::here("sandbox_fe","FE_sae.png"), width=1400, height=1050, res=130)
par(mfrow=c(2,2), mar=c(7,4,3,1))

# A: optimism gap
m <- rbind(insample=S$cor_insample, LOAO_CV=S$cor_loaocv)
barplot(m, beside=TRUE, names.arg=S$cell, col=c("grey70","firebrick"), las=2,
        ylim=c(0,0.8), ylab="cor(area prev, prediction)",
        main="A. Optimism gap: in-sample fit overstates\nout-of-sample area skill")
legend("topright", c("in-sample (what fit reports)","leave-area-out CV (honest)"),
       fill=c("grey70","firebrick"), bty="n", cex=0.75)

# B: CV skill vs grand mean
skill <- 100*(1 - S$rmse_synth/S$rmse_mean)
cols <- ifelse(skill>0,"forestgreen","firebrick")
bp<-barplot(skill, names.arg=S$cell, col=cols, las=2, ylim=c(-20,25),
        ylab="% RMSE reduction vs grand mean", main="B. Out-of-sample proxy skill\n(>0 beats predicting the mean)")
abline(h=0); text(bp, skill+ifelse(skill>0,2,-2), sprintf("%+.0f%%",skill), cex=0.8)

# C: shrinkage vs #clusters (Malawi child_vitA)
d <- res[["Malawi_child_vitA"]]
plot(jitter(d$nclust), d$gamma, pch=19, col="steelblue", xlab="# survey clusters in Admin2",
     ylab="FH shrinkage weight gamma (1=trust direct)", main="C. SAE borrows strength:\nfew-cluster areas pulled to the model",
     ylim=c(0,1), xlim=c(0.5, max(d$nclust)+0.5))
# mean gamma by #clusters (robust to sparse x); connect with a line
agg <- aggregate(gamma~nclust, d, mean)
lines(agg$nclust, agg$gamma, col="red", lwd=2, type="b", pch=17)

# D: uncertainty half-width direct vs FH
m2 <- rbind(direct=S$hw_direct, FH=S$hw_fh)
barplot(m2, beside=TRUE, names.arg=S$cell, col=c("orange","forestgreen"), las=2,
        ylim=c(0,0.25), ylab="median 95% CI half-width",
        main="D. Honest area uncertainty:\ndirect too noisy; FH tightens via borrowing")
legend("topright", c("direct estimate","FH-EBLUP"), fill=c("orange","forestgreen"), bty="n", cex=0.8)
dev.off()
cat("wrote FE_sae.png\n")
