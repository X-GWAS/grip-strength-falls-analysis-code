suppressMessages({library(TwoSampleMR); library(data.table); library(ieugwasr)})
Sys.setenv(OPENGWAS_JWT = trimws(readLines("/tmp/opengwas_jwt", warn=FALSE)[1]))
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")
d <- as.data.table(fread("res/03_反向MR_逐SNP.csv"))
d <- d[结局 == "握力（右手，UKB）"]
d[, wald := beta.outcome/beta.exposure]
d[, wald_se := abs(se.outcome/beta.exposure)]
d[, z := wald/wald_se]
d[, p_wald := 2*pnorm(-abs(z))]
e <- as.data.table(fread("res/12_工具变量表型画像.csv"))
eb <- e[表型=="跟骨骨密度 eBMD", .(SNP, eBMD_beta = beta.outcome, eBMD_p = pval.outcome)]
d <- merge(d, eb, by="SNP", all.x=TRUE)
setorder(d, -wald)
cat("\n===== 每个工具变量对握力的 Wald 比（右手握力）=====\n")
print(d[, .(SNP, effect_allele.exposure, beta.exposure=round(beta.exposure,4),
            wald=round(wald,4), wald_se=round(wald_se,4), p_wald=signif(p_wald,3),
            eBMD_beta=round(eBMD_beta,3), eBMD_p=signif(eBMD_p,3))], row.names=FALSE)
fwrite(d, "res/20_逐SNP_Wald比.csv")

cat("\n\n===== 补充：骨密度(eBMD) → 握力 =====\n")
ee <- extract_instruments("ebi-a-GCST90029004", p1=5e-8, clump=TRUE)
cat("eBMD 工具变量:", nrow(ee), "\n")
od <- extract_outcome_data(snps=ee$SNP, outcomes="ukb-b-10215"); od$samplesize.outcome <- 461089
dd <- harmonise_data(ee, od, action=2); dd <- dd[dd$mr_keep,]
r <- generate_odds_ratios(mr(dd, method_list=c("mr_ivw","mr_weighted_median","mr_egger_regression"))); setDT(r)
print(r[, .(method, nsnp, b=round(b,5), se=round(se,5), pval=signif(pval,3))], row.names=FALSE)
het <- mr_heterogeneity(dd); ple <- mr_pleiotropy_test(dd)
cat(sprintf("Q p=%.3g ; Egger 截距 p=%.3f\n", het$Q_pval[het$method=="Inverse variance weighted"], ple$pval))
fwrite(r, "res/21_eBMD到握力.csv")

cat("\n\n===== 补充：简单跌倒致骨折 → 握力 =====\n")
ef <- extract_instruments("ukb-b-15251", p1=5e-8, clump=TRUE)
cat("骨折工具变量:", if (is.null(ef)) 0 else nrow(ef), "\n")
if (!is.null(ef) && nrow(ef)>0) {
  od2 <- extract_outcome_data(snps=ef$SNP, outcomes="ukb-b-10215"); od2$samplesize.outcome <- 461089
  d2 <- harmonise_data(ef, od2, action=2); d2 <- d2[d2$mr_keep,]
  if (nrow(d2)>=3) {
    r2 <- generate_odds_ratios(mr(d2, method_list=c("mr_ivw","mr_weighted_median"))); setDT(r2)
    print(r2[, .(method, nsnp, b=round(b,5), se=round(se,5), pval=signif(pval,3))], row.names=FALSE)
    fwrite(r2, "res/22_骨折到握力.csv")
  } else cat("SNP 不足\n")
}
