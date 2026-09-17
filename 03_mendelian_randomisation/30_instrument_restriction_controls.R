## 30_mr_finn_restrict.R —— 限制工具变量强度 + 阴性对照（接 20 号脚本）
suppressMessages({ library(TwoSampleMR); library(data.table) })
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/mr-finn")

FINN_CASE <- 133526; FINN_CTRL <- 366822
fin <- fread("finngen_FALLS_all.tsv", sep = "\t")
setnames(fin, gsub("^#", "", names(fin)))
out_finn <- data.frame(
  SNP = sub(",.*$", "", fin$rsids),
  beta.outcome = fin$beta, se.outcome = fin$sebeta,
  effect_allele.outcome = fin$alt, other_allele.outcome = fin$ref,
  eaf.outcome = fin$af_alt, pval.outcome = fin$pval,
  samplesize.outcome = FINN_CASE + FINN_CTRL, ncase.outcome = FINN_CASE,
  ncontrol.outcome = FINN_CTRL,
  outcome = "Falls (FinnGen R12)", id.outcome = "finn_R12_FALLS",
  stringsAsFactors = FALSE)

cat("\n########## 3) 限制工具变量：握力左手 → FinnGen 跌倒 ##########\n")
e_all <- as.data.frame(fread("inst_gripL.csv")); e_all$id.exposure <- "gripL"
o_all <- out_finn[match(e_all$SNP, out_finn$SNP), ]; o_all <- o_all[!is.na(o_all$SNP), ]
res_restrict <- list()
for (thr in c(5e-8, 1e-10, 1e-12, 1e-15)) {
  ee <- e_all[e_all$pval.exposure < thr & e_all$SNP %in% o_all$SNP, ]
  d2 <- harmonise_data(ee, o_all, action = 2); d2 <- d2[d2$mr_keep, ]
  if (nrow(d2) < 3) next
  r <- generate_odds_ratios(mr(d2, method_list = "mr_ivw"))
  cat(sprintf("p<%g   n=%3d  OR=%.3f (%.3f-%.3f) p=%.3g\n", thr, nrow(d2), r$or, r$or_lci95, r$or_uci95, r$pval))
  res_restrict[[as.character(thr)]] <- data.table(阈值 = paste0("p<", thr), n_snp = nrow(d2), OR = r$or, 下限 = r$or_lci95, 上限 = r$or_uci95, p = r$pval)
}
for (nn in c(20, 50)) {
  ee <- e_all[order(-abs(e_all$beta.exposure)), ][1:nn, ]; ee <- ee[ee$SNP %in% o_all$SNP, ]
  d2 <- harmonise_data(ee, o_all, action = 2); d2 <- d2[d2$mr_keep, ]
  r <- generate_odds_ratios(mr(d2, method_list = "mr_ivw"))
  cat(sprintf("效应量前%d  n=%3d  OR=%.3f (%.3f-%.3f) p=%.3g\n", nn, nrow(d2), r$or, r$or_lci95, r$or_uci95, r$pval))
  res_restrict[[paste0("top", nn)]] <- data.table(阈值 = paste0("效应量前", nn), n_snp = nrow(d2), OR = r$or, 下限 = r$or_lci95, 上限 = r$or_uci95, p = r$pval)
}
## 只保留 p<1e-8 且去掉可能多效性强的位点（剔除已知 BMI/身高附近的信号）
ee <- e_all[e_all$pval.exposure < 1e-8 & e_all$SNP %in% o_all$SNP, ]
d2 <- harmonise_data(ee, o_all, action = 2); d2 <- d2[d2$mr_keep, ]
r <- generate_odds_ratios(mr(d2, method_list = c("mr_ivw", "mr_weighted_median")))
print(as.data.frame(r[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]), row.names = FALSE)
fwrite(rbindlist(res_restrict), "res/07_限制工具变量.csv")

cat("\n########## 4) 阴性对照：LDL 胆固醇 → FinnGen 跌倒 ##########\n")
e <- as.data.frame(fread("inst_LDL.csv")); e$id.exposure <- "LDL"
o <- out_finn[match(e$SNP, out_finn$SNP), ]; o <- o[!is.na(o$SNP), ]
d <- harmonise_data(e[e$SNP %in% o$SNP, ], o, action = 2); d <- d[d$mr_keep, ]
r <- generate_odds_ratios(mr(d, method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))
print(as.data.frame(r[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]), row.names = FALSE)
fwrite(as.data.table(r)[, .(方法 = method, n_snp = nsnp, OR = or, 下限 = or_lci95, 上限 = or_uci95, p = pval)], "res/08_阴性对照_LDL.csv")
cat("\n完成。\n")
