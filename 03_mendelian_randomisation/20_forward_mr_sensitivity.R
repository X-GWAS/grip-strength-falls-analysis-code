## 20_mr_finn_sens.R —— 独立结局 MR 的稳健性诊断
## 内容：MR-PRESSO、留一法、限制工具变量强度、阴性对照（LDL→跌倒）

suppressMessages({ library(TwoSampleMR); library(data.table); library(MRPRESSO) })
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/mr-finn")

FINN_CASE <- 133526; FINN_CTRL <- 366822; FINN_N <- FINN_CASE + FINN_CTRL

fin <- fread("finngen_FALLS_all.tsv", sep = "\t")
setnames(fin, gsub("^#", "", names(fin)))
out_finn <- data.frame(
  SNP = sub(",.*$", "", fin$rsids),
  beta.outcome = fin$beta, se.outcome = fin$sebeta,
  effect_allele.outcome = fin$alt, other_allele.outcome = fin$ref,
  eaf.outcome = fin$af_alt, pval.outcome = fin$pval,
  samplesize.outcome = FINN_N, ncase.outcome = FINN_CASE,
  ncontrol.outcome = FINN_CTRL,
  outcome = "Falls (FinnGen R12)", id.outcome = "finn_R12_FALLS",
  stringsAsFactors = FALSE)

load_exp <- function(file) {
  e <- as.data.frame(fread(file)); e$id.exposure <- file
  o <- out_finn[match(e$SNP, out_finn$SNP), ]; o <- o[!is.na(o$SNP), ]
  d <- harmonise_data(e[e$SNP %in% o$SNP, ], o, action = 2)
  d[d$mr_keep, ]
}

## ---- 1) MR-PRESSO ----
cat("\n########## 1) MR-PRESSO（离群点校正）##########\n")
presso_res <- list()
for (f in c("inst_gripL.csv", "inst_gripR.csv", "inst_lowGrip.csv", "inst_lowGripFNIH.csv")) {
  d <- load_exp(f)
  tag <- sub("inst_|\\.csv", "", f)
  if (nrow(d) < 4) { cat(tag, "SNP 不足，跳过\n"); next }
  pr <- tryCatch(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = as.data.frame(d), NbDistribution = 1000, SignifThreshold = 0.05),
                error = function(e) NULL)
  if (is.null(pr)) { cat(tag, "PRESSO 失败\n"); next }
  m1 <- as.data.table(pr$`Main MR results`)[1]
  b <- suppressWarnings(as.numeric(m1$`Causal Estimate`)); p <- suppressWarnings(as.numeric(m1$`P-value`))
  glob <- suppressWarnings(as.numeric(pr$`MR-PRESSO results`$`Global Test`$Pvalue))
  oi <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  n_out <- if (is.null(oi)) 0L else length(oi)
  cat(sprintf("%-16s n=%3d  全局异质性 p=%.3g | 校正 OR=%.3f (%.3f-%.3f) p=%.3g | 离群点 %d 个\n",
              tag, nrow(d), glob, exp(b), exp(b - 1.96 * m1$Sd), exp(b + 1.96 * m1$Sd), p, n_out))
  presso_res[[tag]] <- data.table(暴露 = tag, n_snp = nrow(d), 全局异质性p = glob,
    校正OR = exp(b), 下限 = exp(b - 1.96 * as.numeric(m1$Sd)), 上限 = exp(b + 1.96 * as.numeric(m1$Sd)),
    p = p, 离群点数 = n_out)
}
fwrite(rbindlist(presso_res), "res/05_MRPRESSO.csv")

## ---- 2) 留一法（握力左手）----
cat("\n########## 2) 留一法：握力左手 → FinnGen 跌倒 ##########\n")
d <- load_exp("inst_gripL.csv")
loo <- mr_leaveoneout(d)
loo <- generate_odds_ratios(loo); setDT(loo)
loo[, 是否越界 := ifelse(or_lci95 > 0.9713 | or_uci95 < 0.9713, "是", "否")]
fwrite(loo[, .(SNP, b, se, p, or, or_lci95, or_uci95, 是否越界)], "res/06_留一法.csv")
cat("留一法结果区间跨越总体估计(0.971)的 SNP 数:",
    sum(loo$是否越界 == "否"), "/", nrow(loo), "\n")
cat("留一法后 OR 范围:", sprintf("%.3f - %.3f", min(loo$or), max(loo$or)), "\n")

## ---- 3) 限制工具变量强度 ----
cat("\n########## 3) 限制工具变量：只用更强的 SNP ##########\n")
e_all <- as.data.frame(fread("inst_gripL.csv"))
o_all <- out_finn[match(e_all$SNP, out_finn$SNP), ]; o_all <- o_all[!is.na(o_all$SNP), ]
res_restrict <- list()
for (thr in c(5e-8, 1e-10, 1e-12, 1e-15)) {
  ee <- e_all[e_all$pval.exposure < thr & e_all$SNP %in% o_all$SNP, ]
  d2 <- harmonise_data(ee, o_all, action = 2); d2 <- d2[d2$mr_keep, ]
  if (nrow(d2) < 3) next
  r <- generate_odds_ratios(mr(d2, method_list = "mr_ivw"))
  cat(sprintf("p<%g   n=%3d  OR=%.3f (%.3f-%.3f) p=%.3g\n", thr, nrow(d2),
              r$or, r$or_lci95, r$or_uci95, r$pval))
  res_restrict[[as.character(thr)]] <- data.table(阈值 = thr, n_snp = nrow(d2),
              OR = r$or, 下限 = r$or_lci95, 上限 = r$or_uci95, p = r$pval)
}
## 按效应量绝对值取前 N 个
for (nn in c(20, 50)) {
  ee <- e_all[order(-abs(e_all$beta.exposure))][1:nn, ]
  ee <- ee[ee$SNP %in% o_all$SNP, ]
  d2 <- harmonise_data(ee, o_all, action = 2); d2 <- d2[d2$mr_keep, ]
  r <- generate_odds_ratios(mr(d2, method_list = "mr_ivw"))
  cat(sprintf("效应量前%d  n=%3d  OR=%.3f (%.3f-%.3f) p=%.3g\n", nn, nrow(d2),
              r$or, r$or_lci95, r$or_uci95, r$pval))
  res_restrict[[paste0("top", nn)]] <- data.table(阈值 = paste0("top", nn), n_snp = nrow(d2),
              OR = r$or, 下限 = r$or_lci95, 上限 = r$or_uci95, p = r$pval)
}
fwrite(rbindlist(res_restrict), "res/07_限制工具变量.csv")

## ---- 4) 阴性对照：LDL 胆固醇 → 跌倒（不应有因果效应）----
cat("\n########## 4) 阴性对照：LDL → FinnGen 跌倒 ##########\n")
d_ldl <- load_exp("inst_LDL.csv")
r <- generate_odds_ratios(mr(d_ldl, method_list = c("mr_ivw", "mr_weighted_median")))
print(as.data.frame(r[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")]), row.names = FALSE)
fwrite(as.data.table(r)[, .(方法 = method, n_snp = nsnp, OR = or, 下限 = or_lci95,
                            上限 = or_uci95, p = pval)], "res/08_阴性对照_LDL.csv")

cat("\n完成。\n")
