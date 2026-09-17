## 10_reverse_mr.R —— 反向 MR：跌倒 → 握力
## 主分析（无样本重叠）：FinnGen R12 FALLS（芬兰，133,526 例/366,822 对照）
##                      → UK Biobank 握力（英国）
## 对照（有样本重叠）：UKB「去年跌倒」→ UKB 握力
suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
Sys.setenv(OPENGWAS_JWT = trimws(readLines("/tmp/opengwas_jwt", warn = FALSE)[1]))
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")
dir.create("res", showWarnings = FALSE)
FINN_CASE <- 133526; FINN_CTRL <- 366822; FINN_N <- FINN_CASE + FINN_CTRL

## ---------- 1) 暴露：FinnGen R12 FALLS 的 16 个独立工具变量 ----------
iv <- fread("01_工具变量_剪枝后.csv")
fin <- fread("finn_FALLS_gws.tsv", sep = "\t"); setnames(fin, gsub("^#", "", names(fin)))
fin[, rsid := sub(",.*$", "", rsids)]
setnames(fin, c("beta","sebeta","ref","alt","af_alt","pval", "nearest_genes"),
             c("beta.exp","se.exp","ref.exp","alt.exp","af.exp","pval.exp","gene"))
iv[, pval := NULL]
iv <- merge(iv, fin[, .(rsid, gene, beta.exp, se.exp, ref.exp, alt.exp, af.exp, pval.exp)],
            by = "rsid")
setnames(iv, c("beta.exp","se.exp","ref.exp","alt.exp","af.exp","pval.exp"),
             c("beta.exposure","se.exposure","other_allele.exposure","effect_allele.exposure","eaf.exposure","pval.exposure"))
setnames(iv, "gene", "nearest_genes")
iv[, `:=`(id.exposure = "finn_R12_FALLS",
          exposure = "Genetic liability to falls (FinnGen R12)",
          units.exposure = "log-odds")]
iv[, mr_keep.exposure := TRUE]
setnames(iv, "rsid", "SNP")
cat("暴露工具变量:", nrow(iv), "个\n")
print(iv[, .(SNP, nearest_genes, beta.exposure = round(beta.exposure,4),
             se.exposure = round(se.exposure,4), eaf.exposure = round(eaf.exposure,3),
             pval.exposure = signif(pval.exposure,3))], row.names = FALSE)
## F 统计量（工具变量强度）
iv[, F := (beta.exposure/se.exposure)^2]
cat(sprintf("\n工具变量强度：F 中位 %.1f，最小 %.1f，最大 %.1f；%d/%d 个 F>10\n",
            median(iv$F), min(iv$F), max(iv$F), sum(iv$F>10), nrow(iv)))

## ---------- 2) 结局：UKB 握力（连续 + 低握力二分类）与阳性/阴性对照 ----------
OUTC <- list(
  list(id = "ukb-b-10215",        tag = "握力（右手，UKB）",   type = "cont"),
  list(id = "ukb-b-7478",         tag = "握力（左手，UKB）",   type = "cont"),
  list(id = "ebi-a-GCST90007526", tag = "低握力 EWGSOP（UKB）", type = "bin"),
  list(id = "ebi-a-GCST90007528", tag = "低握力 EWGSOP（UKB，第二版）", type = "bin")
)
main <- list(); diag <- list(); snpinfo <- list()
for (o in OUTC) {
  od <- tryCatch(extract_outcome_data(snps = iv$SNP, outcomes = o$id), error = function(e) NULL)
  if (is.null(od) || nrow(od) == 0) { cat("\n[", o$tag, "] 取不到结局数据，跳过\n"); next }
  names(od)[names(od) == "id.outcome"] <- "id.outcome"
  od$outcome <- o$tag
  d <- harmonise_data(as.data.frame(iv), od, action = 2)
  d <- d[d$mr_keep, ]
  nsnp <- nrow(d)
  cat(sprintf("\n===== 跌倒 → %s | 工具变量 %d 个，协调后 %d 个 =====\n", o$tag, nrow(iv), nsnp))
  if (nsnp < 3) { cat("SNP 太少，跳过\n"); next }
  r <- generate_odds_ratios(mr(d, method_list = c("mr_ivw","mr_egger_regression",
                                                  "mr_weighted_median","mr_weighted_mode")))
  setDT(r); r[, `:=`(结局 = o$tag, 结局类型 = o$type, n_snp = nsnp)]
  main[[o$tag]] <- r[, .(结局, 结局类型, 方法 = method, n_snp, b, se, pval, or, or_lci95, or_uci95)]
  het <- mr_heterogeneity(d); ple <- mr_pleiotropy_test(d)
  st <- tryCatch(directionality_test(d), error = function(e) NULL)
  pr <- tryCatch(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                 SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
                 DISTORTIONtest = TRUE, data = as.data.frame(d), NbDistribution = 2000,
                 SignifThreshold = 0.05), error = function(e) NULL)
  pb <- NA_real_; pp <- NA_real_; no <- NA_integer_
  if (!is.null(pr)) {
    m1 <- as.data.table(pr$`Main MR results`)[1]
    pb <- suppressWarnings(as.numeric(m1$`Causal Estimate`)); pp <- suppressWarnings(as.numeric(m1$`P-value`))
    oi <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    no <- if (is.null(oi)) 0L else length(oi)
  }
  diag[[o$tag]] <- data.table(结局 = o$tag, n_snp = nsnp,
    Q_p = het$Q_pval[het$method == "Inverse variance weighted"],
    Egger截距 = ple$egger_intercept, Egger截距P = ple$pval,
    MRPRESSO_beta = pb, MRPRESSO_p = pp, MRPRESSO离群点 = no,
    Steiger方向正确 = if (!is.null(st)) sprintf("%d/%d", sum(st$correct_causal_direction), nrow(st)) else NA_character_)
  snpinfo[[o$tag]] <- as.data.table(d)[, .(结局 = o$tag, SNP, effect_allele.exposure, beta.exposure,
      se.exposure, eaf.exposure, effect_allele.outcome, beta.outcome, se.outcome,
      eaf.outcome, pval.outcome, mr_keep)]
  fwrite(rbindlist(main, fill=TRUE), "res/01_反向MR_主结果.csv")
  fwrite(rbindlist(diag, fill=TRUE), "res/02_反向MR_诊断.csv")
  ivw <- r[method == "Inverse variance weighted"]
  cat(sprintf("   IVW beta=%+.5f (se %.5f)  p=%.4g\n", ivw$b, ivw$se, ivw$pval))
}
fwrite(rbindlist(main, fill=TRUE), "res/01_反向MR_主结果.csv")
fwrite(rbindlist(diag, fill=TRUE), "res/02_反向MR_诊断.csv")
fwrite(rbindlist(snpinfo, fill=TRUE), "res/03_反向MR_逐SNP.csv")
cat("\n\n===== 反向 MR 主结果汇总 =====\n")
print(rbindlist(main, fill=TRUE)[方法 == "Inverse variance weighted",
      .(结局, n_snp, beta = round(b,5), se = round(se,5), p = signif(pval,3))], row.names = FALSE)
cat("\n===== 诊断 =====\n"); print(rbindlist(diag, fill=TRUE), row.names = FALSE)
cat("\n完成。\n")
