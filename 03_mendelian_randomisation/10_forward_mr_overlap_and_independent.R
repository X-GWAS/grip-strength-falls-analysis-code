## 10_mr_finn.R —— 用独立于 UK Biobank 的跌倒 GWAS 重做握力→跌倒 MR
## 目的：原分析的暴露（ukb-b-7478 握力）与结局（ukb-b-2535 去年跌倒）都来自 UKB，
##       样本重叠会低估标准误。这里把结局换成 FinnGen R12 的 FALLS 端点
##       （133,526 例 / 366,822 对照，芬兰人群），与 UKB 无个体重叠。

suppressMessages({ library(TwoSampleMR); library(data.table) })

setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/mr-finn")
dir.create("res", showWarnings = FALSE)

FINN_CASE  <- 133526
FINN_CTRL  <- 366822
FINN_N     <- FINN_CASE + FINN_CTRL

## ---- 1) 读取 FinnGen 结局 ----
fin <- fread("finngen_FALLS_selected.tsv", sep = "\t")
setnames(fin, gsub("^#", "", names(fin)))
cat("FinnGen FALLS 抽取到变异数:", nrow(fin), "\n")

out_finn <- data.frame(
  SNP                  = fin$rsids,
  beta.outcome         = fin$beta,
  se.outcome           = fin$sebeta,
  effect_allele.outcome = fin$alt,
  other_allele.outcome  = fin$ref,
  eaf.outcome          = fin$af_alt,
  pval.outcome         = fin$pval,
  samplesize.outcome   = FINN_N,
  ncase.outcome        = FINN_CASE,
  ncontrol.outcome     = FINN_CTRL,
  outcome              = "Falls (FinnGen R12, 独立于 UKB)",
  id.outcome           = "finn_R12_FALLS",
  stringsAsFactors     = FALSE
)
## FinnGen 的 rsids 字段可能是逗号分隔的多个 rsID，取第一个用于匹配
out_finn$SNP <- sub(",.*$", "", out_finn$SNP)

## ---- 2) 逐个暴露跑 MR ----
exposures <- list(
  list(file = "inst_gripL.csv",   tag = "握力（左手，UKB)",   presso = TRUE),
  list(file = "inst_gripR.csv",   tag = "握力（右手，UKB)",   presso = TRUE),
  list(file = "inst_lowGrip.csv", tag = "低握力（EWGSOP 定义）", presso = TRUE),
  list(file = "inst_BMI.csv",     tag = "BMI（阳性对照）",     presso = FALSE)
)

main <- list(); diag <- list(); snpinfo <- list()

for (ex in exposures) {
  e <- fread(ex$file)
  e <- as.data.frame(e)
  e$exposure <- ex$tag
  e$id.exposure <- e$exposure

  o <- out_finn[match(e$SNP, out_finn$SNP), ]
  o <- o[!is.na(o$SNP), ]
  e2 <- e[e$SNP %in% o$SNP, ]

  d <- harmonise_data(e2, o, action = 2)
  d <- d[d$mr_keep, ]
  nsnp <- nrow(d)
  cat("\n=====", ex$tag, "| 工具变量", nrow(e), "个，匹配并协调后", nsnp, "个 =====\n")
  flush.console()
  if (nsnp < 3) { cat("SNP 太少，跳过\n"); next }

  r <- mr(d, method_list = c("mr_ivw", "mr_egger_regression",
                             "mr_weighted_median", "mr_weighted_mode"))
  r <- generate_odds_ratios(r)
  setDT(r)
  r[, `:=`(暴露 = ex$tag, 结局 = "FinnGen R12 跌倒", n_snp = nsnp)]
  main[[ex$tag]] <- r[, .(暴露, 方法 = method, n_snp, b, se, pval, or, or_lci95, or_uci95)]

  het  <- mr_heterogeneity(d)
  plei <- mr_pleiotropy_test(d)
  do_presso <- isTRUE(ex$presso) && !nzchar(Sys.getenv("SKIP_PRESSO"))
  pr   <- if (do_presso) tryCatch(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                    BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                    SdExposure = "se.exposure", OUTLIERtest = TRUE,
                    DISTORTIONtest = TRUE, data = as.data.frame(d),
                    NbDistribution = 2000, SignifThreshold = 0.05),
                   error = function(err) NULL) else NULL

  presso_b <- NA_real_; presso_p <- NA_real_; n_out <- NA_integer_
  if (!is.null(pr)) {
    m1 <- as.data.table(pr$`Main MR results`)[1]
    presso_b <- suppressWarnings(as.numeric(m1$`Causal Estimate`))
    presso_p <- suppressWarnings(as.numeric(m1$`P-value`))
    oi <- pr$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    n_out <- if (is.null(oi)) 0L else length(oi)
  }

  steiger <- tryCatch(directionality_test(d), error = function(err) NULL)

  diag[[ex$tag]] <- data.table(
    暴露 = ex$tag,
    n_snp = nsnp,
    Q_p = het$Q_pval[het$method == "Inverse variance weighted"],
    Egger_intercept = plei$egger_intercept,
    Egger_intercept_p = plei$pval,
    MRPRESSO_beta = presso_b,
    MRPRESSO_p = presso_p,
    MRPRESSO_离群点 = n_out,
    Steiger方向正确率 = if (!is.null(steiger)) sprintf("%d/%d",
                       sum(steiger$correct_causal_direction), nrow(steiger)) else NA_character_
  )

  snpinfo[[ex$tag]] <- as.data.table(d)[, .(暴露 = ex$tag, SNP, effect_allele.exposure,
        beta.exposure, se.exposure, eaf.exposure, effect_allele.outcome,
        beta.outcome, se.outcome, eaf.outcome, pval.outcome)]
  ## 边跑边落盘，避免中途失败丢失结果
  fwrite(rbindlist(main),   "res/01_FinnGen_MR_主结果.csv")
  fwrite(rbindlist(diag),   "res/02_FinnGen_MR_诊断.csv")
  fwrite(rbindlist(snpinfo), "res/03_FinnGen_MR_逐SNP.csv")
}

main_dt <- rbindlist(main); diag_dt <- rbindlist(diag)
fwrite(main_dt, "res/01_FinnGen_MR_主结果.csv")
fwrite(diag_dt, "res/02_FinnGen_MR_诊断.csv")
fwrite(rbindlist(snpinfo), "res/03_FinnGen_MR_逐SNP.csv")

cat("\n\n===== 主结果（每 +1 单位暴露的跌倒 OR）=====\n")
print(main_dt[方法 == "Inverse variance weighted"], row.names = FALSE)
cat("\n===== 诊断 =====\n")
print(diag_dt, row.names = FALSE)

## ---- 3) 对照：原有分析（暴露结局同源 UKB，存在样本重叠）----
cat("\n\n===== 对照：UKB 握力 → UKB 跌倒（存在样本重叠的旧分析）=====\n")
suppressMessages(library(ieugwasr))
res_old <- list()
for (cb in list(c("ukb-b-7478", "ukb-b-2535", "握力（左手）→ UKB 去年跌倒"),
                c("ukb-b-7478", "ebi-a-GCST90012857", "握力（左手）→ UKB Falling risk (Trajanoska)"))) {
  e <- extract_instruments(cb[1], p1 = 5e-8, clump = TRUE)
  o <- extract_outcome_data(snps = e$SNP, outcomes = cb[2])
  d <- harmonise_data(e, o); d <- d[d$mr_keep, ]
  r <- generate_odds_ratios(mr(d, method_list = "mr_ivw"))
  cat(sprintf("%-45s SNP=%3d  OR=%.3f (%.3f-%.3f)  p=%.3g\n",
              cb[3], nrow(d), r$or, r$or_lci95, r$or_uci95, r$pval))
  res_old[[cb[3]]] <- data.table(分析 = cb[3], n_snp = nrow(d), or = r$or,
                                 lci = r$or_lci95, uci = r$or_uci95, p = r$pval)
}
fwrite(rbindlist(res_old), "res/04_对照_旧分析_UKB同源.csv")

cat("\n完成。结果写入 res/ 目录。\n")
