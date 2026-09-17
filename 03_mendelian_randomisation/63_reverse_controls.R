suppressMessages({library(TwoSampleMR); library(data.table); library(ieugwasr)})
Sys.setenv(OPENGWAS_JWT = trimws(readLines("/tmp/opengwas_jwt", warn=FALSE)[1]))
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")
res <- list()
do <- function(xid, tag, oid="ukb-b-10215"){
  e <- tryCatch(extract_instruments(xid, p1=5e-8, clump=TRUE), error=function(e) NULL)
  if (is.null(e)||nrow(e)==0) { cat(sprintf("%-28s 无工具变量\n", tag)); return(NULL) }
  od <- extract_outcome_data(snps=e$SNP, outcomes=oid)
  if (is.null(od)||nrow(od)==0) { cat(sprintf("%-28s 无结局数据\n", tag)); return(NULL) }
  d <- harmonise_data(e, od, action=2); d <- d[d$mr_keep,]
  if (nrow(d)<3) { cat(sprintf("%-28s SNP 不足\n", tag)); return(NULL) }
  r <- generate_odds_ratios(mr(d, method_list=c("mr_ivw","mr_weighted_median"))); setDT(r)
  ivw <- r[method=="Inverse variance weighted"]
  ## 手算固定效应以便看清是否被随机效应放大
  bf <- sum(d$beta.outcome*d$beta.exposure/d$se.outcome^2)/sum(d$beta.exposure^2/d$se.outcome^2)
  sef <- sqrt(1/sum(d$beta.exposure^2/d$se.outcome^2))
  het <- mr_heterogeneity(d); qp <- het$Q_pval[het$method=="Inverse variance weighted"]
  cat(sprintf("%-28s SNP=%4d  IVW b=%+.5f (se %.5f) p=%.3g | 固定效应 b=%+.5f (se %.5f) | Q p=%.3g\n",
              tag, nrow(d), ivw$b, ivw$se, ivw$pval, bf, sef, qp))
  data.table(对照=tag, id=xid, n_snp=nrow(d), b_IVW=ivw$b, se_IVW=ivw$se, p_IVW=ivw$pval,
             b_固定=bf,  se_固定=sef, Q_p=qp)
}
res[["身高"]]      <- do("ukb-b-10787", "身高 → 握力（阳性对照）")
res[["ALM"]]       <- do("ebi-a-GCST90000027", "四肢肌肉量 → 握力（阳性对照）")
res[["出生体重"]]  <- do("ieu-b-4811", "出生体重 → 握力（阴性对照）")
res[["BMI"]]       <- do("ieu-b-40", "BMI → 握力")
res[["LDL"]]       <- do("ieu-a-300", "LDL → 握力（阴性对照）")
fwrite(rbindlist(res, fill=TRUE), "res/15b_对照_补充.csv")
