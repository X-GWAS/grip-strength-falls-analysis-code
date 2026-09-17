## 20_sens.R —— 反向 MR 敏感性分析
suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
Sys.setenv(OPENGWAS_JWT = trimws(readLines("/tmp/opengwas_jwt", warn = FALSE)[1]))
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/反向MR-跌倒到握力")
FINN_CASE <- 133526; FINN_CTRL <- 366822; FINN_N <- FINN_CASE + FINN_CTRL
dir.create("res", showWarnings = FALSE)

build_exp <- function() {
  iv <- fread("01_工具变量_剪枝后.csv")
  setnames(iv, "rsid", "SNP")
  fin <- fread("finn_FALLS_gws.tsv", sep = "\t"); setnames(fin, gsub("^#", "", names(fin)))
  fin[, SNP := sub(",.*$", "", rsids)]
  fin <- fin[!duplicated(SNP)]
  ed <- merge(iv[, .(SNP)], fin[, .(SNP, nearest_genes, beta, sebeta, ref, alt, af_alt, pval)], by = "SNP")
  setnames(ed, c("beta","sebeta","ref","alt","af_alt","pval"),
               c("beta.exposure","se.exposure","other_allele.exposure",
                 "effect_allele.exposure","eaf.exposure","pval.exposure"))
  ed[, `:=`(id.exposure = "finn_R12_FALLS",
            exposure = "Genetic liability to falls (FinnGen R12)",
            samplesize.exposure = FINN_N, units.exposure = "log-odds",
            mr_keep.exposure = TRUE, nearest_genes = nearest_genes)]
  ed
}
exp_all <- build_exp()
exp_all[, F := (beta.exposure/se.exposure)^2]

## ---------- S1: Steiger 方向性检验（补 samplesize）----------
cat("\n########## S1 Steiger 方向性检验 ##########\n")
steig_rows <- list()
for (o in list(c("ukb-b-10215","握力（右手，UKB）"), c("ukb-b-7478","握力（左手，UKB）"))) {
  od <- extract_outcome_data(snps = exp_all$SNP, outcomes = o[1])
  od$samplesize.outcome <- if (o[1] == "ukb-b-10215") 461089 else 461026
  d <- harmonise_data(as.data.frame(exp_all), od, action = 2); d <- d[d$mr_keep, ]
  st <- tryCatch(directionality_test(d), error = function(e) {cat("ERR:", conditionMessage(e), "\n"); NULL})
  if (!is.null(st)) {
    stx <- as.data.table(st)
    cat(sprintf("%-18s 方向正确 %d/%d , Steiger p=%.3g , SNP 解释暴露方差 %.4f%% / 结局 %.4f%%\n", o[2],
                sum(stx$correct_causal_direction), nrow(stx), stx$steiger_pval[1],
                100*stx$snp_r2.exposure[1], 100*stx$snp_r2.outcome[1]))
    steig_rows[[o[2]]] <- cbind(data.table(结局=o[2]), stx)
  }
}
fwrite(rbindlist(steig_rows, fill=TRUE), "res/10_Steiger.csv")

## ---------- S2: 留一法 ----------
cat("\n########## S2 留一法（右手握力）##########\n")
od <- extract_outcome_data(snps = exp_all$SNP, outcomes = "ukb-b-10215"); od$samplesize.outcome <- 461089
d <- harmonise_data(as.data.frame(exp_all), od, action = 2); d <- d[d$mr_keep, ]
loo <- mr_leaveoneout(d); setDT(loo)
loo <- merge(loo, as.data.table(d)[, .(SNP, nearest_genes)], by = "SNP", all.x = TRUE)
fwrite(loo, "res/11_留一法.csv")
print(loo[order(-abs(b)), .(SNP, nearest_genes, b = round(b,5), se = round(se,5), p = signif(p,3))][1:6], row.names = FALSE)
cat("留一法 beta 范围:", sprintf("%.5f ~ %.5f", min(loo$b), max(loo$b)), "\n")

## ---------- S3: 工具变量到底捕获了什么表型 ----------
cat("\n########## S3 工具变量对骨密度/骨折/人体测量的效应 ##########\n")
probe <- list(
  list("ebi-a-GCST90029004", "跟骨骨密度 eBMD"),
  list("ukb-b-15251",        "简单跌倒致骨折"),
  list("ieu-b-40",           "BMI"),
  list("ukb-b-10787",        "身高"),
  list("ukb-b-7478",         "握力（左手，回看）")
)
probe_rows <- list()
for (p in probe) {
  od <- tryCatch(extract_outcome_data(snps = exp_all$SNP, outcomes = p[[1]]), error = function(e) NULL)
  if (is.null(od) || nrow(od)==0) { cat("  [",p[[2]],"] 无数据\n"); next }
  d <- harmonise_data(as.data.frame(exp_all), od, action = 2); d <- d[d$mr_keep, ]
  dd <- as.data.table(d)[, .(SNP, nearest_genes, pval.outcome, beta.outcome, se.outcome)]
  dd[, z := beta.outcome/se.outcome]
  dd[, 表型 := p[[2]]]
  probe_rows[[p[[2]]]] <- dd
  cat(sprintf("  %-20s 可比较 SNP=%2d ; 名义显著(p<0.05) %2d 个 ; 校正显著(p<0.05/16) %d 个\n",
              p[[2]], nrow(dd), sum(dd$pval.outcome<0.05), sum(dd$pval.outcome<0.05/16)))
}
PR <- rbindlist(probe_rows, fill=TRUE)
fwrite(PR, "res/12_工具变量表型画像.csv")

## 骨密度通路工具变量（对 eBMD 名义显著）
bone_snps <- PR[表型 == "跟骨骨密度 eBMD" & pval.outcome < 0.05, SNP]
cat("\n对 eBMD 名义显著的位点（", length(bone_snps), "个）:", paste(bone_snps, collapse=", "), "\n")
fwrite(data.table(SNP = bone_snps), "res/13_骨密度通路位点.csv")

## ---------- S4: 限制工具变量集 ----------
cat("\n########## S4 工具变量集拆分 ##########\n")
sets <- list(
  list("全部工具变量",  exp_all$SNP),
  list("排除骨密度位点", setdiff(exp_all$SNP, bone_snps)),
  list("仅骨密度位点",  intersect(exp_all$SNP, bone_snps))
)
res4 <- list()
for (s in sets) {
  ee <- exp_all[SNP %in% s[[2]]]
  od <- extract_outcome_data(snps = ee$SNP, outcomes = "ukb-b-10215"); od$samplesize.outcome <- 461089
  d <- harmonise_data(as.data.frame(ee), od, action = 2); d <- d[d$mr_keep, ]
  if (nrow(d) < 3) { cat(sprintf("  %-16s SNP=%d 太少，跳过\n", s[[1]], nrow(d))); next }
  r <- generate_odds_ratios(mr(d, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median")))
  setDT(r); r[, `:=`(工具变量集 = s[[1]], n_snp = nrow(d))]
  res4[[s[[1]]]] <- r[, .(工具变量集, 方法=method, n_snp, b, se, pval)]
  ivw <- r[method=="Inverse variance weighted"]
  het <- mr_heterogeneity(d)
  cat(sprintf("  %-16s SNP=%2d  IVW beta=%+.5f (se %.5f) p=%.3g ; Q p=%.3g\n",
              s[[1]], nrow(d), ivw$b, ivw$se, ivw$pval, het$Q_pval[het$method=="Inverse variance weighted"]))
}
fwrite(rbindlist(res4, fill=TRUE), "res/14_工具变量集拆分.csv")

## ---------- S5: 阳性/阴性对照 ----------
cat("\n########## S5 阳性(BMI)与阴性(LDL)对照：→握力 ##########\n")
ctl <- list(list("ieu-b-40","BMI（阳性对照）"), list("ieu-a-300","LDL（阴性对照）"))
res5 <- list()
for (c1 in ctl) {
  e <- tryCatch(extract_instruments(c1[[1]], p1=5e-8, clump=TRUE),
                error=function(e) {cat("  [",c1[[2]],"] 取工具变量失败:", conditionMessage(e), "\n"); NULL})
  if (is.null(e)||nrow(e)==0) { cat("  [",c1[[2]],"] 无工具变量\n"); next }
  od <- extract_outcome_data(snps = e$SNP, outcomes = "ukb-b-10215"); od$samplesize.outcome <- 461089
  d <- harmonise_data(e, od, action=2); d <- d[d$mr_keep,]
  r <- generate_odds_ratios(mr(d, method_list="mr_ivw")); setDT(r)
  cat(sprintf("  %-18s SNP=%3d  beta=%+.5f (se %.5f) p=%.3g\n", c1[[2]], nrow(d), r$b, r$se, r$pval))
  res5[[c1[[2]]]] <- data.table(对照=c1[[2]], n_snp=nrow(d), b=r$b, se=r$se, pval=r$pval)
}
fwrite(rbindlist(res5, fill=TRUE), "res/15_对照.csv")

## ---------- S6: 有样本重叠的对照臂（UKB 跌倒 → UKB 握力）----------
cat("\n########## S6 对照臂：UKB 去年跌倒 → UKB 握力（样本重叠）##########\n")
res6 <- list()
for (xid in c("ukb-b-2535","ebi-a-GCST90012857")) {
  e <- tryCatch(extract_instruments(xid, p1=5e-8, clump=TRUE), error=function(e) NULL)
  if (is.null(e)||nrow(e)==0) { cat("  [",xid,"] 无工具变量\n"); next }
  od <- extract_outcome_data(snps = e$SNP, outcomes = "ukb-b-10215"); od$samplesize.outcome <- 461089
  d <- harmonise_data(e, od, action=2); d <- d[d$mr_keep,]
  if (nrow(d)<3) { cat("  [",xid,"] SNP 太少\n"); next }
  r <- generate_odds_ratios(mr(d, method_list=c("mr_ivw","mr_weighted_median")))
  setDT(r); r[, `:=`(暴露=xid, n_snp=nrow(d))]
  res6[[xid]] <- r[, .(暴露, 方法=method, n_snp, b, se, pval)]
  ivw <- r[method=="Inverse variance weighted"]
  cat(sprintf("  %-22s SNP=%3d  IVW beta=%+.5f (se %.5f) p=%.3g\n", xid, nrow(d), ivw$b, ivw$se, ivw$pval))
}
fwrite(rbindlist(res6, fill=TRUE), "res/16_重叠对照臂.csv")

## ---------- S7: 可检出效应量与"若为反向因果应看到多大效应" ----------
cat("\n########## S7 功效与实际可检出效应 ##########\n")
ivwR <- fread("res/01_反向MR_主结果.csv")[结局=="握力（右手，UKB）" & 方法=="Inverse variance weighted"]
se_ivw <- ivwR$se
mde80 <- 2.802*se_ivw; mde90 <- 3.242*se_ivw
GRIP_SD <- 6.2   # ELSA 女性握力 SD 6.2 kg（本课题四队列的中位量级）
## 观察性关联 OR 0.863 / 每 +1 SD 握力；FinnGen FALLS 病例比例
p <- FINN_CASE/FINN_N
b_logit <- log(0.863)
b_lin <- b_logit * p * (1-p)                 # 跌倒概率对握力的斜率（每 SD 握力）
b_rev_prob <- b_lin / (p*(1-p))              # 握力(SD) 对跌倒(0/1) 的斜率
b_rev_logor <- b_rev_prob * p * (1-p)        # 换算成"每 1 单位 log-OR"
cat(sprintf("IVW SE = %.4f SD\n", se_ivw))
cat(sprintf("80%% 功效可检出 %.4f SD（≈ %.2f kg）；90%% 功效 %.4f SD（≈ %.2f kg）\n",
            mde80, mde80*GRIP_SD, mde90, mde90*GRIP_SD))
cat(sprintf("若观察性关联完全由'跌倒→握力'解释，反向效应应约为 %.4f SD（≈ %.2f kg）\n",
            abs(b_rev_logor), abs(b_rev_logor)*GRIP_SD))
cat(sprintf("→ 本反向 MR 的最小可检出效应是上述值的 %.1f 倍：功效不足\n",
            mde80/abs(b_rev_logor)))
fwrite(data.table(项=c("IVW_SE_SD","80%功效_MDE_SD","80%功效_MDE_kg","90%功效_MDE_SD",
                       "观察性蕴含的反向效应_SD","观察性蕴含的反向效应_kg","MDE与蕴含效应之比"),
                  值=c(se_ivw, mde80, mde80*GRIP_SD, mde90, abs(b_rev_logor),
                       abs(b_rev_logor)*GRIP_SD, mde80/abs(b_rev_logor))),
       "res/17_功效换算.csv")
cat("\n完成。\n")
