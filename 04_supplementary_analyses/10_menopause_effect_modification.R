## 80_anm_effect_modification.R —— CHARLS：绝经年龄（ANM）对「握力→跌倒」关联的效应修饰
##
## 样本与协变量对齐正文主分析：
##   fall_1518 ~ scale(grip_2011) + age3 + factor(edu) + married + rural + smoke + drink +
##                hibp + diab + bmi + adl + cesd
##   设计：ids=~communityID, strata=~stratum, weights=~wt3, nest=TRUE
##
## 两套口径并列，因为 CHARLS 现行脚本比其他三个队列多一项 scale(slope_yr)：
##   A0 不含速率 —— 与 ELSA/HRS/SHARE 主模型同模板，也是正文 Methods 描述的口径
##   A1 含速率   —— CHARLS 现行脚本口径（正文未描述）
##
## ANM 来自 work/ANM-骨折与体能-2015随访/cohort_2011_2015.csv
## （2011 年 da028_1/da028_2，已按台账 C022 合并 zda028_2 携带值，范围限 11–70）。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")

B   <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹"
OUT <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/妇科-ANM效应修饰"
dir.create(file.path(OUT, "res"), showWarnings = FALSE, recursive = TRUE)

## ---------- 1. 组装样本 ----------
D <- as.data.table(readRDS(file.path(B, "build_D.rds")))
D[, ID := as.character(ID)]
anm_link <- fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/ANM-骨折与体能-2015随访/cohort_2011_2015.csv",
                  select = c("ID", "anm", "early45", "meno_ever"))
anm_link[, ID := as.character(ID)]
D <- merge(D, unique(anm_link, by = "ID"), by = "ID", all.x = TRUE)

W <- D[ragender == 2 & !is.na(age) & age >= 45]
lg <- melt(W, id.vars = "ID", measure.vars = c("grip_2011", "r2gripsum", "grip_2015"),
           variable.factor = FALSE, value.name = "grip")
lg[, t := c(0, 2, 4)[match(variable, c("grip_2011", "r2gripsum", "grip_2015"))]]
lg <- lg[!is.na(grip)][, if (.N >= 2) .SD, by = ID]
sl <- lg[, .(slope_yr = as.numeric(coef(lm(grip ~ t))[2])), by = ID]
W <- merge(W, sl, by = "ID", all.x = TRUE)
W <- W[!is.na(wt3) & wt3 > 0]

fml_cov <- "age3 + factor(edu) + married + rural + smoke + drink + hibp + diab + bmi + adl + cesd"
need <- c("fall_1518", "grip_2011", "age3", "edu", "married", "rural", "smoke", "drink",
          "hibp", "diab", "bmi", "adl", "cesd", "wt3", "communityID", "stratum", "anm")
cov_cols <- setdiff(need, "anm")     # data.table 的 ..x 不做内联运算，先算好

## ---------- 2. 先复现正文那条 CHARLS 估计（证明速率项的影响）----------
run_glm <- function(dat, rhs, label) {
  dat <- as.data.frame(dat)
  des <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3, data = dat, nest = TRUE)
  m <- svyglm(as.formula(paste("fall_1518 ~", rhs)), design = des, family = quasibinomial())
  co <- summary(m)$coefficients
  tm <- grep("scale\\(grip_2011\\)", rownames(co), value = TRUE)[1]
  b <- co[tm, 1]; se <- co[tm, 2]
  data.table(模型 = label, n = length(m$y), 事件 = sum(dat$fall_1518 == 1, na.rm = TRUE),
             OR = exp(b), 下限 = exp(b - 1.96 * se), 上限 = exp(b + 1.96 * se), p = co[tm, 4])
}
base_noslope <- W[complete.cases(W[, ..cov_cols])]
base_slope   <- base_noslope[!is.na(slope_yr)]
rep_tab <- rbind(
  run_glm(base_noslope, paste("scale(grip_2011) +", fml_cov), "主口径（不含速率，A0）"),
  run_glm(base_slope,   paste("scale(grip_2011) + scale(slope_yr) +", fml_cov), "含速率（现行脚本，A1）")
)
cat("===== 复现正文 CHARLS 估计 =====\n"); print(rep_tab)
fwrite(rep_tab, file.path(OUT, "res", "01_主口径复现.csv"))

## ---------- 3. ANM 效应修饰：两套口径并列 ----------
A0 <- W[complete.cases(W[, ..need])]
A1 <- A0[!is.na(slope_yr)]

run_anm <- function(A, label) {
  A <- copy(A)
  A[, anm_grp := cut(anm, c(-Inf, 45, 50, 53, Inf),
                     labels = c("<45", "45–49", "50–52", "≥53"), right = FALSE)]
  A[, early := fifelse(anm < 45, "早绝经(<45)", "非早绝经(≥45)")]
  A[, grip_z := scale(grip_2011)[, 1]]

  cat("\n\n##########", label, "| n =", nrow(A), "##########\n")
  print(A[, .(n = .N, 事件 = sum(fall_1518 == 1), 事件率 = round(mean(fall_1518 == 1), 3),
              握力均值 = round(mean(grip_2011), 1), 年龄均值 = round(mean(age), 1)),
          by = anm_grp][order(anm_grp)])
  print(A[, .(n = .N, 事件 = sum(fall_1518 == 1),
              事件率 = round(mean(fall_1518 == 1), 3)), by = early])

  desA <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3,
                    data = as.data.frame(A), nest = TRUE)
  pick <- function(m, pat) {
    co <- summary(m)$coefficients
    tm <- grep(pat, rownames(co), value = TRUE)
    if (!length(tm)) return(NULL)
    data.table(项 = tm, OR = exp(co[tm, 1]), 下限 = exp(co[tm, 1] - 1.96 * co[tm, 2]),
               上限 = exp(co[tm, 1] + 1.96 * co[tm, 2]), p = co[tm, 4])
  }
  m_cont  <- svyglm(as.formula(paste("fall_1518 ~ grip_z * I(anm/5) +", fml_cov)),
                    design = desA, family = quasibinomial())
  m_early <- svyglm(as.formula(paste("fall_1518 ~ grip_z * factor(early) +", fml_cov)),
                    design = desA, family = quasibinomial())
  m_grp   <- svyglm(as.formula(paste("fall_1518 ~ grip_z * factor(anm_grp) +", fml_cov)),
                    design = desA, family = quasibinomial())

  cat("\n-- 交互：ANM 连续（每 5 岁）× 握力 --\n"); print(pick(m_cont, "grip_z"))
  cat("\n-- 交互：早绝经(<45) × 握力 --\n");        print(pick(m_early, "grip_z|early"))
  cat("\n-- 交互：四分组 × 握力 --\n");             print(pick(m_grp, "grip_z"))
  wt <- regTermTest(m_grp, ~grip_z:factor(anm_grp))
  we <- regTermTest(m_early, ~grip_z:factor(early))
  cat("\n四分组交互联合 Wald：P =", format.pval(wt$p, digits = 3), "\n")
  cat("早绝经交互联合 Wald：P =", format.pval(we$p, digits = 3), "\n")

  strata_or <- rbindlist(lapply(split(A, A$anm_grp), function(s) {
    if (nrow(s) < 50 || sum(s$fall_1518 == 1) < 10) return(NULL)
    d <- as.data.frame(s)
    des <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3, data = d, nest = TRUE)
    m <- svyglm(as.formula(paste("fall_1518 ~ scale(grip_2011) +", fml_cov)),
                design = des, family = quasibinomial())
    co <- summary(m)$coefficients; tm <- "scale(grip_2011)"
    data.table(分层 = as.character(s$anm_grp[1]), n = length(m$y), 事件 = sum(d$fall_1518 == 1),
               OR = exp(co[tm, 1]), 下限 = exp(co[tm, 1] - 1.96 * co[tm, 2]),
               上限 = exp(co[tm, 1] + 1.96 * co[tm, 2]), p = co[tm, 4])
  }))
  cat("\n-- 分层 OR（握力每 +1 SD）--\n"); print(strata_or)

  early_or <- rbindlist(lapply(split(A, A$early), function(s) {
    d <- as.data.frame(s)
    des <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3, data = d, nest = TRUE)
    m <- svyglm(as.formula(paste("fall_1518 ~ scale(grip_2011) +", fml_cov)),
                design = des, family = quasibinomial())
    co <- summary(m)$coefficients; tm <- "scale(grip_2011)"
    data.table(分层 = as.character(s$early[1]), n = length(m$y), 事件 = sum(d$fall_1518 == 1),
               OR = exp(co[tm, 1]), 下限 = exp(co[tm, 1] - 1.96 * co[tm, 2]),
               上限 = exp(co[tm, 1] + 1.96 * co[tm, 2]), p = co[tm, 4])
  }))
  cat("\n-- 分层 OR：早绝经 vs 非早绝经 --\n"); print(early_or)
  list(A = A, strata = strata_or, early = early_or,
       inter_p = c(连续每5岁 = pick(m_cont, "grip_z:I")$p,
                   早绝经二分类 = we$p,
                   四分组联合 = wt$p))
}

res0 <- run_anm(A0, "口径 A0：不含速率（与其余三队列同模板）")
res1 <- run_anm(A1, "口径 A1：含速率（CHARLS 现行脚本）")

cat("\n\n===== 两口径交互 P 对比 =====\n")
print(rbind(data.table(口径 = "A0 不含速率", t(res0$inter_p)),
            data.table(口径 = "A1 含速率",   t(res1$inter_p)), fill = TRUE))

fwrite(res0$strata, file.path(OUT, "res", "02_分层OR_四组_A0.csv"))
fwrite(res1$strata, file.path(OUT, "res", "02_分层OR_四组_A1.csv"))
fwrite(res0$early,  file.path(OUT, "res", "03_分层OR_早绝经_A0.csv"))
fwrite(res1$early,  file.path(OUT, "res", "03_分层OR_早绝经_A1.csv"))
fwrite(res0$A[, .(ID, anm, early45, early, anm_grp, age, age3, edu, married, rural, smoke, drink,
                  hibp, diab, bmi, adl, cesd, grip_2011, fall_1518, wt3, communityID, stratum)],
       file.path(OUT, "res", "04_分析样本_含ANM_A0.csv"))
fwrite(res1$A[, .(ID, anm, early45, early, anm_grp, age, age3, edu, married, rural, smoke, drink,
                  hibp, diab, bmi, adl, cesd, grip_2011, slope_yr, fall_1518, wt3, communityID, stratum)],
       file.path(OUT, "res", "04_分析样本_含ANM_A1.csv"))
cat("\n已输出：res/01–04\n")
