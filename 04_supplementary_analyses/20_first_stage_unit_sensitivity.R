## 02_补齐两条敏感性.R —— 重算两条正文里用到、但原脚本未留存的数据
##   ① 省×城乡（province × urban_nbs）作为分层单位时，CHARLS 主模型标准误的变化
##   ② CHARLS 中 BMI 薄端（<18.5 kg/m²）人群里 BMI 与跌倒的关联
## 口径与 60_charls_设计效应.R 的主模型一致（已删速率项），样本为删速率后的 t1/样本_CHARLS.rds。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")

B   <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹"
OUT <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/投稿-删速率后补充"
dir.create(file.path(OUT, "res"), showWarnings = FALSE, recursive = TRUE)

d <- as.data.table(readRDS(file.path(B, "t1", "样本_CHARLS.rds")))
fml_cov <- "age3 + factor(edu) + married + rural + smoke + drink + hibp + diab + bmi + adl + cesd"
need <- c("fall_1518", "grip_2011", "age3", "edu", "married", "rural", "smoke", "drink",
          "hibp", "diab", "bmi", "adl", "cesd", "wt3", "communityID", "stratum", "province", "urban_nbs")
d <- d[complete.cases(d[, ..need])]
d[, strata_prov := paste(province, urban_nbs, sep = "_")]

## 注：样本里的 stratum 本身就是「省_城乡」（52 层），所以主分析的分层已经是省×城乡；
##     正文所说「换更粗的第一阶段单位」，指的是把 PSU 从社区换成省×城乡。
fit2 <- function(dat, ids_var, strata_var, label) {
  dat <- as.data.frame(dat)
  fml <- if (is.null(strata_var)) NULL else as.formula(paste0("~", strata_var))
  des <- svydesign(ids = as.formula(paste0("~", ids_var)), strata = fml,
                   weights = ~wt3, data = dat, nest = TRUE)
  m <- svyglm(as.formula(paste("fall_1518 ~ scale(grip_2011) +", fml_cov)),
              design = des, family = quasibinomial())
  co <- summary(m)$coefficients["scale(grip_2011)", ]
  data.table(设定 = label, n = length(m$y), OR = exp(co[1]), SE_logOR = co[2],
             下限 = exp(co[1] - 1.96 * co[2]), 上限 = exp(co[1] + 1.96 * co[2]), p = co[4])
}
cat("===== ① 第一阶段单位对比 =====\n")
cmp <- rbind(
  fit2(d, "communityID", "stratum",       "主分析：PSU=社区，分层=省×城乡"),
  fit2(d, "strata_prov", NULL,            "敏感性：PSU=省×城乡（无分层）"),
  fit2(d, "strata_prov", "province",      "敏感性：PSU=省×城乡，分层=省")
)
print(cmp)
se0 <- cmp[1]$SE_logOR
for (i in 2:3) cat(sprintf("%s：标准误变化 %+.1f%%（%.4f → %.4f），点估计 %.3f → %.3f\n",
                           cmp[i]$设定, 100*(cmp[i]$SE_logOR/se0 - 1), se0, cmp[i]$SE_logOR, cmp[1]$OR, cmp[i]$OR))
fwrite(cmp, file.path(OUT, "res", "01_第一阶段单位敏感性.csv"))

cat("\n===== ② BMI 薄端（<18.5）=====\n")
thin <- d[bmi < 18.5]
cat("薄端样本 n =", nrow(thin), " 事件 =", sum(thin$fall_1518 == 1), "\n")
if (nrow(thin) > 60 && sum(thin$fall_1518 == 1) > 10) {
  thin2 <- as.data.frame(thin)
  des2 <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3, data = thin2, nest = TRUE)
  ## 薄端内以 BMI 为暴露（每 +5 kg/m²），协变量不含 bmi
  m2 <- svyglm(fall_1518 ~ I(bmi / 5) + age3 + factor(edu) + married + rural + smoke + drink +
                 hibp + diab + adl + cesd, design = des2, family = quasibinomial())
  co2 <- summary(m2)$coefficients["I(bmi/5)", ]
  res2 <- data.table(分析 = "CHARLS BMI<18.5 内 BMI(每+5) → 跌倒", n = length(m2$y),
                     事件 = sum(thin2$fall_1518 == 1), OR = exp(co2[1]),
                     下限 = exp(co2[1] - 1.96 * co2[2]), 上限 = exp(co2[1] + 1.96 * co2[2]), p = co2[4])
  print(res2)
  fwrite(res2, file.path(OUT, "res", "02_BMI薄端.csv"))
} else {
  cat("薄端样本太小，无法建模\n")
}
