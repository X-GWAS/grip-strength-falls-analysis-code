## 80_Table1_基线特征.R —— 四队列女性分析样本基线特征表（Table 1）
## 输入：各队列 t1/样本_*.rds（主模型完整案例样本）
## 口径：基线行为"每人首次进入主模型分析的暴露波"，按各队列抽样设计加权；
##       结局行按全部分析记录（同一人可贡献多个窗口）统计。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")
B <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
T1 <- file.path(B, "外部复现-ELSA与HRS/t1"); dir.create(T1, showWarnings = FALSE, recursive = TRUE)

prep <- function(rds) {
  d <- as.data.table(readRDS(rds))
  if ("ID" %in% names(d) && !"id" %in% names(d)) setnames(d, "ID", "id")
  if (!"rnd" %in% names(d)) d[, rnd := 1]
  d[, minrnd := min(rnd), by = id]
  d[, baseline := rnd == minrnd]
  d
}
des_of <- function(d, ids, strata, wt) {
  x <- as.data.frame(d)
  if (is.null(strata)) svydesign(ids = as.formula(paste0("~", ids)), weights = as.formula(paste0("~", wt)), data = x)
  else svydesign(ids = as.formula(paste0("~", ids)), strata = as.formula(paste0("~", strata)),
                 weights = as.formula(paste0("~", wt)), data = x, nest = TRUE)
}

## 教育三分类（各队列自定义，口径见文末脚注）
edu3 <- function(x, cohort) switch(cohort,
  CHARLS = fifelse(is.na(x), NA_integer_, fifelse(x <= 4, 1L, fifelse(x <= 7, 2L, 3L))),
  ELSA   = fifelse(is.na(x), NA_integer_, fifelse(x <= 1, 1L, fifelse(x == 3, 2L, 3L))),
  HRS    = as.integer(x), SHARE = as.integer(x))

coh <- list(
  CHARLS = list(d = prep(file.path(B, "新暴露-体力活动与握力轨迹/t1/样本_CHARLS.rds")),
                ids = "communityID", strata = "stratum", wt = "wt3",
                grip = "grip_z", grip_kg = "grip_2011",
                fall = "fall_1518", inj = "fall_inj_1518_fix", fx = "fx_1518",
                dep = "cesd", extra = c(married = "married", rural = "rural", drink = "drink")),
  ELSA   = list(d = prep(file.path(T1, "样本_ELSA.rds")), ids = "id", strata = NULL, wt = "wt",
                grip = "grip_z", grip_kg = "grip",
                fall = "fall", inj = "fallinj", fx = "fxinc",
                dep = "cesd", extra = c(mstat = "mstat", drink = "drink", hearte = "hearte", stroke = "stroke")),
  HRS    = list(d = prep(file.path(T1, "样本_HRS.rds")), ids = "psu", strata = "strat", wt = "wt",
                grip = "grip_z", grip_kg = "grip",
                fall = "fall", inj = "fallinj", fx = "fxinc",
                dep = NULL, extra = c(iadl = "iadl")),
  SHARE  = list(d = prep(file.path(T1, "样本_SHARE.rds")), ids = "psu", strata = "strat", wt = "wt",
                grip = "grip_z", grip_kg = "grip",
                fall = "fall", inj = "fallinj", fx = "fxinc",
                dep = "eurod", extra = c(mstat = "mstat", hearte = "hearte", stroke = "stroke"))
)

res <- list()
for (cn in names(coh)) {
  z <- coh[[cn]]; d <- z$d
  d[, edu3 := edu3(get("edu"), cn)]
  db <- d[baseline == TRUE]
  db[, `:=`(edu_low = as.integer(edu3 == 1), edu_mid = as.integer(edu3 == 2),
            edu_high = as.integer(edu3 == 3))]
  des <- des_of(db, z$ids, z$strata, z$wt)
  wm <- function(v) { f <- as.formula(paste0("~", v)); as.numeric(coef(svymean(f, des, na.rm = TRUE)))[1] }
  wp <- function(v) { f <- as.formula(paste0("~", v)); as.numeric(coef(svymean(f, des, na.rm = TRUE)))[1] }
  row <- function(变量, 值) res[[length(res) + 1]] <<- data.table(变量 = 变量, 队列 = cn, 值 = 值)
  row("分析记录数（人×窗口）", nrow(d))
  row("去重人数", uniqueN(d$id))
  row("基线平均随访窗口数", round(nrow(d) / uniqueN(d$id), 2))
  row("年龄，均值 (SD)", sprintf("%.1f (%.1f)", wm("age"), sqrt(as.numeric(svyvar(~age, des, na.rm = TRUE)))))
  row("教育：低（%）", round(100 * wp("edu_low"), 1))
  row("教育：中（%）", round(100 * wp("edu_mid"), 1))
  row("教育：高（%）", round(100 * wp("edu_high"), 1))
  row(if (cn == "HRS") "目前吸烟（≥1 支/日，%）" else "曾经吸烟（%）", round(100 * wp("smoke"), 1))
  row("BMI，均值 (SD)", sprintf("%.1f (%.1f)", wm("bmi"), sqrt(as.numeric(svyvar(~bmi, des, na.rm = TRUE)))))
  row("高血压（%）", round(100 * wp("hibp"), 1))
  row("糖尿病（%）", round(100 * wp("diab"), 1))
  row("ADL 受限项数，均值", round(wm("adl"), 2))
  if (!is.null(z$dep)) row("抑郁症状评分，均值（队列特异量表）", round(wm(z$dep), 2))
  for (nm in names(z$extra)) {
    v <- z$extra[[nm]]; lab <- c(mstat = "有配偶/伴侣（%）", married = "有配偶（%）", rural = "农村居住（%）",
                                 drink = "饮酒（%）", hearte = "心脏病（%）", stroke = "卒中（%）", iadl = "IADL 受限项数，均值")[[nm]]
    val <- if (nm %in% c("iadl")) round(wm(v), 2) else round(100 * wp(v), 1)
    row(lab, val)
  }
  row("握力，均值 (SD) kg", sprintf("%.1f (%.1f)", wm(z$grip_kg),
      sqrt(as.numeric(svyvar(as.formula(paste0("~", z$grip_kg)), des, na.rm = TRUE)))))
  row("分析样本新发跌倒，n", sum(d[[z$fall]] == 1, na.rm = TRUE))
  row("分析样本新发跌倒，%", round(100 * mean(d[[z$fall]] == 1, na.rm = TRUE), 1))
  if (!is.null(d[[z$inj]])) { row("分析样本新发跌倒受伤，n", sum(d[[z$inj]] == 1, na.rm = TRUE)) }
  if (!is.null(d[[z$fx]]))  { row("分析样本新发髋骨骨折，n", sum(d[[z$fx]] == 1, na.rm = TRUE)) }
}
T1L <- rbindlist(res)
W <- dcast(T1L, 变量 ~ 队列, value.var = "值")
ord <- c("分析记录数（人×窗口）","去重人数","基线平均随访窗口数","年龄，均值 (SD)",
         "教育：低（%）","教育：中（%）","教育：高（%）","有配偶/伴侣（%）","有配偶（%）","农村居住（%）",
         "曾经吸烟（%）","目前吸烟（≥1 支/日，%）","饮酒（%）","BMI，均值 (SD)","高血压（%）","糖尿病（%）","心脏病（%）","卒中（%）",
         "ADL 受限项数，均值","IADL 受限项数，均值","抑郁症状评分，均值（队列特异量表）","握力，均值 (SD) kg",
         "分析样本新发跌倒，n","分析样本新发跌倒，%","分析样本新发跌倒受伤，n","分析样本新发髋骨骨折，n")
W <- W[ord[ord %in% W$变量]]
fwrite(W, file.path(T1, "Table1_四队列女性基线特征.csv"))
cat("\n===== Table 1（四队列女性分析样本）=====\n"); print(W)
