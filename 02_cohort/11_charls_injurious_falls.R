## 12_charls_fallinj.R —— 补齐 CHARLS「握力 → 跌倒受伤」结局，用于与 ELSA/HRS 三方对照
## 关键坑：da024（跌倒后是否受伤）只对跌倒者提问，非跌倒者全为缺失。
##         若直接建模，样本只剩跌倒者，变成"病例内部比较"，得到的是假结果。
##         正确做法：fall==0 → 受伤记 0；fall==1 → 取 da024。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹"

D <- as.data.table(readRDS(file.path(out, "build_D.rds")))
D[, ID := as.character(ID)]
old <- fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/ANM-骨折与体能-2015随访/cohort_2011_2015.csv",
             select = c("ID","anm","early45","meno_ever"))
old[, ID := as.character(ID)]
D <- merge(D, unique(old, by = "ID"), by = "ID", all.x = TRUE)
W <- D[ragender == 2 & !is.na(age) & age >= 45]

lg <- melt(W, id.vars = "ID", measure.vars = c("grip_2011","r2gripsum","grip_2015"),
           variable.factor = FALSE, value.name = "grip")
lg[, wave := match(variable, c("grip_2011","r2gripsum","grip_2015"))]
lg[, t := c(0, 2, 4)[wave]]
lg <- lg[!is.na(grip)][, if (.N >= 2) .SD, by = ID]
sl <- lg[, .(slope_yr = as.numeric(coef(lm(grip ~ t))[2]), n_wave = .N), by = ID]
W <- merge(W, sl, by = "ID", all.x = TRUE)
sd_slope <- sd(W$slope_yr, na.rm = TRUE); sd_grip <- sd(W$grip_2011, na.rm = TRUE)

## ---------- 修正版跌倒受伤定义 ----------
fix_inj <- function(fall, inj_raw) {
  fv <- fall
  iv <- as.integer(inj_raw >= 1)            # da024 >= 1 记 1（受伤即算）
  iv[is.na(fv)] <- NA_integer_              # 跌倒状态未知 → 结局未知
  iv[fv == 0] <- 0L                         # 关键：未跌倒 = 未受伤
  iv
}
W[, fall_inj_1518_fix := fix_inj(fall_1518, w4_da024)]
W[, fall_inj_1315_fix := fix_inj(fall_1315, w3_da024)]
W[, inj_inc_1115 := {
  a <- fall_inj_1315_fix; b <- fall_inj_1518_fix
  any_yes <- (!is.na(a) & a == 1) | (!is.na(b) & b == 1)
  both <- !is.na(a) & !is.na(b)
  fifelse(both, as.integer(any_yes), fifelse(any_yes, 1L, NA_integer_))
}]
cat("修正后 2015→2018 跌倒受伤：事件", sum(W$fall_inj_1518_fix == 1, na.rm = TRUE),
    "/ 非缺失", sum(!is.na(W$fall_inj_1518_fix)), "\n")
cat("对照：修正前（只含跌倒者）事件", sum(W$fall_inj_1518 == 1, na.rm = TRUE),
    "/ 非缺失", sum(!is.na(W$fall_inj_1518)), "\n")

## ---------- 抽样设计（与 10_analyze_A.R 完全一致） ----------
Wc <- W[!is.na(wt3) & wt3 > 0]
des <- svydesign(ids = ~communityID, strata = ~stratum, weights = ~wt3,
                 data = as.data.frame(Wc), nest = TRUE)
fml_cov <- "age3 + factor(edu) + married + rural + smoke + drink + hibp + diab + bmi + adl + cesd"
ext <- function(m, term, lab, dat, oc) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(NULL)
  b <- co[term, 1]; se <- co[term, 2]
  data.table(项 = lab, OR = exp(b), 下限 = exp(b - 1.96*se), 上限 = exp(b + 1.96*se),
             p = co[term, 4], n = length(m$y), 事件 = sum(m$y == 1))
}

## 主：2011 握力 → 2015→2018 新发跌倒受伤
d <- subset(des, !is.na(fall_inj_1518_fix) & !is.na(slope_yr) & !is.na(grip_2011) & !is.na(age3))
m <- svyglm(as.formula(paste("fall_inj_1518_fix ~ scale(grip_2011) + scale(slope_yr) +", fml_cov)),
            design = d, family = quasibinomial())
co <- summary(m)$coefficients["scale(grip_2011)", ]
cat(sprintf("\n[CHARLS] 握力→跌倒受伤 n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
                        length(m$y), sum(m$y == 1), exp(co[1]),
            exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4]))
R1 <- ext(m, "scale(grip_2011)", paste0("CHARLS：2011 握力（每 SD≈", round(sd_grip,1), "kg）→ 2015→2018 新发跌倒受伤"),
          d, "fall_inj_1518_fix")

## 副：2011→2015 窗口
d2 <- subset(des, !is.na(inj_inc_1115) & !is.na(grip_2011) & !is.na(age3))
m2 <- svyglm(as.formula(paste("inj_inc_1115 ~ scale(grip_2011) +", fml_cov)),
             design = d2, family = quasibinomial())
co2 <- summary(m2)$coefficients["scale(grip_2011)", ]
cat(sprintf("[CHARLS] 握力→跌倒受伤（2011→2015）n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(d2), sum(d2$variables$inj_inc_1115 == 1), exp(co2[1]),
            exp(co2[1]-1.96*co2[2]), exp(co2[1]+1.96*co2[2]), co2[4]))
R2 <- ext(m2, "scale(grip_2011)", "CHARLS：2011 握力 → 2011→2015 新发跌倒受伤", d2, "inj_inc_1115")

## 对照：不做修正（只含跌倒者）会得到什么
d0 <- subset(des, !is.na(fall_inj_1518) & !is.na(slope_yr) & !is.na(grip_2011) & !is.na(age3))
m0 <- svyglm(as.formula(paste("fall_inj_1518 ~ scale(grip_2011) + scale(slope_yr) +", fml_cov)),
             design = d0, family = quasibinomial())
co0 <- summary(m0)$coefficients["scale(grip_2011)", ]
cat(sprintf("[对照·错误口径] 仅在跌倒者内部比较 n=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(d0), exp(co0[1]), exp(co0[1]-1.96*co0[2]), exp(co0[1]+1.96*co0[2]), co0[4]))

fwrite(rbindlist(list(R1, R2), fill = TRUE), file.path(out, "res", "11_CHARLS_跌倒受伤.csv"))
cat("\n已输出 res/11_CHARLS_跌倒受伤.csv\n")
