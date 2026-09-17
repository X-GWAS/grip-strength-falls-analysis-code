## 60_charls_设计效应.R —— CHARLS 主分析的设计效应（design effect）**
## 口径与 10_analyze_A.R / 12_charls_fallinj.R 完全一致，仅额外报告设计效应与有效样本量。
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹"
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/00_deff_helper.R")

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

fix_inj <- function(fall, inj_raw) {
  fv <- fall; iv <- as.integer(inj_raw >= 1)
  iv[is.na(fv)] <- NA_integer_; iv[fv == 0] <- 0L; iv
}
W[, fall_inj_1518_fix := fix_inj(fall_1518, w4_da024)]
Wc <- W[!is.na(wt3) & wt3 > 0]
cat("暴露波加权样本:", nrow(Wc), " 社区(PSU):", uniqueN(Wc$communityID),
    " 层:", uniqueN(Wc$stratum), " 权重CV:",
    round(sqrt(nrow(Wc)*sum(Wc$wt3^2)/sum(Wc$wt3)^2 - 1), 3), "\n")

## 【2026-09-17 投稿口径修订】原主模型额外含 scale(slope_yr)（2011–2015 握力年变化速率）。
##   该协变量在 ELSA/HRS/SHARE 主模型中并不存在，且属暴露测量之后、结局窗口期间的事后变量，
##   在主模型里调整会造成过调整。经确认后从主模型、受伤跌倒、髋骨骨折三个模型一并删除，
##   四队列自此真正使用同一分析模板。备份见 work/备份-删速率前-20260917/。
fml_cov <- "age3 + factor(edu) + married + rural + smoke + drink + hibp + diab + bmi + adl + cesd"
cc <- function(d, oc) d[!is.na(d[[oc]]) & !is.na(d$grip_2011) & !is.na(d$age3)]

DEFF <- rbindlist(list(
  survey_deff(as.data.frame(cc(copy(Wc), "fall_1518")),
              paste("fall_1518 ~ scale(grip_2011) +", fml_cov),
              "scale(grip_2011)", "CHARLS 跌倒（主模型）", "wt3",
              psu = "communityID", strata = "stratum", 事件列 = "fall_1518"),
  survey_deff(as.data.frame(cc(copy(Wc), "fall_inj_1518_fix")),
              paste("fall_inj_1518_fix ~ scale(grip_2011) +", fml_cov),
              "scale(grip_2011)", "CHARLS 跌倒受伤", "wt3",
              psu = "communityID", strata = "stratum", 事件列 = "fall_inj_1518_fix"),
  survey_deff(as.data.frame(cc(copy(Wc), "fx_1518")),
              paste("fx_1518 ~ scale(grip_2011) +", fml_cov),
              "scale(grip_2011)", "CHARLS 新发髋骨骨折", "wt3",
              psu = "communityID", strata = "stratum", 事件列 = "fx_1518")
), fill = TRUE)

fwrite(DEFF, file.path(out, "res", "12_CHARLS_设计效应.csv"))
cat("\n===== CHARLS 设计效应 =====\n")
print(DEFF[, .(项, n, 事件, OR_完整设计, 下限, 上限, SE_完整设计, SE_仅加权, SE_简单随机,
               DEFF_总体, DEFF_加权, DEFF_聚类, 权重CV, DEFF_权重Kish, 有效样本量, PSU数, 层数)])
## ---------- 导出主模型分析样本（Table 1 用；2026-09-16 新增）----------
t1dir <- file.path(out, "t1"); dir.create(t1dir, showWarnings = FALSE, recursive = TRUE)
mv <- all.vars(as.formula(paste0("~", fml_cov)))
need <- unique(c("ID","wt3","communityID","stratum","age", all.vars(as.formula(paste0("fall_1518 ~ scale(grip_2011) +", fml_cov)))))
dd <- as.data.frame(Wc)
dd <- dd[complete.cases(dd[, need]), ]
dd$grip_z <- dd$grip_2011 / sd(W$grip_2011, na.rm = TRUE)
saveRDS(dd, file.path(t1dir, "样本_CHARLS.rds"))
cat("\n[Table1] CHARLS 模型分析样本 n =", nrow(dd), " 社区 =", uniqueN(dd$communityID), "\n")
print(table(dd$edu, useNA = "ifany"))
