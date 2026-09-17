## 20_hrs.R —— HRS 平行复现：握力 → 跌倒 / 跌倒受伤 / 髋骨骨折
## 注意：本 harmonized 文件缺人口学变量，性别与出生年从 Cross-Wave Tracker 合并
suppressMessages({
  source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R")
  library(haven); library(data.table); library(survey)
})
options(survey.lonely.psu = "adjust")
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS"

grip_w <- 8:13; fall_w <- 9:14          # 暴露 w8-w13，结局 w9-w14
## w7 不用：该波 mbmi（BMI）仅 515 人非缺失，加进模型后样本塌到 251
wave_year <- function(w) 1990 + 2 * w   # HRS 每两年一波，w1=1992

vars <- c("hhidpn", "raeducl",
          paste0("r", grip_w, "lgrip1"), paste0("r", grip_w, "rgrip1"),
          paste0("r", fall_w, "fall"), paste0("r", fall_w, "fallinj"),
          paste0("r", grip_w, "fall"),   # 暴露波跌倒，用于"仅基线未跌倒"敏感性分析
          paste0("r", grip_w, "hipe"), paste0("r", fall_w, "hipe"),
          paste0("r", grip_w, "mbmi"), paste0("r", grip_w, "smokef"),
          paste0("r", grip_w, "adlfive"),
          paste0("r", grip_w, "iadlfour"), paste0("r", grip_w, "rxhibp"),
          paste0("r", grip_w, "rxdiab"),
          ## 权重取【暴露波】同波的体检子样本个体权重；取结局波权重会与暴露错位，样本近乎清零
          paste0("r", grip_w, "nwtresp"))

cat("读取 HRS harmonized ...\n")
H <- as.data.table(read_dta(DB$file$hrs_h, col_select = any_of(vars)))
cat("行数:", nrow(H), " 列数:", ncol(H), "\n")

TRK <- "/Volumes/TSD302/七大数据库/3. HRS  美国/HRS_美国/Raw_data/Cross-Wave Tracker File/trk2020tr_r.dta"
## Tracker 无合并好的 hhidpn，主键拆成 HHID + PN，按 HRS 规则 hhidpn = HHID*1000 + PN
## SECU = 抽样误差计算单元（PSU），STRATUM = 抽样层；二者 100% 覆盖，是 HRS 官方推荐的设计变量
T <- as.data.table(read_dta(TRK, col_select = any_of(c("HHID","PN","GENDER","BIRTHYR","SECU","STRATUM"))))
T[, hhidpn := as.numeric(HHID) * 1000 + as.numeric(PN)]
T <- T[!is.na(hhidpn) & !duplicated(hhidpn)]
cat("Tracker 行数:", nrow(T), " 变量:", paste(names(T), collapse=", "), "\n")
cat("主键交集:", sum(H$hhidpn %in% T$hhidpn), "/", nrow(H), "\n")
H <- merge(H, T, by = "hhidpn", all.x = TRUE)
cat("合并后女性数:", sum(H$GENDER == 2, na.rm = TRUE), " / 非缺失性别", sum(!is.na(H$GENDER)), "\n")

## ---------- 长表 ----------
lg <- rbindlist(lapply(seq_along(grip_w), function(k) {
  gw <- grip_w[k]; fw <- fall_w[k]
  a <- suppressWarnings(as.numeric(H[[paste0("r", gw, "lgrip1")]]))
  b <- suppressWarnings(as.numeric(H[[paste0("r", gw, "rgrip1")]]))
  a[a <= 0 | a > 100] <- NA; b[b <= 0 | b > 100] <- NA
  grip <- rowMeans(cbind(a, b), na.rm = TRUE); grip[is.nan(grip)] <- NA
  fv <- fifelse(is.na(H[[paste0("r", fw, "fall")]]), NA_integer_, as.integer(H[[paste0("r", fw, "fall")]] == 1))
  iv <- as.integer(H[[paste0("r", fw, "fallinj")]] == 1)
  fi <- fifelse(is.na(fv), NA_integer_, fifelse(fv == 0L, 0L, iv))
  h0 <- as.integer(H[[paste0("r", gw, "hipe")]] == 1)
  h1 <- as.integer(H[[paste0("r", fw, "hipe")]] == 1)
  ## 严格口径：必须基线 hipe==0 且随访 hipe 非缺失；不把"基线缺失但随访阳性"算作新发
  fx <- fifelse(!is.na(h0) & !is.na(h1) & h0 == 0, h1, NA_integer_)
  f0 <- fifelse(is.na(H[[paste0("r", gw, "fall")]]), NA_integer_, as.integer(H[[paste0("r", gw, "fall")]] == 1))
  data.table(id = H$hhidpn, rnd = k, grip_wave = gw, fall_wave = fw, grip = grip,
             fall = fv, fallinj = fi, fxinc = fx,
             fall0 = f0,
             female = as.integer(H$GENDER == 2), birthyr = H$BIRTHYR,
             edu = as.numeric(H$raeducl),
             bmi = as.numeric(H[[paste0("r", gw, "mbmi")]]),
             ## 修正（2026-09-16）：HRS harmonized 的 rXsmokef 是"每日吸烟支数"，不是吸烟与否。
             ## 原写法 ==1 只会命中"每天恰好 1 支"（13,098 人中仅 42 人），等于一个无意义变量。
             ## 改为 >=1（当前每日吸烟），分析样本中约 15%，与文献一致。
             smoke = as.integer(H[[paste0("r", gw, "smokef")]] >= 1),
             adl = as.numeric(H[[paste0("r", gw, "adlfive")]]),
             iadl = as.numeric(H[[paste0("r", gw, "iadlfour")]]),
             hibp = as.integer(H[[paste0("r", gw, "rxhibp")]] == 1),
             diab = as.integer(H[[paste0("r", gw, "rxdiab")]] == 1),
             strat = as.numeric(H$STRATUM),
             psu   = as.numeric(H$SECU),
             wt = as.numeric(H[[paste0("r", gw, "nwtresp")]]))
}))
lg[, age := wave_year(grip_wave) - birthyr]
lg <- lg[female == 1 & !is.na(age) & age >= 50]
cat("\n窗口记录数:", nrow(lg), " 人数:", uniqueN(lg$id), "\n")
print(lg[, .(n = .N, 握力非缺失 = sum(!is.na(grip)), 握力均值 = round(mean(grip, na.rm = TRUE), 2),
             跌倒非缺失 = sum(!is.na(fall)), 跌倒事件 = sum(fall == 1, na.rm = TRUE),
             受伤事件 = sum(fallinj == 1, na.rm = TRUE), 骨折事件 = sum(fxinc == 1, na.rm = TRUE)), by = rnd])

sd_g <- sd(lg[rnd == 1]$grip, na.rm = TRUE)
lg[, grip_z := grip / sd_g]
cat("握力 SD (w7):", round(sd_g, 2), "kg\n")

## 完整抽样设计：Tracker 的 STRATUM（80 层）+ SECU（PSU）+ 暴露波体检子样本权重
## 注意：harmonized 文件里的 rXsecure 是"工作是否有保障"的就业题，不是 PSU，务必不要误用
lg <- lg[!is.na(wt) & wt > 0 & !is.na(strat) & !is.na(psu)]
des <- svydesign(ids = ~psu, strata = ~strat, weights = ~wt, data = as.data.frame(lg), nest = TRUE)
cat("HRS 设计：PSU", length(unique(lg$psu)), "个，层", length(unique(lg$strat)), "层\n")
## HRS harmonized 里 rxhearte / rxstroke / mstat 在 w7 以后完全缺失，drinkcut 在体检子样本里近乎全缺，均不纳入
cov  <- "age + factor(edu) + smoke + hibp + diab + bmi + adl + iadl + factor(rnd)"
cov1 <- "age + factor(edu) + smoke + hibp + diab + bmi + adl + iadl"
res <- list()
for (oc in c("fall","fallinj","fxinc")) {
  d <- subset(des, !is.na(get(oc)) & !is.na(grip) & !is.na(age))
  m <- svyglm(as.formula(paste(oc, "~ grip_z +", cov)), design = d, family = quasibinomial())
  co <- summary(m)$coefficients["grip_z", ]
  nm <- c(fall="跌倒（任意）", fallinj="跌倒受伤", fxinc="新发髋骨骨折")[[oc]]
  cat(sprintf("\n[HRS] %-10s n=%6d 事件=%5d  每 SD OR=%.3f (%.3f-%.3f) p=%.4f\n",
              oc, length(m$y), sum(m$y == 1), exp(co[1]), exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4]))
  res[[oc]] <- data.table(项 = paste0("HRS：握力（每 SD=", round(sd_g,1), "kg）→ ", nm),
                          OR = exp(co[1]), 下限 = exp(co[1]-1.96*co[2]), 上限 = exp(co[1]+1.96*co[2]),
                          p = co[4], n = length(m$y), 事件 = sum(m$y == 1))
}
per <- rbindlist(lapply(seq_along(grip_w), function(k) {
  d <- subset(des, rnd == k & !is.na(fall) & !is.na(grip) & !is.na(age))
  m <- tryCatch(svyglm(as.formula(paste("fall ~ grip_z +", cov1)), design = d, family = quasibinomial()),
                error = function(e) { cat("  窗口", k, "失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients["grip_z", ]
  data.table(窗口 = paste0("w", grip_w[k], "→w", fall_w[k]), n = length(m$y),
                          事件 = sum(m$y == 1), OR = exp(co[1]),
             下限 = exp(co[1]-1.96*co[2]), 上限 = exp(co[1]+1.96*co[2]), p = co[4])
}))
print(per)
## ---------- 阳性对照：结局是否真的能被已知危险因素抓到 ----------
cat("\n--- 阳性对照（pooled fall 模型，同一样本）---\n")
dc <- subset(des, !is.na(fall) & !is.na(grip) & !is.na(age))
mp <- svyglm(fall ~ grip_z + age + factor(edu) + smoke + hibp + diab + bmi + adl + iadl + factor(rnd),
             design = dc, family = quasibinomial())
print(round(summary(mp)$coefficients[c("age","diab","adl","bmi"), ], 4))

## ---------- 敏感性：只保留暴露波未跌倒者（真"新发"跌倒） ----------
cat("\n--- 敏感性分析：排除基线已跌倒者 ---\n")
sen <- svydesign(ids = ~psu, strata = ~strat, weights = ~wt, data = as.data.frame(lg[fall0 == 0]), nest = TRUE)
ds <- subset(sen, !is.na(fall) & !is.na(grip) & !is.na(age))
ms <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = ds, family = quasibinomial())
cs <- summary(ms)$coefficients["grip_z", ]
cat(sprintf("[HRS-敏感性] fall n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(ds), sum(ds$variables$fall == 1), exp(cs[1]), exp(cs[1]-1.96*cs[2]), exp(cs[1]+1.96*cs[2]), cs[4]))
fwrite(data.table(分析 = "HRS-仅基线未跌倒者", n = nrow(ds), 事件 = sum(ds$variables$fall == 1),
                  OR = exp(cs[1]), 下限 = exp(cs[1]-1.96*cs[2]), 上限 = exp(cs[1]+1.96*cs[2]), p = cs[4]),
       file.path(out, "res", "05_HRS_敏感性.csv"))
fwrite(rbindlist(res, fill = TRUE), file.path(out, "res", "03_HRS_主结果.csv"))
fwrite(per, file.path(out, "res", "04_HRS_分窗口.csv"))
cat("\n===== HRS 汇总 =====\n"); print(rbindlist(res, fill = TRUE)); print(per)

## ---------- 设计效应（2026-09-16 新增）----------
source(file.path(out, "00_deff_helper.R"))
cc <- function(d, oc) d[!is.na(d[[oc]]) & !is.na(d$grip) & !is.na(d$age) & !is.na(d$wt) & d$wt > 0, ]
DEFF <- rbindlist(list(
  survey_deff(as.data.frame(cc(copy(lg), "fall")), paste("fall ~ grip_z +", cov),
              "grip_z", "HRS 跌倒（主模型）", "wt", psu = "psu", strata = "strat", 事件列 = "fall"),
  survey_deff(as.data.frame(cc(copy(lg), "fallinj")), paste("fallinj ~ grip_z +", cov),
              "grip_z", "HRS 跌倒受伤", "wt", psu = "psu", strata = "strat", 事件列 = "fallinj"),
  survey_deff(as.data.frame(cc(copy(lg), "fxinc")), paste("fxinc ~ grip_z +", cov),
              "grip_z", "HRS 新发髋骨骨折", "wt", psu = "psu", strata = "strat", 事件列 = "fxinc")
), fill = TRUE)
fwrite(DEFF, file.path(out, "res", "31_HRS_设计效应.csv"))
cat("\n===== HRS 设计效应 =====\n")
print(DEFF[, .(项, n, 事件, OR_完整设计, SE_完整设计, SE_简单随机, DEFF_总体, 有效样本量, 含聚类, PSU数, 层数)])
## ---------- 导出主模型分析样本（Table 1 用；2026-09-16 新增）----------
t1dir <- file.path(out, "t1"); dir.create(t1dir, showWarnings = FALSE, recursive = TRUE)
mv <- all.vars(as.formula(paste0("~", cov)))
need <- unique(c("id","rnd","grip","wt","strat","psu", all.vars(as.formula(paste0("fall ~ grip_z +", cov)))))
dd <- as.data.frame(des$variables)
dd <- dd[complete.cases(dd[, need]), ]
saveRDS(dd, file.path(t1dir, "样本_HRS.rds"))
cat("\n[Table1] HRS 模型分析样本 n =", nrow(dd), " 去重人数 =", uniqueN(dd$id),
    " 平均窗口数 =", round(nrow(dd)/uniqueN(dd$id), 2), "\n")
print(table(dd$edu, useNA = "ifany"))
