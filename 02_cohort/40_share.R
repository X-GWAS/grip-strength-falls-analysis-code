## 50_share.R —— SHARE（欧洲）平行复现：握力 → 跌倒 / 新发髋骨骨折
## 设计：暴露在前（握力），结局在后（下一可用波次的跌倒）；女性 ≥50 岁
## 注意：SHARE 的跌倒题是"近 6 个月内是否因跌倒就医/受困扰"（rXfall_s），
##       回忆窗口比 CHARLS/ELSA/HRS 的 1–2 年窄，事件率更低，属已知口径差异。
##       SHARE 没有跌倒受伤变量（全库检索 fallinj 为 0 个），故该结局缺席。
suppressMessages({
  source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R")
  library(haven); library(data.table); library(survey)
})
options(survey.lonely.psu = "adjust")
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS"

grip_w <- c(1, 2, 4, 5, 6)      # 暴露波（w7 因 eurod/smokev/权重覆盖骤降而不用）
fall_w <- c(2, 4, 5, 6, 7)      # 对应结局波（严格在后）
vars <- c("mergeid", "ragender", "raeducl", "country",
          ## 抽样设计：rapsu1/rapsu2 为初级抽样单元，raestrat1/raestrat2 为分层变量（100% 覆盖，空字符串表示该国未使用该套编码）
          "rapsu1", "rapsu2", "raestrat1", "raestrat2",
          paste0("r", c(1,2,4,5,6,7,8), "lgrip1"), paste0("r", c(1,2,4,5,6,7,8), "rgrip1"),
          paste0("r", c(1,2,4,5,6,7,8), "hipe"),
          paste0("r", c(2,4,5,6,7,8), "hip"),      # 自上次访问以来髋骨骨折（更干净的新发定义）
          paste0("r", unique(c(grip_w, fall_w)), "fall_s"),
          paste0("r", grip_w, "agey"), paste0("r", grip_w, "bmi"),
          paste0("r", grip_w, "adlwa"), paste0("r", grip_w, "eurod"),
          paste0("r", grip_w, "smokev"), paste0("r", grip_w, "mstat"),
          paste0("r", grip_w, "hibpe"), paste0("r", grip_w, "diabe"),
          paste0("r", grip_w, "hearte"), paste0("r", grip_w, "stroke"),
          paste0("r", grip_w, "wtresp"))

cat("读取 SHARE harmonized ...\n")
D <- as.data.table(read_dta(DB$file$share_h, col_select = any_of(vars)))
cat("行数:", nrow(D), " 列数:", ncol(D), "\n")

for (w in unique(c(grip_w, fall_w))) {
  a <- paste0("r", w, "lgrip1"); b <- paste0("r", w, "rgrip1")
  av <- if (a %in% names(D)) as.numeric(D[[a]]) else rep(NA_real_, nrow(D))
  bv <- if (b %in% names(D)) as.numeric(D[[b]]) else rep(NA_real_, nrow(D))
  av[av <= 0 | av > 100] <- NA; bv[bv <= 0 | bv > 100] <- NA
  g <- rowMeans(cbind(av, bv), na.rm = TRUE); g[is.nan(g)] <- NA
  D[[paste0("grip", w)]] <- g
}
cat("rXhip（自上次访问以来）非缺失:", sapply(c(2,4,5,6,7,8), function(w) {
  v <- paste0("r", w, "hip"); sum(!is.na(as.numeric(D[[v]])))
}), "\n")
for (w in unique(c(grip_w, fall_w))) {
  v <- paste0("r", w, "fall_s")
  D[[paste0("fall", w)]] <- if (v %in% names(D)) fifelse(is.na(D[[v]]), NA_integer_, as.integer(D[[v]] == 1)) else NA_integer_
}
## 新发髋骨骨折：暴露波 0 → 结局波 1
for (k in seq_along(grip_w)) {
  a <- paste0("r", grip_w[k], "hipe"); b <- paste0("r", fall_w[k], "hipe")
  av <- if (a %in% names(D)) as.integer(D[[a]] == 1) else NA_integer_
  bv <- if (b %in% names(D)) as.integer(D[[b]] == 1) else NA_integer_
  ## 严格口径：必须基线 hipe==0 且随访 hipe 非缺失
  D[[paste0("fxinc", fall_w[k])]] <- fifelse(!is.na(av) & !is.na(bv) & av == 0, bv, NA_integer_)
}

lg <- rbindlist(lapply(seq_along(grip_w), function(k) {
  gw <- grip_w[k]; fw <- fall_w[k]
  data.table(id = D$mergeid, rnd = k, grip_wave = gw, fall_wave = fw,
             grip = D[[paste0("grip", gw)]],
             fall = D[[paste0("fall", fw)]],
             fall0 = D[[paste0("fall", gw)]],     # 暴露波跌倒，用于基线未跌倒者敏感性分析
             fxinc = D[[paste0("fxinc", fw)]],
             ## 直接新发定义：基线无髋骨骨折 且 本波报告"自上次访问以来"发生
             fx_sw = fifelse(is.na(as.integer(D[[paste0("r", grip_w[k], "hipe")]] == 1)) |
                             as.integer(D[[paste0("r", grip_w[k], "hipe")]] == 1) == 1, NA_integer_,
                             as.integer(D[[paste0("r", fw, "hip")]] == 1)),
             age = as.numeric(D[[paste0("r", gw, "agey")]]),
             female = as.integer(D$ragender == 2),
             edu = as.numeric(D$raeducl),
             bmi = as.numeric(D[[paste0("r", gw, "bmi")]]),
             adl = as.numeric(D[[paste0("r", gw, "adlwa")]]),
             eurod = as.numeric(D[[paste0("r", gw, "eurod")]]),
             smoke = as.integer(D[[paste0("r", gw, "smokev")]] == 1),
             mstat = as.integer(D[[paste0("r", gw, "mstat")]] %in% c(1, 2)),
             hibp = as.integer(D[[paste0("r", gw, "hibpe")]] == 1),
             diab = as.integer(D[[paste0("r", gw, "diabe")]] == 1),
             hearte = as.integer(D[[paste0("r", gw, "hearte")]] == 1),
             stroke = as.integer(D[[paste0("r", gw, "stroke")]] == 1),
             country = as.numeric(D$country),
             wt = as.numeric(D[[paste0("r", gw, "wtresp")]]))   # 权重取暴露波
}))
lg <- lg[female == 1 & !is.na(age) & age >= 50]

## ---------- 抽样设计变量构造（2026-09-16 新增）----------
## SHARE 各国使用 1 套或 2 套设计编码，未使用者存放空字符串。规则：
##   优先取 raestrat1/rapsu1，为空则取 raestrat2/rapsu2；
##   两者皆缺的记录（约 15% 来自 8 个国家）单独成层并各给唯一 PSU（相当于不加聚类信息，避免误设相关结构）。
okc <- function(x) !(is.na(x) | x %in% c("", "NA", "."))
v1 <- as.character(D$raestrat1)[match(lg$id, D$mergeid)]
v2 <- as.character(D$raestrat2)[match(lg$id, D$mergeid)]
p1 <- as.character(D$rapsu1)[match(lg$id, D$mergeid)]
p2 <- as.character(D$rapsu2)[match(lg$id, D$mergeid)]
lg[, strat := fifelse(okc(v1), paste0(country, "_", v1),
               fifelse(okc(v2), paste0(country, "_", v2), paste0(country, "_nostrat")))]
lg[, psu := fifelse(okc(p1), paste0(country, "_", p1),
             fifelse(okc(p2), paste0(country, "_", p2), paste0("uid_", .I)))]
cat("SHARE 设计：可用真实 PSU 的记录", sum(!grepl("^uid_", lg$psu)), "/", nrow(lg),
    "；层数", uniqueN(lg$strat), "；PSU 数", uniqueN(lg$psu), "\n")
cat("\n窗口记录数:", nrow(lg), " 人数:", uniqueN(lg$id), "\n")
print(lg[, .(n = .N, 握力非缺失 = sum(!is.na(grip)), 握力均值 = round(mean(grip, na.rm = TRUE), 2),
             跌倒非缺失 = sum(!is.na(fall)), 跌倒事件 = sum(fall == 1, na.rm = TRUE),
             骨折事件 = sum(fxinc == 1, na.rm = TRUE)), by = rnd])

sd_g <- sd(lg[rnd == 1]$grip, na.rm = TRUE)
lg[, grip_z := grip / sd_g]
cat("握力 SD (w1):", round(sd_g, 2), "kg\n")
cat("暴露波跌倒状态非缺失:", sum(!is.na(lg$fall0)), "\n")

lg <- lg[!is.na(wt) & wt > 0]
des <- svydesign(ids = ~psu, strata = ~strat, weights = ~wt, data = as.data.frame(lg), nest = TRUE)
cov  <- "age + factor(edu) + mstat + smoke + hibp + diab + hearte + stroke + bmi + adl + eurod + factor(rnd)"
cov1 <- "age + factor(edu) + mstat + smoke + hibp + diab + hearte + stroke + bmi + adl + eurod"
ext <- function(m, term, lab, dat, oc) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(NULL)
  b <- co[term, 1]; se <- co[term, 2]
  data.table(项 = lab, OR = exp(b), 下限 = exp(b - 1.96*se), 上限 = exp(b + 1.96*se),
             p = co[term, 4], n = length(m$y), 事件 = sum(m$y == 1, na.rm = TRUE))
}

res <- list()
for (oc in c("fall", "fxinc")) {
  d <- subset(des, !is.na(get(oc)) & !is.na(grip) & !is.na(age))
  m <- svyglm(as.formula(paste(oc, "~ grip_z +", cov)), design = d, family = quasibinomial())
  nm <- c(fall = "跌倒（近 6 个月）", fxinc = "新发髋骨骨折")[[oc]]
  res[[oc]] <- ext(m, "grip_z", paste0("SHARE：握力（每 SD=", round(sd_g, 1), "kg）→ ", nm), d, oc)
  cat(sprintf("\n[SHARE] %-8s n=%6d 事件=%5d OR=%.3f (%.3f-%.3f) p=%.4f\n", oc, res[[oc]]$n,
              res[[oc]]$事件, res[[oc]]$OR, res[[oc]]$下限, res[[oc]]$上限, res[[oc]]$p))
}
per <- rbindlist(lapply(seq_along(grip_w), function(k) {
  d <- subset(des, rnd == k & !is.na(fall) & !is.na(grip) & !is.na(age))
  m <- tryCatch(svyglm(as.formula(paste("fall ~ grip_z +", cov1)), design = d, family = quasibinomial()),
                error = function(e) { cat("  窗口", k, "失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients["grip_z", ]
  data.table(窗口 = paste0("w", grip_w[k], "→w", fall_w[k]), n = length(m$y),
                          事件 = sum(m$y == 1), OR = exp(co[1]),
             下限 = exp(co[1] - 1.96*co[2]), 上限 = exp(co[1] + 1.96*co[2]), p = co[4])
}))
print(per)

cat("\n--- 阳性对照（pooled fall 模型）---\n")
dc <- subset(des, !is.na(fall) & !is.na(grip) & !is.na(age))
mp <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = dc, family = quasibinomial())
ct <- summary(mp)$coefficients
print(round(ct[intersect(c("age", "bmi", "adl", "eurod", "diab", "stroke"), rownames(ct)), ], 4))

RES <- rbindlist(res, fill = TRUE)
## ---------- 敏感性：只保留暴露波未跌倒者 ----------
cat("\n--- 敏感性分析：排除基线已跌倒者 ---\n")
sen <- svydesign(ids = ~psu, strata = ~strat, weights = ~wt, data = as.data.frame(lg[fall0 == 0]), nest = TRUE)
ds <- subset(sen, !is.na(fall) & !is.na(grip) & !is.na(age))
ms <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = ds, family = quasibinomial())
cs <- summary(ms)$coefficients["grip_z", ]
cat(sprintf("[SHARE-敏感性] fall n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(ds), sum(ds$variables$fall == 1), exp(cs[1]), exp(cs[1]-1.96*cs[2]), exp(cs[1]+1.96*cs[2]), cs[4]))
## ---------- 敏感性 2：骨折改用"自上次访问以来"的直接新发定义 ----------
cat("\n--- 敏感性 2：骨折改用 rXhip（自上次访问以来）---\n")
d2 <- subset(des, !is.na(fx_sw) & !is.na(grip) & !is.na(age))
m2 <- svyglm(as.formula(paste("fx_sw ~ grip_z +", cov)), design = d2, family = quasibinomial())
c2 <- summary(m2)$coefficients["grip_z", ]
cat(sprintf("[SHARE-敏感性2] fx_sw n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(d2), sum(d2$variables$fx_sw == 1), exp(c2[1]),
            exp(c2[1]-1.96*c2[2]), exp(c2[1]+1.96*c2[2]), c2[4]))

SEN <- rbind(
  data.table(分析 = "SHARE-仅基线未跌倒者", 结局 = "跌倒", n = nrow(ds), 事件 = sum(ds$variables$fall == 1),
             OR = exp(cs[1]), 下限 = exp(cs[1]-1.96*cs[2]), 上限 = exp(cs[1]+1.96*cs[2]), p = cs[4]),
  data.table(分析 = "SHARE-骨折改用直接新发定义", 结局 = "新发髋骨骨折", n = nrow(d2),
             事件 = sum(d2$variables$fx_sw == 1), OR = exp(c2[1]),
             下限 = exp(c2[1]-1.96*c2[2]), 上限 = exp(c2[1]+1.96*c2[2]), p = c2[4]))
fwrite(SEN, file.path(out, "res", "13_SHARE_敏感性.csv"))

fwrite(RES, file.path(out, "res", "11_SHARE_主结果.csv"))
fwrite(per, file.path(out, "res", "12_SHARE_分窗口.csv"))
cat("\n===== SHARE 汇总 =====\n"); print(RES); print(per)

## ---------- 设计效应（2026-09-16 新增）----------
source(file.path(out, "00_deff_helper.R"))
cc <- function(d, oc) d[!is.na(d[[oc]]) & !is.na(d$grip) & !is.na(d$age) & !is.na(d$wt) & d$wt > 0, ]
DEFF <- rbindlist(list(
  survey_deff(as.data.frame(cc(copy(lg), "fall")), paste("fall ~ grip_z +", cov),
              "grip_z", "SHARE 跌倒（主模型）", "wt", psu = "psu", strata = "strat", 事件列 = "fall"),
  survey_deff(as.data.frame(cc(copy(lg), "fxinc")), paste("fxinc ~ grip_z +", cov),
              "grip_z", "SHARE 新发髋骨骨折", "wt", psu = "psu", strata = "strat", 事件列 = "fxinc"),
  survey_deff(as.data.frame(cc(copy(lg), "fx_sw")), paste("fx_sw ~ grip_z +", cov),
              "grip_z", "SHARE 骨折（直接新发定义）", "wt", psu = "psu", strata = "strat", 事件列 = "fx_sw")
), fill = TRUE)
fwrite(DEFF, file.path(out, "res", "32_SHARE_设计效应.csv"))
cat("\n===== SHARE 设计效应 =====\n")
print(DEFF[, .(项, n, 事件, OR_完整设计, SE_完整设计, SE_简单随机, DEFF_总体, 有效样本量, 含聚类, PSU数, 层数)])

## ---------- 敏感性 3：权重在国内标准化（2026-09-16 新增）----------
## 动机：SHARE 各版权重按国家人口规模标定，直接把 28 国权重放在一起，
##       相当于让德/法/意等大国主导合并估计，权重 CV 高达 1.55 → DEFF≈3.4。
##       本敏感性把权重在【国家内】标准化（各国权重和=该国样本量），
##       即"跨国按样本量加权"，用于检验结论是否被国家权重结构驱动。
## 注：n / 事件 / 权重CV 均按 svyglm 实际使用样本（含协变量完整）统计。
lg[, wt_norm := wt / sum(wt) * .N, by = country]
dn <- svydesign(ids = ~psu, strata = ~strat, weights = ~wt_norm,
                data = as.data.frame(lg), nest = TRUE)
dn2 <- subset(dn, !is.na(fall) & !is.na(grip) & !is.na(age))
mn2 <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = dn2, family = quasibinomial())
cn2 <- summary(mn2)$coefficients["grip_z", ]
## 2026-09-17 修正：n 必须取模型实际样本量（svyglm 会静默丢弃协变量缺失行）
## 2026-09-17 修正：n 必须取模型实际样本量（svyglm 会静默丢弃协变量缺失行）
##   注意不能用 dn2$variables$wt_norm[keep2] —— keep2 长度短于设计行数会被 R 循环补齐
wn <- as.numeric(mn2$prior.weights); n2 <- length(wn)
ev2 <- sum(mn2$y == 1)
kish2 <- n2 * sum(wn^2) / sum(wn)^2
cat(sprintf("\n[SHARE-敏感性3·权重国内标准化] fall n=%d OR=%.3f (%.3f-%.3f) p=%.4f；权重CV=%.3f DEFF_Kish=%.2f\n",
            n2, exp(cn2[1]), exp(cn2[1]-1.96*cn2[2]), exp(cn2[1]+1.96*cn2[2]), cn2[4],
            sqrt(kish2 - 1), kish2))
fwrite(data.table(分析 = "SHARE-权重国内标准化（跨国按样本量加权）", 结局 = "跌倒", n = n2,
                  事件 = ev2, OR = exp(cn2[1]),
                  下限 = exp(cn2[1]-1.96*cn2[2]), 上限 = exp(cn2[1]+1.96*cn2[2]), p = cn2[4],
                  权重CV = sqrt(kish2 - 1), DEFF_权重Kish = kish2),
       file.path(out, "res", "33_SHARE_权重标准化敏感性.csv"))
## ---------- 导出主模型分析样本（Table 1 用；2026-09-16 新增）----------
t1dir <- file.path(out, "t1"); dir.create(t1dir, showWarnings = FALSE, recursive = TRUE)
mv <- all.vars(as.formula(paste0("~", cov)))
need <- unique(c("id","rnd","country","grip","wt","strat","psu", all.vars(as.formula(paste0("fall ~ grip_z +", cov)))))
dd <- as.data.frame(des$variables)
dd <- dd[complete.cases(dd[, need]), ]
saveRDS(dd, file.path(t1dir, "样本_SHARE.rds"))
cat("\n[Table1] SHARE 模型分析样本 n =", nrow(dd), " 去重人数 =", uniqueN(dd$id),
    " 平均窗口数 =", round(nrow(dd)/uniqueN(dd$id), 2), "\n")
print(table(dd$edu, useNA = "ifany"))
