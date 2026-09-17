## 10_elsa.R —— ELSA 平行复现：握力 → 跌倒 / 跌倒受伤 / 髋骨骨折
## 设计：暴露在前（w2/4/6/8 握力），结局在后（w3/5/7/9 跌倒）；女性 ≥50 岁
suppressMessages({
  source("/Users/mac/.codex/knowledge/数据库资产/db_paths.R")
  library(haven); library(data.table); library(survey)
})
options(survey.lonely.psu = "adjust")

out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS"
dir.create(file.path(out, "res"), showWarnings = FALSE, recursive = TRUE)

## ---------- 变量清单 ----------
grip_w <- c(2, 4, 6, 8); fall_w <- c(3, 5, 7, 9)
vars <- c("idauniq", "ragender", "raeduc_e", "raeducl",
          paste0("r", 1:9, "agey"),
          paste0("r", grip_w, "lgrip1"), paste0("r", grip_w, "lgrip2"),
          paste0("r", grip_w, "rgrip1"), paste0("r", grip_w, "rgrip2"),
          paste0("r", grip_w, "gripsum"),
          paste0("r", fall_w, "fall"), paste0("r", fall_w, "fallinj"),
          paste0("r", c(2,4,6,8), "fall"), "r4fall1y",   # 暴露波跌倒，用于基线未跌倒者敏感性分析
          "r4fall1y", "r4fallinj",
          paste0("r", 1:9, "hipe"),
          paste0("r", c(2,4,6,8), "mbmi"),
          paste0("r", c(2,4,6,8), "smokev"),
          paste0("r", c(2,4,6,8), "drink"),
          paste0("r", c(2,4,6,8), "adlwa"),
          paste0("r", c(2,4,6,8), "cesd"),
          paste0("r", c(2,4,6,8), "hibpe"), paste0("r", c(2,4,6,8), "diabe"),
          paste0("r", c(2,4,6,8), "hearte"), paste0("r", c(2,4,6,8), "stroke"),
          paste0("r", c(2,4,6,8), "mstat"),
          paste0("r", c(2,4,6,8), "cwtresp"), paste0("r", c(3,5,7,9), "cwtresp"))

cat("读取 ELSA harmonized ...\n")
D <- as.data.table(read_dta(DB$file$elsa_h, col_select = any_of(vars)))
cat("行数:", nrow(D), " 列数:", ncol(D), "\n")

## 握力：取优势手（若无则左右手第 1 次均值）
for (w in grip_w) {
  g1 <- paste0("r", w, "lgrip1"); g2 <- paste0("r", w, "rgrip1")
  a <- if (g1 %in% names(D)) as.numeric(D[[g1]]) else rep(NA_real_, nrow(D))
  b <- if (g2 %in% names(D)) as.numeric(D[[g2]]) else rep(NA_real_, nrow(D))
  a[a <= 0 | a > 100] <- NA; b[b <= 0 | b > 100] <- NA
  D[[paste0("grip", w)]] <- rowMeans(cbind(a, b), na.rm = TRUE)
  D[is.nan(get(paste0("grip", w))), (paste0("grip", w)) := NA_real_]
}
for (w in fall_w) {
  f <- paste0("r", w, "fall"); fi <- paste0("r", w, "fallinj")
  D[[paste0("fall", w)]] <- if (f %in% names(D)) fifelse(is.na(D[[f]]), NA_integer_, as.integer(D[[f]] == 1)) else NA_integer_
  ## fallinj 只对跌倒者提问：未跌倒记 0，跌倒者取其值
  fv <- D[[paste0("fall", w)]]
  iv <- if (fi %in% names(D)) as.integer(D[[fi]] == 1) else rep(NA_integer_, nrow(D))
  D[[paste0("fallinj", w)]] <- fifelse(is.na(fv), NA_integer_, fifelse(fv == 0L, 0L, iv))
}
## 注：w4 的跌倒题是"过去一年"（r4fall1y），与其余波"过去两年"口径不同，
##     本分析不使用 w4 作结局，故不纳入（避免口径混用）。

## 髋骨骨折：ever → 新发（严格口径：必须前一波 hipe==0 且后一波 hipe 非缺失）
## 修正记录（2026-09-16）：原写法把"基线缺失但随访 hipe==1"也算成新发，
## 会把早年的陈旧骨折当成本次随访新发，SHARE 上实测污染 3,405 例 vs 真新发 2,224 例。
for (k in seq_along(grip_w)) {
  a <- paste0("r", grip_w[k], "hipe"); b <- paste0("r", fall_w[k], "hipe")
  av <- if (a %in% names(D)) as.integer(D[[a]] == 1) else NA_integer_
  bv <- if (b %in% names(D)) as.integer(D[[b]] == 1) else NA_integer_
  D[[paste0("fxinc", fall_w[k])]] <- fifelse(!is.na(av) & !is.na(bv) & av == 0, bv, NA_integer_)
}

## ---------- 长表：四个暴露-结局窗口 ----------
lg <- rbindlist(lapply(seq_along(grip_w), function(k) {
  gw <- grip_w[k]; fw <- fall_w[k]
  data.table(
    id     = D$idauniq,
    rnd    = k,
    grip_wave = gw, fall_wave = fw,
    grip   = D[[paste0("grip", gw)]],
    fall   = D[[paste0("fall", fw)]],
    fallinj= D[[paste0("fallinj", fw)]],
    fxinc  = D[[paste0("fxinc", fw)]],
    ## 暴露波是否已跌倒（w4 用 r4fall1y，口径为"过去一年"，略窄）
    fall0  = if (gw == 4) fifelse(is.na(D$r4fall1y), NA_integer_, as.integer(D$r4fall1y == 1))
             else fifelse(is.na(D[[paste0("r", gw, "fall")]]), NA_integer_, as.integer(D[[paste0("r", gw, "fall")]] == 1)),
    age    = D[[paste0("r", gw, "agey")]],
    female = as.integer(D$ragender == 2),
    edu    = as.numeric(D$raeduc_e),
    bmi    = as.numeric(D[[paste0("r", gw, "mbmi")]]),
    smoke  = as.integer(D[[paste0("r", gw, "smokev")]] == 1),
    drink  = as.integer(D[[paste0("r", gw, "drink")]] == 1),
    adl    = as.numeric(D[[paste0("r", gw, "adlwa")]]),
    cesd   = as.numeric(D[[paste0("r", gw, "cesd")]]),
    hibp   = as.integer(D[[paste0("r", gw, "hibpe")]] == 1),
    diab   = as.integer(D[[paste0("r", gw, "diabe")]] == 1),
    hearte = as.integer(D[[paste0("r", gw, "hearte")]] == 1),
    stroke = as.integer(D[[paste0("r", gw, "stroke")]] == 1),
    mstat  = as.integer(D[[paste0("r", gw, "mstat")]] %in% c(1, 2)),
    ## 权重统一取【暴露波】横断面个体权重（与 HRS 一致）；结局波权重版本另作敏感性
    wt     = as.numeric(D[[paste0("r", gw, "cwtresp")]]),
    wt_out = as.numeric(D[[paste0("r", fw, "cwtresp")]])
  )
}))
lg <- lg[female == 1 & !is.na(age) & age >= 50]
cat("\n窗口记录数:", nrow(lg), " 人数:", uniqueN(lg$id), "\n")
print(lg[, .(n = .N, 握力非缺失 = sum(!is.na(grip)), 握力均值 = round(mean(grip, na.rm = TRUE), 2),
             跌倒非缺失 = sum(!is.na(fall)), 跌倒事件 = sum(fall == 1, na.rm = TRUE),
             受伤事件 = sum(fallinj == 1, na.rm = TRUE), 骨折事件 = sum(fxinc == 1, na.rm = TRUE)), by = rnd])

## 标准化握力（以第 1 轮 SD 为准，便于与 CHARLS 对照）
sd_g <- sd(lg[rnd == 1]$grip, na.rm = TRUE)
lg[, grip_z := grip / sd_g]
cat("握力 SD (w2):", round(sd_g, 2), "kg\n")

## ---------- 模型 ----------
## 抽样设计（2026-09-16 核查结论）：
##   ELSA harmonized 文件里的 rXsecure 实为"工作是否有保障"（一道就业题），不是初级抽样单元；
##   rXstrat（hse stratification variable）覆盖率仅 47%–58%，且各波编码不一致（w2 为 1–100 的层号，
##   w4 起变成 1 万以上的高基数 ID），无法作为可用的分层变量。
##   因此 ELSA 无法做完整设计校正。HRS 与 SHARE 另有可用的分层与 PSU，见各自脚本。
##
## 【2026-09-17 v2 修正】原来把 ELSA 写成 ids = ~1（仅权重），等于把同一名女性在不同随访
##   窗口的多条记录当作互相独立的观测。ELSA 每人平均贡献 2.15 个窗口，必须处理人内相关。
##   实跑对照：ids = ~1 时 SE = 0.03150；以参与者为最后阶段聚类单位时 SE = 0.03322（+5.5%）。
##   同时对 HRS/SHARE 做了同样测试：按 PSU 聚类的 SE 等于或大于按人聚类（HRS 0.03059 vs
##   0.02925），说明有 PSU 的队列其设计校正已涵盖人内相关，无需改动。
##   故 v2 起 ELSA 统一采用 ids = ~id（等价于按人聚类的稳健标准误）。
##   仍需在局限性说明：ELSA 真实的分层与 PSU 未随数据发布，该修正只是最小必要修正，
##   其设计效应仍被低估。
lg <- lg[!is.na(wt) & wt > 0]
des <- svydesign(ids = ~id, weights = ~wt, data = as.data.frame(lg))
cov <- "age + factor(edu) + mstat + smoke + drink + hibp + diab + hearte + stroke + bmi + adl + cesd + factor(rnd)"
cov_1 <- "age + factor(edu) + mstat + smoke + drink + hibp + diab + hearte + stroke + bmi + adl + cesd"
ext <- function(m, term, lab) {
  co <- summary(m)$coefficients
  if (!term %in% rownames(co)) return(NULL)
  data.table(项 = lab, OR = exp(co[term,1]), 下限 = exp(co[term,1]-1.96*co[term,2]),
             上限 = exp(co[term,1]+1.96*co[term,2]), p = co[term,4])
}

res <- list()
for (oc in c("fall", "fallinj", "fxinc")) {
  d <- subset(des, !is.na(get(oc)) & !is.na(grip) & !is.na(age))
  m <- svyglm(as.formula(paste(oc, "~ grip_z +", cov)), design = d, family = quasibinomial())
  co <- summary(m)$coefficients["grip_z", ]
  nm <- c(fall = "跌倒（任意）", fallinj = "跌倒受伤", fxinc = "新发髋骨骨折")[[oc]]
  res[[oc]] <- data.table(项 = paste0("ELSA：握力（每 SD=", round(sd_g,1), "kg）→ ", nm),
                          OR = exp(co[1]), 下限 = exp(co[1]-1.96*co[2]), 上限 = exp(co[1]+1.96*co[2]),
                          p = co[4], n = length(m$y), 事件 = sum(m$y == 1))
  cat(sprintf("\n[ELSA] %-10s n=%6d 事件=%5d  每 SD OR=%.3f (%.3f-%.3f) p=%.4f\n",
              oc, length(m$y), sum(m$y == 1), exp(co[1]),
              exp(co[1]-1.96*co[2]), exp(co[1]+1.96*co[2]), co[4]))
}

## 分轮次
per <- rbindlist(lapply(unique(lg$rnd), function(k) {
  d <- subset(des, rnd == k & !is.na(fall) & !is.na(grip) & !is.na(age))
  m <- tryCatch(svyglm(as.formula(paste("fall ~ grip_z +", cov_1)), design = d, family = quasibinomial()),
                error = function(e) { cat("  窗口", k, "拟合失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients["grip_z", ]
  data.table(窗口 = paste0("w", grip_w[k], "→w", fall_w[k]), n = length(m$y),
                          事件 = sum(m$y == 1), OR = exp(co[1]),
             下限 = exp(co[1]-1.96*co[2]), 上限 = exp(co[1]+1.96*co[2]), p = co[4])
}))
print(per)

## ---------- 阳性对照：结局能否被已知危险因素抓到 ----------
cat("\n--- 阳性对照（pooled fall 模型，同一样本）---\n")
dc <- subset(des, !is.na(fall) & !is.na(grip) & !is.na(age))
mp <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = dc, family = quasibinomial())
print(round(summary(mp)$coefficients[intersect(c("age","bmi","adl","diab","cesd"), rownames(summary(mp)$coefficients)), ], 4))

## ---------- 敏感性 1：只保留暴露波未跌倒者 ----------
cat("\n--- 敏感性 1：排除基线已跌倒者 ---\n")
sen <- svydesign(ids = ~id, weights = ~wt, data = as.data.frame(lg[!is.na(wt) & wt > 0 & fall0 == 0]))
ds <- subset(sen, !is.na(fall) & !is.na(grip) & !is.na(age))
ms <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = ds, family = quasibinomial())
cs <- summary(ms)$coefficients["grip_z", ]
cat(sprintf("[ELSA-敏感性1] fall n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(ds), sum(ds$variables$fall == 1), exp(cs[1]), exp(cs[1]-1.96*cs[2]), exp(cs[1]+1.96*cs[2]), cs[4]))

## ---------- 敏感性 2：改用结局波权重 ----------
cat("\n--- 敏感性 2：改用结局波权重 ---\n")
des2 <- svydesign(ids = ~id, weights = ~wt_out, data = as.data.frame(lg[!is.na(wt_out) & wt_out > 0]))
d2 <- subset(des2, !is.na(fall) & !is.na(grip) & !is.na(age))
m2 <- svyglm(as.formula(paste("fall ~ grip_z +", cov)), design = d2, family = quasibinomial())
c2 <- summary(m2)$coefficients["grip_z", ]
cat(sprintf("[ELSA-敏感性2] fall n=%d 事件=%d OR=%.3f (%.3f-%.3f) p=%.4f\n",
            nrow(d2), sum(d2$variables$fall == 1), exp(c2[1]), exp(c2[1]-1.96*c2[2]), exp(c2[1]+1.96*c2[2]), c2[4]))

SEN <- rbind(
  data.table(分析 = "ELSA-仅基线未跌倒者", 结局 = "跌倒", n = nrow(ds), 事件 = sum(ds$variables$fall == 1),
             OR = exp(cs[1]), 下限 = exp(cs[1]-1.96*cs[2]), 上限 = exp(cs[1]+1.96*cs[2]), p = cs[4]),
  data.table(分析 = "ELSA-结局波权重", 结局 = "跌倒", n = nrow(d2), 事件 = sum(d2$variables$fall == 1),
             OR = exp(c2[1]), 下限 = exp(c2[1]-1.96*c2[2]), 上限 = exp(c2[1]+1.96*c2[2]), p = c2[4]))
fwrite(SEN, file.path(out, "res", "06_ELSA_敏感性.csv"))

RES <- rbindlist(res, fill = TRUE)
fwrite(RES, file.path(out, "res", "01_ELSA_主结果.csv"))
fwrite(per, file.path(out, "res", "02_ELSA_分窗口.csv"))

## ---------- 设计效应（2026-09-16 新增）----------
source(file.path(out, "00_deff_helper.R"))
cc <- function(d, oc) d[!is.na(d[[oc]]) & !is.na(d$grip) & !is.na(d$age) & !is.na(d$wt) & d$wt > 0, ]
base_dat <- lg[!is.na(lg$wt) & lg$wt > 0, ]
DEFF <- rbindlist(list(
  survey_deff(as.data.frame(cc(copy(base_dat), "fall")),
              paste("fall ~ grip_z +", cov), "grip_z", "ELSA 跌倒（主模型）", "wt", psu = "id", 事件列 = "fall"),
  survey_deff(as.data.frame(cc(copy(base_dat), "fallinj")),
              paste("fallinj ~ grip_z +", cov), "grip_z", "ELSA 跌倒受伤", "wt", psu = "id", 事件列 = "fallinj"),
  survey_deff(as.data.frame(cc(copy(base_dat), "fxinc")),
              paste("fxinc ~ grip_z +", cov), "grip_z", "ELSA 新发髋骨骨折", "wt", psu = "id", 事件列 = "fxinc")
), fill = TRUE)
fwrite(DEFF, file.path(out, "res", "30_ELSA_设计效应.csv"))
cat("\n===== ELSA 设计效应 =====\n")
print(DEFF[, .(项, n, 事件, OR_完整设计, SE_完整设计, SE_简单随机, DEFF_总体, 有效样本量, 含聚类, PSU数, 层数)])
cat("\n===== ELSA 汇总 =====\n"); print(RES); print(per)
## ---------- 导出主模型分析样本（Table 1 用；2026-09-16 新增）----------
t1dir <- file.path(out, "t1"); dir.create(t1dir, showWarnings = FALSE, recursive = TRUE)
mv <- all.vars(as.formula(paste0("~", cov)))
need <- unique(c("id","rnd","grip","wt", all.vars(as.formula(paste0("fall ~ grip_z +", cov)))))
dd <- as.data.frame(des$variables)
dd <- dd[complete.cases(dd[, need]), ]
saveRDS(dd, file.path(t1dir, "样本_ELSA.rds"))
cat("\n[Table1] ELSA 模型分析样本 n =", nrow(dd), " 去重人数 =", uniqueN(dd$id),
    " 平均窗口数 =", round(nrow(dd)/uniqueN(dd$id), 2), "\n")
print(table(dd$edu, useNA = "ifany"))
