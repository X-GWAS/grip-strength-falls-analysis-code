## 83_Table2_亚组.R —— Table 2：主结果 ＋ 亚组分析（四队列统一协变量口径）
## 设计要点：
##   1) 亚组模型用【四队列统一核心协变量】age + factor(edu) + smoke + hibp + diab + bmi + adl
##      （ELSA/HRS/SHARE 另加 factor(rnd) 固定波次效应；CHARLS 为单窗口、无 rnd），
##      分层时把该分层变量本身从协变量里去掉；
##   2) 每层内先做队列特异的 survey 设计校正模型，再按 DerSimonian–Laird 随机效应合并；
##   3) 交互检验双轨：各队列内用 regTermTest 做设计校正 Wald 检验（报告逐队列 P），
##      跨队列用各层合并估计之间的 Q 检验（报告合并 P for interaction）。
suppressMessages({ library(data.table); library(survey); library(ggplot2); library(ragg) })
options(survey.lonely.psu = "adjust")
B  <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
T1 <- file.path(B, "外部复现-ELSA与HRS/t1")
RES<- file.path(B, "外部复现-ELSA与HRS/res")

prep <- function(rds) { d <- as.data.table(readRDS(rds)); if ("ID" %in% names(d) && !"id" %in% names(d)) setnames(d, "ID", "id"); d }
coh <- list(
  CHARLS = list(d = prep(file.path(B, "新暴露-体力活动与握力轨迹/t1/样本_CHARLS.rds")),
                ids = "communityID", strata = "stratum", wt = "wt3", lab = "CHARLS 中国", y = "fall_1518", rnd = FALSE),
  ELSA   = list(d = prep(file.path(T1, "样本_ELSA.rds")), ids = "id",    strata = NULL,     wt = "wt", lab = "ELSA 英国",   y = "fall",      rnd = TRUE),
  HRS    = list(d = prep(file.path(T1, "样本_HRS.rds")),  ids = "psu",  strata = "strat",  wt = "wt", lab = "HRS 美国",    y = "fall",      rnd = TRUE),
  SHARE  = list(d = prep(file.path(T1, "样本_SHARE.rds")),ids = "psu",  strata = "strat",  wt = "wt", lab = "SHARE 欧洲",  y = "fall",      rnd = TRUE)
)

## ---------- 分层变量定义（四队列统一口径）----------
## 教育：各队列「最低受教育层级」= 低（CHARLS raeduc_c<=2、ELSA raeduc_e==1、HRS/SHARE raeducl==1）
make_strata <- function(d, cn) {
  d <- copy(d)
  d[, age_g := factor(cut(age, c(-Inf, 65, 75, Inf), labels = c("<65", "65–74", "≥75")), levels = c("<65","65–74","≥75"))]
  d[, bmi_g := factor(cut(bmi, c(-Inf, 25, 30, Inf), labels = c("<25", "25–29.9", "≥30")), levels = c("<25","25–29.9","≥30"))]
  d[, cmd := fifelse(is.na(hibp) | is.na(diab), NA_integer_, hibp + diab)]
  d[, cmd_g := factor(fifelse(cmd >= 2, "≥2", as.character(cmd)), levels = c("0", "1", "≥2"))]
  d[, adl_g := factor(fifelse(adl >= 1, "有 ADL 受限", "无 ADL 受限"), levels = c("无 ADL 受限", "有 ADL 受限"))]
  d[, smoke_g := factor(fifelse(smoke == 1, "曾/现吸烟", "从不吸烟"), levels = c("从不吸烟", "曾/现吸烟"))]
  lo <- if (cn == "CHARLS") 2 else if (cn == "ELSA") 1 else 1
  d[, edu_g := factor(fifelse(edu <= lo & !is.na(edu), "低教育", "中/高教育"), levels = c("中/高教育", "低教育"))]
  d
}
STRAT <- list(
  "年龄" = list(v = "age_g", drop = "age"),
  "BMI" = list(v = "bmi_g", drop = "bmi"),
  "心血管代谢共病（高血压/糖尿病）" = list(v = "cmd_g", drop = c("hibp", "diab")),
  "ADL 受限" = list(v = "adl_g", drop = "adl"),
  "吸烟" = list(v = "smoke_g", drop = "smoke"),
  "教育" = list(v = "edu_g", drop = "edu"))
CORE <- c("age", "edu", "smoke", "hibp", "diab", "bmi", "adl")

## ---------- 随机效应合并（DL）----------
pf <- function(OR, LO, HI) {
  b <- log(OR); se <- (log(HI) - log(LO)) / (2 * 1.96)
  w <- 1 / se^2; bf <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bf)^2); df <- length(b) - 1
  C <- sum(w) - sum(w^2) / sum(w); tau2 <- if (df > 0) max(0, (Q - df) / C) else 0
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0 && df > 0) max(0, (Q - df) / Q) * 100 else 0
  list(or = exp(br), lo = exp(br - 1.96 * ser), hi = exp(br + 1.96 * ser),
       p = 2 * pnorm(-abs(br / ser)), I2 = I2, Q = Q)
}

rows <- list(); irows <- list()
for (cn in names(coh)) {
  z <- coh[[cn]]; d <- make_strata(z$d, cn)
  cat("\n##########", z$lab, "  n =", nrow(d), "##########\n")
  keep <- c(z$y, "grip_z", CORE, unlist(lapply(STRAT, `[[`, "v")), if (z$rnd) "rnd")
  dd <- d[complete.cases(d[, ..keep])]
  des <- if (is.null(z$strata)) {
    svydesign(ids = as.formula(paste0("~", z$ids)), weights = as.formula(paste0("~", z$wt)), data = as.data.frame(dd))
  } else {
    svydesign(ids = as.formula(paste0("~", z$ids)), strata = as.formula(paste0("~", z$strata)),
              weights = as.formula(paste0("~", z$wt)), data = as.data.frame(dd), nest = TRUE)
  }
  cov_full <- paste(c(CORE, if (z$rnd) "factor(rnd)"), collapse = " + ")
  ## 总体（统一协变量口径）
  m0 <- svyglm(as.formula(paste(z$y, "~ grip_z +", cov_full)), design = des, family = quasibinomial())
  co <- summary(m0)$coefficients
  cat(sprintf("总体（统一协变量）OR=%.3f (%.3f-%.3f) p=%.3g n=%d\n",
              exp(co["grip_z",1]), exp(co["grip_z",1]-1.96*co["grip_z",2]), exp(co["grip_z",1]+1.96*co["grip_z",2]),
              co["grip_z",4], length(m0$y)))
  rows[[paste(cn,"总体")]] <- data.table(队列 = z$lab, 分层变量 = "总体", 层 = "全部", n = length(m0$y),
    事件 = sum(m0$y), OR = exp(co["grip_z",1]), 下限 = exp(co["grip_z",1]-1.96*co["grip_z",2]),
    上限 = exp(co["grip_z",1]+1.96*co["grip_z",2]), p = co["grip_z",4])
  ## 逐分层变量
  for (sn in names(STRAT)) {
    sv <- STRAT[[sn]]$v; dropv <- STRAT[[sn]]$drop
    covs <- paste(c(setdiff(CORE, dropv), if (z$rnd) "factor(rnd)"), collapse = " + ")
    levs <- levels(dd[[sv]])
    for (lv in levs) {
      sub <- subset(des, eval(parse(text = paste0("dd$", sv, " == '", lv, "'"))))
      if (nrow(sub$variables) < 60) { cat(sprintf("  [跳过] %s %s 层 n=%d 太少\n", sn, lv, nrow(sub$variables))); next }
      m <- tryCatch(svyglm(as.formula(paste(z$y, "~ grip_z +", covs)), design = sub, family = quasibinomial()),
                    error = function(e) NULL)
      if (is.null(m)) { cat(sprintf("  [失败] %s %s\n", sn, lv)); next }
      c2 <- summary(m)$coefficients
      if (!"grip_z" %in% rownames(c2)) next
      ev <- sum(m$y)
      cat(sprintf("  %-28s %-10s OR=%.3f (%.3f-%.3f) p=%.3g n=%d 事件=%d\n", sn, lv,
                  exp(c2["grip_z",1]), exp(c2["grip_z",1]-1.96*c2["grip_z",2]), exp(c2["grip_z",1]+1.96*c2["grip_z",2]),
                  c2["grip_z",4], length(m$y), ev))
      rows[[paste(cn, sn, lv)]] <- data.table(队列 = z$lab, 分层变量 = sn, 层 = lv, n = length(m$y), 事件 = ev,
        OR = exp(c2["grip_z",1]), 下限 = exp(c2["grip_z",1]-1.96*c2["grip_z",2]),
        上限 = exp(c2["grip_z",1]+1.96*c2["grip_z",2]), p = c2["grip_z",4])
    }
    ## 队列内交互检验（设计校正 Wald）
    mi <- tryCatch(svyglm(as.formula(paste(z$y, "~ grip_z *", sv, "+", covs)), design = des, family = quasibinomial()),
                   error = function(e) NULL)
    if (!is.null(mi)) {
      tt <- tryCatch(regTermTest(mi, as.formula(paste0("~ grip_z:", sv))), error = function(e) NULL)
      if (!is.null(tt)) {
        cat(sprintf("  【交互】%s 设计校正 Wald P = %.4f\n", sn, tt$p))
        irows[[paste(cn, sn)]] <- data.table(队列 = z$lab, 分层变量 = sn, 交互P = tt$p)
      }
    }
  }
}
R <- rbindlist(rows); I <- rbindlist(irows)
fwrite(R, file.path(RES, "60_Table2_分层估计.csv"))
fwrite(I, file.path(RES, "61_Table2_队列内交互P.csv"))

## ---------- 按层合并 ----------
pool <- R[分层变量 != "总体", {
  r <- pf(OR, 下限, 上限)
  .(队列数 = .N, 合并OR = r$or, 下限 = r$lo, 上限 = r$hi, p = r$p, I2 = r$I2, Q = r$Q)
}, by = .(分层变量, 层)]
## 总体合并（统一协变量口径）＋ 主模型口径对照
tot <- pf(R[分层变量 == "总体"]$OR, R[分层变量 == "总体"]$下限, R[分层变量 == "总体"]$上限)
cat(sprintf("\n总体合并（统一协变量口径）：OR=%.3f (%.3f-%.3f) p=%.3g  I2=%.1f%%\n", tot$or, tot$lo, tot$hi, tot$p, tot$I2))

## 跨层交互（合并层间 Q 检验）
inter <- pool[, {
  b <- log(合并OR); se <- (log(上限) - log(下限)) / (2 * 1.96)
  w <- 1/se^2; bf <- sum(w*b)/sum(w); Q <- sum(w*(b-bf)^2); df <- .N - 1
  .(层数 = .N, 层间Q = Q, 层间P = pchisq(Q, df, lower.tail = FALSE))
}, by = 分层变量]
print(inter)
fwrite(pool, file.path(RES, "62_Table2_层合并.csv"))
fwrite(inter, file.path(RES, "63_Table2_层间交互P.csv"))
cat("\n完成：res/60_Table2_分层估计.csv / 61_队列内交互P / 62_层合并 / 63_层间交互P\n")
