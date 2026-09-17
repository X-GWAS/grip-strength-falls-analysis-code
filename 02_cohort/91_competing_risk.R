## 91_竞争风险分析.R —— 死亡竞争风险对握力—跌倒关联的影响
## 四套口径（保持协变量、权重、设计完全一致，只改"死亡怎么处理"）：
##   A 原因别风险（cause-specific）＝现有主分析口径：只用观察到结局的人
##   B 亚分布／累计发生率（single-interval Fine-Gray 等价）：把结局前死亡者保留在风险集、记为"未跌倒"
##   C 最坏情形：把结局前死亡者记为"已跌倒"
##   D/E 极端边界：把"其他缺失"也算进来，全部记为未跌倒（D，乐观）／全部记为已跌倒（E，悲观）
suppressMessages({ library(data.table); library(survey) })
options(survey.lonely.psu = "adjust")
B  <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/竞争风险死亡"
OUT<- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res"

coh <- list(
  CHARLS = list(f = "elig_CHARLS.rds", lab = "CHARLS 中国", sd = 7.5, ids = "psu", strata = "strat", wt = "wt", rnd = FALSE),
  ELSA   = list(f = "elig_ELSA.rds",   lab = "ELSA 英国",   sd = 6.2, ids = "id",   strata = NULL,     wt = "wt", rnd = TRUE),
  HRS    = list(f = "elig_HRS.rds",    lab = "HRS 美国",    sd = 6.0, ids = "psu", strata = "strat",  wt = "wt", rnd = TRUE),
  ## 【2026-09-17 v2】SHARE 在本竞争风险样本里只保留了 id 与 wt（无 psu/strat），
  ##   原来用 ids = "1" 会把同一名女性的多个随访窗口当成独立观测、低估标准误；
  ##   与 ELSA 同理改为按参与者聚类（ids = "id"），这是数据可得性下的最小必要修正。
  SHARE  = list(f = "elig_SHARE.rds",  lab = "SHARE 欧洲",  sd = 7.1, ids = "id",  strata = NULL,     wt = "wt", rnd = TRUE)
)
CORE <- c("age", "edu", "smoke", "hibp", "diab", "bmi", "adl")

pf <- function(OR, LO, HI) {
  b <- log(OR); se <- (log(HI) - log(LO)) / (2 * 1.96)
  w <- 1 / se^2; bf <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bf)^2); df <- length(b) - 1
  C <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (Q - df) / C)
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(or = exp(br), lo = exp(br - 1.96 * ser), hi = exp(br + 1.96 * ser),
       p = 2 * pnorm(-abs(br / ser)), I2 = I2)
}

rows <- list(); desc <- list()
for (cn in names(coh)) {
  z <- coh[[cn]]
  d <- as.data.table(readRDS(file.path(B, "t1", z$f)))
  d <- d[!is.na(grip) & !is.na(wt) & wt > 0]
  d[, grip_z := grip / z$sd]
  ## 三分类状态
  d[, y := fall]
  d[is.na(fall) & died == 1, y_death := 1L]
  cat(sprintf("\n########## %s ##########\n符合条件 %d；观测 %d；死亡 %d；其他缺失 %d\n",
              z$lab, nrow(d), sum(d$obs), sum(d$died), sum(d$othermiss)))
  desc[[cn]] <- data.table(队列 = z$lab, 符合条件 = nrow(d), 观察到结局 = sum(d$obs),
                           结局前死亡 = sum(d$died), 其他缺失 = sum(d$othermiss),
                           死亡占比 = sprintf("%.1f%%", 100 * sum(d$died) / nrow(d)))
  covs <- paste(c(CORE, if (z$rnd) "factor(rnd)"), collapse = " + ")
  keep <- c("y", "grip_z", CORE, if (z$rnd) "rnd")
  fit <- function(x, tag) {
    x <- x[complete.cases(x[, ..keep])]
    if (nrow(x) < 50) return(NULL)
    des <- if (is.null(z$strata)) svydesign(ids = as.formula(paste0("~", z$ids)),
                                            weights = as.formula(paste0("~", z$wt)), data = as.data.frame(x))
           else svydesign(ids = as.formula(paste0("~", z$ids)), strata = as.formula(paste0("~", z$strata)),
                          weights = as.formula(paste0("~", z$wt)), data = as.data.frame(x), nest = TRUE)
    m <- tryCatch(svyglm(as.formula(paste("y ~ grip_z +", covs)), design = des, family = quasibinomial()),
                  error = function(e) { cat("   [失败]", tag, conditionMessage(e), "\n"); NULL })
    if (is.null(m)) return(NULL)
    c2 <- summary(m)$coefficients["grip_z", ]
    cat(sprintf("  %-34s OR=%.3f (%.3f-%.3f) p=%.3g  n=%d 事件=%d\n", tag,
                exp(c2[1]), exp(c2[1] - 1.96 * c2[2]), exp(c2[1] + 1.96 * c2[2]), c2[4],
                length(m$y), sum(m$y)))
    data.table(队列 = z$lab, 口径 = tag, n = length(m$y), 事件 = sum(m$y),
               OR = exp(c2[1]), 下限 = exp(c2[1] - 1.96 * c2[2]), 上限 = exp(c2[1] + 1.96 * c2[2]), p = c2[4])
  }
  xA <- copy(d); xA[is.na(fall), y := NA_integer_]
  rows[[paste(cn, "A")]] <- fit(xA, "A 原因别风险（仅存活观察者）")
  xB <- copy(d); xB[, y := fifelse(!is.na(fall), as.integer(fall == 1), fifelse(died == 1, 0L, NA_integer_))]
  rows[[paste(cn, "B")]] <- fit(xB, "B 亚分布（死亡记为未跌倒）")
  xC <- copy(d); xC[, y := fifelse(!is.na(fall), as.integer(fall == 1), fifelse(died == 1, 1L, NA_integer_))]
  rows[[paste(cn, "C")]] <- fit(xC, "C 最坏情形（死亡记为跌倒）")
  xD <- copy(d); xD[, y := fifelse(!is.na(fall), as.integer(fall == 1), 0L)]
  rows[[paste(cn, "D")]] <- fit(xD, "D 全部缺失记为未跌倒（乐观）")
  xE <- copy(d); xE[, y := fifelse(!is.na(fall), as.integer(fall == 1), 1L)]
  rows[[paste(cn, "E")]] <- fit(xE, "E 全部缺失记为跌倒（悲观）")
}
R <- rbindlist(rows); D <- rbindlist(desc)
fwrite(R, file.path(OUT, "70_竞争风险_逐队列.csv"))
fwrite(D, file.path(OUT, "71_竞争风险_人数构成.csv"))

cat("\n\n===== 合并（DerSimonian–Laird 随机效应）=====\n")
P <- R[, { r <- pf(OR, 下限, 上限)
           .(队列数 = .N, 合并OR = r$or, 下限 = r$lo, 上限 = r$hi, p = r$p, I2 = r$I2) }, by = 口径]
print(P)
fwrite(P, file.path(OUT, "72_竞争风险_合并.csv"))
cat("\n完成：res/70_竞争风险_逐队列.csv / 71_人数构成 / 72_合并\n")
