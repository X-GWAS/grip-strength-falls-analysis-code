## 70_设计效应汇总与合并.R —— 四队列 survey 设计校正后的设计效应汇总 + 合并估计重跑
## 输入：30_ELSA / 31_HRS / 32_SHARE 设计效应表；CHARLS res/12_CHARLS_设计效应.csv
## 输出：res/40_四队列设计效应汇总.csv、41_设计校正前后对照.csv、42_合并估计_设计校正.csv
suppressMessages({ library(data.table) })
W  <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
out <- file.path(W, "外部复现-ELSA与HRS")

R <- function(f) fread(f)
e <- R(file.path(out, "res", "30_ELSA_设计效应.csv"));  e[, 队列 := "ELSA（英国）"]
h <- R(file.path(out, "res", "31_HRS_设计效应.csv"));   h[, 队列 := "HRS（美国）"]
s <- R(file.path(out, "res", "32_SHARE_设计效应.csv")); s[, 队列 := "SHARE（欧洲）"]
c <- R(file.path(W, "新暴露-体力活动与握力轨迹/res", "12_CHARLS_设计效应.csv")); c[, 队列 := "CHARLS（中国）"]

DEFF <- rbindlist(list(c, e, h, s), fill = TRUE)
DEFF[, 结局 := fifelse(grepl("跌倒受伤", 项), "跌倒受伤",
               fifelse(grepl("骨折", 项), "新发髋骨骨折", "跌倒（任意）"))]
DEFF <- DEFF[!(队列 == "SHARE（欧洲）" & 结局 == "新发髋骨骨折" & grepl("直接新发", 项))]
DEFF[, `:=`(OR_设计校正 = OR_完整设计, 下限 = exp(log(OR_完整设计) - 1.96 * SE_完整设计),
            上限 = exp(log(OR_完整设计) + 1.96 * SE_完整设计),
            OR_仅加权 = OR_仅加权,   ## 【2026-09-17 修正】原写成 exp(OR_仅加权)，把已是 OR 的值再取一次指数，展示列错位（不影响 lower_w/upper_w 与合并结果）
            lower_w = exp(log(OR_仅加权) - 1.96 * SE_仅加权),
            upper_w = exp(log(OR_仅加权) + 1.96 * SE_仅加权),
            SE膨胀 = SE_完整设计 / SE_仅加权,
            方法 = fifelse(含聚类, "分层+PSU+权重", "仅权重（无可用分层/PSU）"))]
## CHARLS 的设计效应表由 Python 脚本产出，没有"有分层"这一列；CHARLS 确实是分层+PSU，按含聚类补齐
if (!"有分层" %in% names(DEFF)) DEFF[, 有分层 := NA]
DEFF[is.na(有分层), 有分层 := 含聚类]
DEFF[, 方法 := fifelse(含聚类 & 有分层, "分层+PSU+权重",
                fifelse(含聚类, "PSU（参与者）+权重（无分层）", "仅权重（无可用分层/PSU）"))]
DEFF[, 有分层 := NULL]
setcolorder(DEFF, c("队列","结局","项","n","事件","OR_设计校正","下限","上限","p_完整设计",
                    "DEFF_总体","有效样本量","DEFF_加权","DEFF_聚类","权重CV","DEFF_权重Kish",
                    "SE膨胀","方法","含聚类","PSU数","层数","OR_仅加权","SE_仅加权","SE_简单随机"))
setorder(DEFF, 结局, 队列)
fwrite(DEFF, file.path(out, "res", "40_四队列设计效应汇总.csv"))

## ---------- 设计校正前后对照 ----------
CMP <- DEFF[, .(队列, 结局, n, 事件,
                OR = OR_设计校正,
                旧_下限 = lower_w, 旧_上限 = upper_w, 旧_p = p_仅加权,
                新_下限 = 下限, 新_上限 = 上限, 新_p = p_完整设计,
                DEFF = DEFF_总体, 有效n = 有效样本量)]
CMP[, `:=`(显著性变化 = fifelse(旧_p < 0.05 & 新_p < 0.05, "仍显著",
                        fifelse(旧_p < 0.05 & 新_p >= 0.05, "校正后转为不显著",
                        fifelse(旧_p >= 0.05 & 新_p < 0.05, "校正后转为显著", "仍不显著"))))]
fwrite(CMP, file.path(out, "res", "41_设计校正前后对照.csv"))

## ---------- 合并估计（DerSimonian–Laird 随机效应），分别用旧/新 SE ----------
pf <- function(OR, LO, HI) {
  b <- log(OR); se <- (log(HI) - log(LO)) / (2 * 1.96)
  w <- 1 / se^2; bf <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bf)^2); df <- length(b) - 1
  C <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (Q - df) / C)
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  data.table(合并OR = exp(br), 下限 = exp(br - 1.96 * ser), 上限 = exp(br + 1.96 * ser),
             p = 2 * pnorm(-abs(br / ser)), I2 = I2, Q_p = pchisq(Q, df, lower.tail = FALSE),
             队列数 = df + 1, tau2 = tau2)
}
M <- rbindlist(lapply(c("跌倒（任意）","跌倒受伤","新发髋骨骨折"), function(oc) {
  d <- DEFF[结局 == oc]
  a <- pf(d$OR_设计校正, d$lower_w, d$upper_w); a[, `:=`(结局 = oc, SE口径 = "校正前（仅加权）")]
  b <- pf(d$OR_设计校正, d$下限, d$上限);       b[, `:=`(结局 = oc, SE口径 = "survey 设计校正后")]
  rbind(a, b)
}), fill = TRUE)
setcolorder(M, c("结局","SE口径","队列数","合并OR","下限","上限","p","I2","tau2","Q_p"))
setorder(M, 结局, SE口径)
fwrite(M, file.path(out, "res", "42_合并估计_设计校正.csv"))

cat("\n===== 四队列设计效应汇总 =====\n")
print(DEFF[, .(队列, 结局, n, 事件, OR_设计校正, 下限, 上限, p_完整设计,
               DEFF_总体, 有效样本量, 权重CV, 方法)])
cat("\n===== 设计校正前后对照 =====\n"); print(CMP)
cat("\n===== 合并估计 =====\n"); print(M)
