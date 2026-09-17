## 30_多队列对照与合并.R —— CHARLS / ELSA / HRS / SHARE 四队列 ＋ MR 对照与合并估计
## 【2026-09-16 更新】全部数值改为 survey 设计校正后的结果（svydesign + svyglm）：
##   CHARLS/ELSA/HRS/SHARE 均使用各自可得的最强设计（分层 + PSU + 权重；ELSA 无可用分层变量 → 仅权重）。
##   校正前（仅加权）的旧数值完整备份在 res_v1_仅权重/，对照见 res/41_设计校正前后对照.csv。
suppressMessages({ library(data.table) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS"

## ---------- 1. 主结果：跌倒（每 +1 SD 握力） ----------
## 【2026-09-17 v2】原来本节把四队列的 OR/CI/p 硬编码在脚本里，重跑上游脚本不会传导，
##   违反"数字必须由脚本现算"的规则（台账 C104）。现改为直接读取各队列的结果文件。
CHARLS_RES <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/新暴露-体力活动与握力轨迹/res/12_CHARLS_设计效应.csv"
MR_RES     <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/mr-finn/res/04_对照_旧分析_UKB同源.csv"

rd1 <- function(p, pat, orc = "OR", loc = "下限", hic = "上限", pc = "p") {
  d <- as.data.table(fread(p, encoding = "UTF-8"))
  r <- d[grepl(pat, 项)][1]
  if (nrow(r) == 0) stop("未匹配到: ", pat, " in ", p)
  list(OR = as.numeric(r[[orc]]), L = as.numeric(r[[loc]]), U = as.numeric(r[[hic]]),
       p = as.numeric(r[[pc]]),
       n = if ("n" %in% names(r)) as.numeric(r[["n"]]) else NA_real_,
       ev = if ("事件" %in% names(r)) as.numeric(r[["事件"]]) else NA_real_)
}
C  <- function(pat) rd1(CHARLS_RES, pat, orc = "OR_完整设计", pc = "p_完整设计")
Cf <- function(pat) rd1(CHARLS_RES, pat, orc = "OR_完整设计", loc = "下限",
                        hic = "上限", pc = "p_完整设计")
E  <- function(pat) rd1(file.path(out, "res", "01_ELSA_主结果.csv"), pat)
H  <- function(pat) rd1(file.path(out, "res", "03_HRS_主结果.csv"), pat)
S  <- function(pat) rd1(file.path(out, "res", "11_SHARE_主结果.csv"), pat)
MRo <- as.data.table(fread(MR_RES, encoding = "UTF-8"))

c_fall <- C("跌倒（主模型）"); c_inj <- C("跌倒受伤"); c_fx <- C("髋骨")
e_fall <- E("跌倒（任意）"); e_inj <- E("跌倒受伤"); e_fx <- E("髋骨")
h_fall <- H("跌倒（任意）"); h_inj <- H("跌倒受伤"); h_fx <- H("髋骨")
s_fall <- S("跌倒");        s_fx <- S("髋骨")

D <- data.table(
  证据源 = c("CHARLS 中国 2011→2018", "ELSA 英国 w2–w8", "HRS 美国 w8–w13",
             "SHARE 欧洲 w1–w6", "MR：遗传预测握力 → 跌倒（UK Biobank）"),
  类型   = c("同源发现队列", "外部复现", "外部复现", "外部复现", "因果证据"),
  握力SD = c(7.5, 6.2, 6.0, 7.1, NA),
  n      = c(c_fall$n, e_fall$n, h_fall$n, s_fall$n, NA),
  事件   = c(c_fall$ev, e_fall$ev, h_fall$ev, s_fall$ev, NA),
  OR     = c(c_fall$OR, e_fall$OR, h_fall$OR, s_fall$OR, MRo[1]$or),
  下限   = c(c_fall$L,  e_fall$L,  h_fall$L,  s_fall$L,  MRo[1]$lci),
  上限   = c(c_fall$U,  e_fall$U,  h_fall$U,  s_fall$U,  MRo[1]$uci),
  p      = c(c_fall$p,  e_fall$p,  h_fall$p,  s_fall$p,  MRo[1]$p)
)
## MHAS（墨西哥）本轮不纳入：harmonized 里握力只有 w3、非缺失仅 257 人（女性 112），
## 做不了暴露-结局配对，属数据可得性问题，不是结果问题。

## ---------- 2. 分结局三方/四方对照 ----------
CMP <- data.table(
  结局 = c(rep("跌倒（任意）", 4), rep("跌倒受伤", 3), rep("新发髋骨骨折", 4)),
  队列 = c("CHARLS", "ELSA", "HRS", "SHARE",
           "CHARLS", "ELSA", "HRS",
           "CHARLS", "ELSA", "HRS", "SHARE"),
  n    = c(c_fall$n, e_fall$n, h_fall$n, s_fall$n,
           c_inj$n,  e_inj$n,  h_inj$n,
           c_fx$n,   e_fx$n,   h_fx$n,   s_fx$n),
  事件 = c(c_fall$ev, e_fall$ev, h_fall$ev, s_fall$ev,
           c_inj$ev,  e_inj$ev,  h_inj$ev,
           c_fx$ev,   e_fx$ev,   h_fx$ev,   s_fx$ev),
  OR   = c(c_fall$OR, e_fall$OR, h_fall$OR, s_fall$OR,
           c_inj$OR,  e_inj$OR,  h_inj$OR,
           c_fx$OR,   e_fx$OR,   h_fx$OR,   s_fx$OR),
  下限 = c(c_fall$L, e_fall$L, h_fall$L, s_fall$L,
           c_inj$L,  e_inj$L,  h_inj$L,
           c_fx$L,   e_fx$L,   h_fx$L,   s_fx$L),
  上限 = c(c_fall$U, e_fall$U, h_fall$U, s_fall$U,
           c_inj$U,  e_inj$U,  h_inj$U,
           c_fx$U,   e_fx$U,   h_fx$U,   s_fx$U),
  p    = c(c_fall$p, e_fall$p, h_fall$p, s_fall$p,
           c_inj$p,  e_inj$p,  h_inj$p,
           c_fx$p,   e_fx$p,   h_fx$p,   s_fx$p))
## SHARE 无 fallinj 变量（全库检索为 0），故跌倒受伤为三队列；
## SHARE 跌倒口径为"近 6 个月"（rXfall_s），其余为 1–2 年，属口径差异而非方法差异。

## ---------- 3. 合并估计（DerSimonian–Laird 随机效应） ----------
pf <- function(OR, LO, HI) {
  b <- log(OR); se <- (log(HI) - log(LO)) / (2 * 1.96)
  w <- 1 / se^2; bf <- sum(w * b) / sum(w)
  Q <- sum(w * (b - bf)^2); df <- length(b) - 1
  C <- sum(w) - sum(w^2) / sum(w); tau2 <- max(0, (Q - df) / C)
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(or = exp(br), lo = exp(br - 1.96 * ser), hi = exp(br + 1.96 * ser),
       p = 2 * pnorm(-abs(br / ser)), I2 = I2, Qp = pchisq(Q, df, lower.tail = FALSE), n_cohort = df + 1)
}
meta_all <- function(ocn, lab) {
  d <- CMP[结局 == ocn]; m <- pf(d$OR, d$下限, d$上限)
  data.table(结局 = lab, 队列数 = m$n_cohort, 总n = sum(d$n), 总事件 = sum(d$事件),
             合并OR = m$or, 下限 = m$lo, 上限 = m$hi, p = m$p,
             Q检验p = m$Qp, I2_百分比 = m$I2)
}
M <- rbindlist(list(meta_all("跌倒（任意）", "跌倒（任意）"),
                    meta_all("跌倒受伤", "跌倒受伤"),
                    meta_all("新发髋骨骨折", "新发髋骨骨折")))
cat("\n=== 主结果（跌倒）===\n"); print(D)
cat("\n=== 分结局对照 ===\n"); print(CMP)
cat("\n=== 合并估计 ===\n"); print(M)

fwrite(D, file.path(out, "res", "07_多队列对照_跌倒.csv"))
fwrite(CMP, file.path(out, "res", "09_分结局对照.csv"))
fwrite(M, file.path(out, "res", "08_多队列合并.csv"))

## ---------- 4. 森林图（三结局 x 多队列 + 合并 + MR） ----------
PNG <- file.path(out, "res", paste0("10_森林图", FSFX, ".png"))
ragg::agg_png(PNG, width = 11.6, height = 7.0,
              units = "in", res = 300, background = "white")
## 英文标签比中文长得多（"Pooled (4 cohorts)" vs "四队列合并"），左侧页边距需相应加宽，否则会被裁切
par(mfrow = c(1, 3), mar = c(4.4, if (EN) 11.6 else 7.0, 4.2, 0.8), family = "Arial Unicode MS", las = 1)
BLUE <- "steelblue"; RED <- "darkred"; GREEN <- "forestgreen"

## 结局英文名映射（图内文字）
OMAP <- c("跌倒（任意）" = "Any fall", "跌倒受伤" = "Injurious fall", "新发髋骨骨折" = "Incident hip fracture")

plot_one <- function(ocn) {
  d <- CMP[结局 == ocn]; m <- pf(d$OR, d$下限, d$上限)
  nrow_co <- nrow(d)
  yco <- seq(nrow_co + 2, 3)          # 队列行（自上而下）
  yp  <- 2                            # 合并行
  orv <- c(d$OR, m$or); lo <- c(d$下限, m$lo); hi <- c(d$上限, m$hi)
  yy <- c(yco, yp); cols <- c(rep(BLUE, nrow_co), RED)
  lab <- c(CLAB(d$队列), TT(paste0(nrow_co, " 队列合并"), paste0("Pooled (", nrow_co, " cohorts)")))
  ylim_top <- nrow_co + 3.1
  if (ocn == "跌倒（任意）") {
    orv <- c(orv, D[类型 == "因果证据", OR]); lo <- c(lo, D[类型 == "因果证据", 下限])
    hi <- c(hi, D[类型 == "因果证据", 上限]); yy <- c(yy, 1)
    cols <- c(cols, GREEN); lab <- c(lab, TT("MR（UKB）", "MR (UK Biobank)"))
  }
  plot(NA, xlim = c(0.6, 1.5), ylim = c(0.4, ylim_top), log = "x",
       yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
  axis(1, at = c(0.7, 0.85, 1, 1.2, 1.45), labels = c("0.7", "0.85", "1.0", "1.2", "1.45"), cex.axis = 1.0)
  abline(v = 1, lty = 2, col = "grey35")
  segments(lo, yy, hi, yy, lwd = 2.2, col = cols)
  points(orv, yy, pch = 15, cex = 1.45, col = cols)
  axis(2, at = c(yy), labels = c(lab), tick = FALSE, cex.axis = if (EN) 0.92 else 1.0)
  text(orv, yy + 0.32, sprintf("%.3f (%.3f-%.3f)", orv, lo, hi), cex = 0.8, pos = 3)
  ## 注意：ragg + Arial Unicode MS 下 font=2 会吞掉数字与拉丁字母，只能用普通字重
  onm <- if (EN) unname(OMAP[ocn]) else ocn
  mtext(paste0(onm, TT("\n合并 OR ", "\nPooled OR "), sprintf("%.3f", m$or),
               "  I²=", sprintf("%.0f", m$I2), "%"), side = 3, line = 0.8, cex = 1.0)
  mtext(TT("OR（每 +1 SD 握力）", "OR per 1-SD higher grip strength"), side = 1, line = 2.4, cex = 0.95)
}
for (ocn in c("跌倒（任意）", "跌倒受伤", "新发髋骨骨折")) plot_one(ocn)
invisible(dev.off())
cat("\n森林图已输出：", PNG, "\n")
