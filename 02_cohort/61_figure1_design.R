## 81_Figure1_研究设计时间轴.R —— Figure 1：四队列研究设计时间轴
## 注意：ragg 出图时若用粗体（font=2），Arial Unicode MS 会吞掉阿拉伯数字与拉丁字母，故全程普通字重。
## 【2026-09-17】加入语言开关：中文版文章出中文图，英文版文章出英文图。
##   Rscript 81_Figure1_研究设计时间轴.R            → res/11_Figure1_研究设计时间轴.png（中文）
##   FIGLANG=en Rscript 81_Figure1_研究设计时间轴.R → res/11_Figure1_研究设计时间轴_EN.png（英文）
suppressMessages({ library(data.table); library(ragg) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")
out <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS"

## 各队列：暴露波（握力）与结局波（跌倒）年份
CO <- list(
  list(name = TT("CHARLS（中国）", "CHARLS (China)"),
       ## 【2026-09-17 删速率项后】CHARLS 主模型样本由 4,724 增至 5,007（原按三波握力速率筛过）
       n = TT("5,007 人（女性 ≥45 岁）", "5,007 women (aged ≥45 years)"),
       recall = TT("回忆窗口：1–3 年（2015、2018 两轮合并为「自上次访问以来」）",
                   "Recall: 1–3 years (the 2015 and 2018 waves combined, 'since last interview')"),
       ## 【2026-09-17 修订】主分析只用 2015 与 2018 两波合并的结局（见 82_Figure3_稳健性.R），
       ##   2013 波不进入主分析，故不再在图上画 2013 的结局标记。
       exp = c(2011), out = c(2015, 2018), use = 2018, xr = c(2011, 2018.6)),
  list(name = TT("ELSA（英国）", "ELSA (England)"),
       n = TT("8,925 人次 / 4,153 人（女性 ≥50 岁）", "8,925 observations / 4,153 women (aged ≥50 years)"),
       recall = TT("回忆窗口：2 年（w4 为 1 年，其跌倒题不用于结局）",
                   "Recall: 2 years (1 year at wave 4; that falls item was not used)"),
       exp = c(2004.5, 2008.5, 2012.5, 2016.5), out = c(2006.5, 2010.5, 2014.5, 2018.5),
       use = c(2006.5, 2010.5, 2014.5, 2018.5), xr = c(2003.9, 2019.1)),
  list(name = TT("HRS（美国）", "HRS (United States)"),
       n = TT("13,098 人次 / 7,140 人（女性 ≥50 岁）", "13,098 observations / 7,140 women (aged ≥50 years)"),
       recall = TT("回忆窗口：2 年（相邻两波配对）", "Recall: 2 years (adjacent waves paired)"),
       ## 【2026-09-17 修订】HRS 实际用的是 w8→w9 至 w13→w14 共 6 个窗口（见 20_hrs.R 第 10 行），
       ##   原图只画了 3 个，与正文「6 in HRS」及 51_逐窗口.csv 不一致，此处补齐。
       exp = c(2006, 2008, 2010, 2012, 2014, 2016), out = c(2008, 2010, 2012, 2014, 2016, 2018),
       use = c(2008, 2010, 2012, 2014, 2016, 2018), xr = c(2005.4, 2018.4)),
  list(name = TT("SHARE（欧洲）", "SHARE (Europe)"),
       n = TT("87,099 人次 / 41,074 人（女性 ≥50 岁）", "87,099 observations / 41,074 women (aged ≥50 years)"),
       recall = TT("回忆窗口：6 个月（口径较窄，事件率偏低）", "Recall: 6 months (narrower window; lower event rate)"),
       exp = c(2004, 2006.5, 2010.5, 2013, 2015), out = c(2006.5, 2010.5, 2013, 2015, 2017),
       use = c(2006.5, 2010.5, 2013, 2015, 2017), xr = c(2003.4, 2017.6))
)

PNG <- file.path(out, "res", paste0("11_Figure1_研究设计时间轴", FSFX, ".png"))
ragg::agg_png(PNG, width = 11.4, height = 6.4, units = "in", res = 320, background = "white")
par(family = "Arial Unicode MS", mar = c(3.2, 0.6, 3.0, 0.6), las = 1)
plot(NA, xlim = c(2000.0, 2019.8), ylim = c(0.35, 4.75), axes = FALSE, xlab = "", ylab = "")
abline(v = seq(2004, 2018, 2), col = "grey90", lwd = 0.8)
XL <- 2000.5            # 左侧文字栏起点
EXPC <- "#1F6FB2"; OUTC <- "#C0392B"; LINK <- "#9AA7B4"

for (i in seq_along(CO)) {
  z <- CO[[i]]; y <- 5 - i
  rng <- range(c(z$exp, z$out))
  lines(rng, c(y, y), col = "grey55", lwd = 1.4)
  ## 暴露 → 结局 连线
  for (k in seq_along(z$exp)) {
    j <- min(k, length(z$use))
    if (!is.na(z$use[j]) && z$use[j] > z$exp[k])
      segments(z$exp[k], y, z$use[j], y, col = LINK, lwd = 1.1, lty = 3)
  }
  points(z$exp, rep(y + 0.13, length(z$exp)), pch = 22, bg = EXPC, col = EXPC, cex = 1.5)
  points(z$out, rep(y - 0.13, length(z$out)), pch = 24, bg = OUTC, col = OUTC, cex = 1.5)
  text(z$exp, rep(y + 0.30, length(z$exp)), sprintf("%.1f", z$exp), cex = 0.62, col = EXPC)
  text(XL, y + 0.30, z$name, adj = 0, cex = 0.98, col = "#1F3A5F")
  text(XL, y + 0.02, z$n, adj = 0, cex = 0.70, col = "grey30")
  text(XL, y - 0.26, z$recall, adj = 0, cex = 0.66, col = "grey45")
}
axis(1, at = seq(2004, 2018, 2), labels = seq(2004, 2018, 2), cex.axis = 0.85, col = "grey40", pos = 0.35)
mtext(TT("年份", "Calendar year"), side = 1, line = 1.9, cex = 0.85)
legend("topright", bty = "n", cex = 0.85, pch = c(22, 24), pt.bg = c(EXPC, OUTC),
       col = c(EXPC, OUTC),
       legend = TT(c("握力（暴露）测量波", "跌倒（结局）判定波"),
                   c("Wave with grip strength (exposure)", "Wave with falls (outcome)")))
mtext(TT("Figure 1  四队列研究设计时间轴：暴露在前、结局在后，全部为女性前瞻性配对",
         "Figure 1  Design of the four cohort analyses: exposure always precedes outcome, all prospective pairings in women"),
      side = 3, line = 1.0, cex = 1.0, col = "#1F3A5F")
invisible(dev.off())
cat("已输出：", PNG, "\n")
