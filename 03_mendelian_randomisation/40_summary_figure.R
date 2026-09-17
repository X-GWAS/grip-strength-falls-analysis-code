## 40_summary_figure.R —— 汇总表 + 对比森林图
suppressMessages({ library(data.table); library(ragg) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
setwd("/Users/mac/Documents/Codex/2026-08-31/xian/work/mr-finn")
dir.create("res", showWarnings = FALSE)

## 观察性合并估计与独立结局 MR 一律由结果文件现算（台账 C104）
POOL <- as.numeric(fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res/08_多队列合并.csv",
                         encoding = "UTF-8")[结局 == "跌倒（任意）", 合并OR][1])
FG   <- fread("res/01_FinnGen_MR_主结果.csv")[方法 == "Inverse variance weighted"]
g1 <- FG[grepl("握力（左手", 暴露)][1]; g2 <- FG[grepl("握力（右手", 暴露)][1]
g3 <- FG[grepl("低握力", 暴露)][1];       g4 <- FG[grepl("BMI", 暴露)][1]
g5 <- fread("res/08_阴性对照_LDL.csv")[方法 == "Inverse variance weighted"][1]
OL <- fread("res/04_对照_旧分析_UKB同源.csv")

d <- data.table(
  序号 = 1:8,
  分析 = TT(c(
    "观察性：四国女性队列合并（握力每 +1 SD）",
    "MR 旧：UKB 握力（左）→ UKB 跌倒（样本重叠）",
    "MR 旧：UKB 握力（左）→ UKB Falling risk（样本重叠）",
    "MR 新：UKB 握力（左）→ FinnGen 跌倒（独立）",
    "MR 新：UKB 握力（右）→ FinnGen 跌倒（独立）",
    "MR 新：低握力 EWGSOP → FinnGen 跌倒（独立）",
    "阳性对照：BMI → FinnGen 跌倒（独立）",
    "阴性对照：LDL → FinnGen 跌倒（独立）"),
    c(
    "Observational: pooled four cohorts of women (per 1-SD grip strength)",
    "Previous MR: UKB grip (left) -> UKB falls (overlapping samples)",
    "Previous MR: UKB grip (left) -> UKB falling risk (overlapping samples)",
    "New MR: UKB grip (left) -> FinnGen falls (independent)",
    "New MR: UKB grip (right) -> FinnGen falls (independent)",
    "New MR: low grip strength (EWGSOP) -> FinnGen falls (independent)",
    "Positive control: BMI -> FinnGen falls (independent)",
    "Negative control: LDL -> FinnGen falls (independent)")),
  组 = c("观察性证据", "MR 存在样本重叠", "MR 存在样本重叠",
         "MR 独立结局", "MR 独立结局", "MR 独立结局", "质控", "质控"),
  OR   = c(POOL, OL$or[1], OL$or[2], g1$or, g2$or, g3$or, g4$or, g5$OR),
  下限 = c(as.numeric(fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res/08_多队列合并.csv",
                            encoding = "UTF-8")[结局 == "跌倒（任意）", 下限][1]),
           OL$lci[1], OL$lci[2], g1$or_lci95, g2$or_lci95, g3$or_lci95, g4$or_lci95, g5$下限),
  上限 = c(as.numeric(fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res/08_多队列合并.csv",
                            encoding = "UTF-8")[结局 == "跌倒（任意）", 上限][1]),
           OL$uci[1], OL$uci[2], g1$or_uci95, g2$or_uci95, g3$or_uci95, g4$or_uci95, g5$上限),
  p = c(as.numeric(fread("/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res/08_多队列合并.csv",
                         encoding = "UTF-8")[结局 == "跌倒（任意）", p][1]),
        OL$p[1], OL$p[2], g1$pval, g2$pval, g3$pval, g4$pval, g5$p)
)
fwrite(d, "res/09_汇总表.csv")

fam <- "Arial Unicode MS"
agg_png(paste0("res/10_对比森林图", FSFX, ".png"), width = 2600, height = 1500, res = 220)
par(mar = c(5, 26, 4, 6), family = fam)
cols <- c("观察性证据" = "#7B1FA2", "MR 存在样本重叠" = "#C62828",
          "MR 独立结局" = "#1565C0", "质控" = "#2E7D32")
plot(NA, xlim = c(0.70, 1.25), ylim = c(0.5, 8.7), log = "x", yaxt = "n",
     xlab = TT("每 +1 SD 握力的跌倒 OR（对数刻度）", "Odds ratio for falls per 1-SD higher grip strength (log scale)"), ylab = "",
     main = TT("握力与跌倒：观察性、同源 MR 与独立结局 MR 的对比",
               "Grip strength and falls: observational, overlapping MR and independent-outcome MR"), cex.main = 1.25)
abline(v = 1, lty = 2, col = "grey45")
abline(v = POOL, lty = 3, col = "#7B1FA2")
y <- 8:1
segments(d$下限, y, d$上限, y, col = cols[d$组], lwd = 3)
points(d$OR, y, pch = 19, cex = 1.5, col = cols[d$组])
lab <- sprintf("%s\nOR %.3f (%.3f-%.3f)", d$分析, d$OR, d$下限, d$上限)
axis(2, at = y, labels = lab, las = 1, cex.axis = 0.72, tick = FALSE, line = -0.6)
mtext(TT(sprintf("紫色虚线 = 四队列观察性合并估计 %.3f", POOL),
         sprintf("Purple dashed line = pooled observational estimate from the four cohorts (%.3f)", POOL)),
      side = 3, line = 0.3, cex = 0.72, col = "#7B1FA2")
dev.off()
cat("图与表已写入 res/\n")
