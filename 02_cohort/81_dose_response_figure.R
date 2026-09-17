## 87_RCS_图.R —— RCS 剂量反应图（分队列 kg 曲线 / 合并 z 曲线 / 阈值点图）
suppressMessages({ library(data.table); library(ggplot2); library(ragg) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")   # 中英图文字开关
RES <- "/Users/mac/Documents/Codex/2026-08-31/xian/work/外部复现-ELSA与HRS/res"
O   <- readRDS(file.path(RES, "85_RCS_对象.rds"))
TH  <- fread(file.path(RES, "81_RCS_上穿1的握力.csv"))
PC  <- fread(file.path(RES, "83_RCS_合并曲线.csv"))
NL  <- fread(file.path(RES, "84_RCS_非线性检验.csv"))

QUEUE_RAW <- c("CHARLS 中国","ELSA 英国","HRS 美国","SHARE 欧洲")
CLV <- TT(QUEUE_RAW, c("CHARLS (China)","ELSA (England)","HRS (USA)","SHARE (Europe)"))
COL <- setNames(c("#C0504D", "#2E74B5", "#4E7B4E", "#8064A2"), CLV)

## ---------- 图 1：四队列绝对握力（kg）RCS ----------
A <- rbindlist(lapply(O$A, function(r) data.table(
  队列 = r$lab, 握力kg = r$grid, OR = r$or, 下限 = r$lo, 上限 = r$hi,
  中位数 = r$medg, 阈值 = TH[队列 == r$lab, 主阈值kg_OR达1.25],
  阈值OR = TH[队列 == r$lab, 该点OR])))
## 注意：内部一律用中文原始队列名做键（names(COL) 在英文模式是英文，比较会全部落空）
A[, 队列 := factor(队列, levels = QUEUE_RAW)]
for (r in O$A) {
  A[队列 == r$lab, 面板 := TT(
    sprintf("%s（n = %s，事件 %s）\n中位数 %.1f kg；OR≥1.25 阈值 %.1f kg；P非线性 %.2f",
            r$lab, format(r$n, big.mark = ","), format(r$ev, big.mark = ","), r$medg,
            TH[队列 == r$lab, 主阈值kg_OR达1.25], r$pnl),
    sprintf("%s (n = %s; events %s)\nmedian %.1f kg; OR>=1.25 threshold %.1f kg; P non-linear %.2f",
            CLAB(r$lab), format(r$n, big.mark = ","), format(r$ev, big.mark = ","), r$medg,
            TH[队列 == r$lab, 主阈值kg_OR达1.25], r$pnl))]
}
A[, 面板 := factor(面板, levels = unique(面板[order(match(as.character(队列), QUEUE_RAW))]))]

g1 <- ggplot(A, aes(x = 握力kg, y = OR)) +
  geom_ribbon(aes(ymin = 下限, ymax = 上限), fill = "#9DB8D2", alpha = 0.35) +
  geom_line(colour = "#1F3A5F", linewidth = 0.85) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_vline(aes(xintercept = 中位数), linetype = 3, colour = "grey35", linewidth = 0.45) +
  geom_vline(xintercept = 16, linetype = 4, colour = "#B03A2E", linewidth = 0.5) +
  geom_vline(xintercept = 18, linetype = 4, colour = "#1F6FB2", linewidth = 0.5) +
  geom_point(data = A[!is.na(阈值) & abs(握力kg - 阈值) < 0.06], colour = "#B03A2E", size = 1.9) +
  facet_wrap(~ 面板, ncol = 2, scales = "free_y") +
  scale_y_log10(breaks = c(0.6, 0.8, 1.0, 1.25, 1.5, 2.0)) +
  coord_cartesian(xlim = c(8, 40)) +
  labs(title = TT("图 A  握力与跌倒的剂量反应：四队列女性限制性立方样条（RCS，4 结点）",
                  "Panel A  Dose-response between grip strength and falls: restricted cubic splines (4 knots) in women from four cohorts"),
       subtitle = TT("参照点 = 各队列握力中位数；红色虚线 = EWGSOP2 女性低握力界值 16 kg，蓝色虚线 = AWGS 2019 女性界值 18 kg；灰点线 = 中位数。\n模型统一校正年龄、教育、吸烟、高血压、糖尿病、BMI、ADL（多波次队列另校正波次固定效应），全部经抽样设计校正（svydesign + svyglm）。",
                     "Reference = cohort-specific median grip strength; red dashed line = EWGSOP2 low-grip threshold for women (16 kg); blue dashed line = AWGS 2019 threshold (18 kg); grey dotted line = median.\nAll models adjusted for age, education, smoking, hypertension, diabetes, BMI and ADL (plus wave fixed effects in multi-wave cohorts) and used design-based standard errors (svydesign + svyglm)."),
       x = TT("握力（kg）", "Grip strength (kg)"),
       y = TT("跌倒 OR（对数刻度，参照该队列中位数）", "Odds ratio for falls (log scale, referenced to the cohort median)")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13.5, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 8.7, colour = "#5A5A5A", lineheight = 1.2),
        strip.text = element_text(size = 9, colour = "#1F3A5F"),
        axis.text = element_text(size = 9), panel.grid.minor = element_blank())
agg_png(file.path(RES, paste0("90_RCS_四队列kg曲线", FSFX, ".png")), width = 2600, height = 1900, res = 200, background = "white")
print(g1); dev.off()

## ---------- 图 2：合并剂量反应（队列内 z，四队列共用结点） ----------
B <- rbindlist(lapply(O$B, function(r) data.table(
  队列 = r$lab, z = r$grid, OR = r$or, 下限 = r$lo, 上限 = r$hi)))
B <- B[z >= -3 & z <= 2]
P <- copy(PC)[z >= -3 & z <= 2]
PL <- TT("四队列合并（DL 随机效应）", "Pooled (DL random effects)")
B[, 队列 := CLAB(队列)]
LG <- data.table(队列 = PL, z = P$z, OR = P$OR, 下限 = P$下限, 上限 = P$上限)
BB <- rbind(B, LG, fill = TRUE)
BB[, 队列 := factor(队列, levels = c(CLV, PL))]
BB[, 层级 := fifelse(队列 == PL, TT("合并", "Pooled"), TT("单队列", "Single cohort"))]

g2 <- ggplot(BB, aes(x = z, y = OR, colour = 队列, fill = 队列)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey35") +
  geom_ribbon(data = BB[层级 == "合并"], aes(ymin = 下限, ymax = 上限), alpha = 0.22, colour = NA) +
  geom_line(aes(linewidth = 层级)) +
  scale_colour_manual(values = c(COL, setNames("#111111", PL))) +
  scale_fill_manual(values = c(COL, setNames("#111111", PL))) +
  guides(fill = "none", colour = guide_legend(nrow = 1)) +
  scale_linewidth_manual(values = setNames(c(0.65, 1.2), c(TT("单队列","Single cohort"), TT("合并","Pooled"))), guide = "none") +
  scale_x_continuous(breaks = -3:2, labels = c("-3", "-2", "-1", TT("0（参照）", "0 (reference)"), "+1", "+2")) +
  scale_y_log10(breaks = c(0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5)) +
  labs(title = TT("图 B  合并剂量反应曲线：握力每降低 1 个队列内标准差，跌倒风险的模式",
                  "Panel B  Pooled dose-response curve for grip strength and falls on the within-cohort standardised scale"),
       subtitle = TT(paste0("四队列共用同一组结点（z = -2.5 / -1.5 / -0.5 / +0.5 / +1.5，参照 z = 0），故样条基函数完全一致、系数可直接合并；\n",
         "合并用逐系数 DerSimonian–Laird 随机效应（对角 τ²，忽略队列内交叉协方差）。合并非线性检验 Wald = ",
         formatC(O$Wnl, format = "f", digits = 2), "，df = ", O$dfnl, "，P = ", formatC(O$Pnl, format = "f", digits = 4),
         " —— 低握力端的上升比高握力端的下降更陡。"),
         paste0("All four cohorts share the same knots (z = -2.5 / -1.5 / -0.5 / +0.5 / +1.5, reference z = 0), so the spline bases are identical and coefficients can be pooled directly.\n",
         "Pooling used coefficient-wise DerSimonian-Laird random effects (diagonal tau-squared).\nPooled test for non-linearity: Wald = ",
         formatC(O$Wnl, format = "f", digits = 2), ", df = ", O$dfnl, ", P = ", formatC(O$Pnl, format = "f", digits = 4),
         " - the rise at the low end is steeper than the fall at the high end.")),
       x = TT("队列内握力 z 值（0 = 该队列女性均值）", "Within-cohort grip strength z score (0 = cohort mean in women)"),
       y = TT("跌倒 OR（对数刻度）", "Odds ratio for falls (log scale)"), colour = NULL, fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13.5, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 8.7, colour = "#5A5A5A", lineheight = 1.2),
        axis.text = element_text(size = 9), legend.position = "bottom",
        legend.text = element_text(size = 9), panel.grid.minor = element_blank())
agg_png(file.path(RES, paste0("91_RCS_合并z曲线", FSFX, ".png")), width = 2400, height = 1550, res = 200, background = "white")
print(g2); dev.off()

## ---------- 图 3：阈值点图（kg） ----------
T2 <- copy(TH); T2[, 队列 := factor(队列, levels = rev(QUEUE_RAW))]
T2[, 队列 := factor(CLAB(as.character(队列)), levels = rev(CLV))]
T2[, 标签 := TT(sprintf("%.1f kg（OR %.2f，95%%CI %.2f–%.2f）", 主阈值kg_OR达1.25, 该点OR, 该点CI下限, 该点CI上限),
                sprintf("%.1f kg (OR %.2f, 95%%CI %.2f-%.2f)", 主阈值kg_OR达1.25, 该点OR, 该点CI下限, 该点CI上限))]
g3 <- ggplot(T2, aes(x = 主阈值kg_OR达1.25, y = 队列, colour = 队列)) +
  annotate("rect", xmin = 16, xmax = 18, ymin = -Inf, ymax = Inf, fill = "#F2C9C4", alpha = 0.45) +
  geom_segment(aes(x = 握力中位数, xend = 主阈值kg_OR达1.25, yend = 队列), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(size = 事件)) +
  geom_point(aes(x = 握力中位数), shape = 124, size = 7, colour = "grey25") +
  geom_hline(yintercept = 0.5 + seq_len(nrow(T2)), colour = "grey94") +
  geom_text(aes(x = 31.6, label = 标签), hjust = 1, size = 3.1, colour = "#333333", show.legend = FALSE) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_size_continuous(range = c(2.2, 6), guide = "none") +
  scale_x_continuous(limits = c(11, 32), breaks = seq(12, 32, 2), expand = expansion(mult = 0.015)) +
  labs(title = TT("图 C  风险开始明显升高的握力阈值（OR ≥ 1.25，95%CI 下限 > 1）",
                  "Panel C  Grip strength at which the odds of falls first exceed 1.25 with a lower confidence limit above 1"),
       subtitle = TT("圆点 = 阈值位置（相对该队列中位数 OR = 1.25）；灰色线段向右延伸至该队列握力中位数（“|”标记）；粉色带宽 = EWGSOP2 女性低握力界值区间 16–18 kg",
                     "Dot = threshold (odds ratio 1.25 relative to that cohort's median); grey segment extends to the cohort median (marked '|'); pink band = EWGSOP2 low-grip range for women (16-18 kg)"),
       x = TT("握力（kg）", "Grip strength (kg)"), y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 8.8, colour = "#5A5A5A"),
        axis.text = element_text(size = 9.5), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
agg_png(file.path(RES, paste0("92_RCS_阈值点图", FSFX, ".png")), width = 2200, height = 900, res = 200, background = "white")
print(g3); dev.off()

cat("图已生成：res/90_RCS_四队列kg曲线", FSFX, ".png, 91_RCS_合并z曲线", FSFX, ".png, 92_RCS_阈值点图", FSFX, ".png\n", sep = "")
cat("\n阈值表：\n"); print(TH)
cat("\n非线性检验：\n"); print(NL)
