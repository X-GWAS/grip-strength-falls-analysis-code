## 50_Figure4_三角验证.R —— Figure 4：观察性 / 同源 MR / 独立结局 MR 三角验证总图
## 三个面板：
##   A 森林图：四种证据层级同图对照（四队列 + 合并 / 样本重叠 MR / 独立结局 MR / 质控对照）
##   B 精度对比：95%CI 宽度（对数刻度）——重叠分析精度虚高的可视化
##   C 兼容性：独立结局 MR 的 95%CI 区间带，叠加观察性与重叠 MR 的点估计
## 所有数值一律由结果 CSV 读出后计算，不在图中手写任何数字（台账 C104）。
## 【2026-09-17】加入语言开关：中文版文章出中文图，英文版文章出英文图。
##   Rscript 50_Figure4_三角验证.R            → res/50_Figure4_三角验证总图.png（中文）
##   FIGLANG=en Rscript 50_Figure4_三角验证.R → res/50_Figure4_三角验证总图_EN.png（英文）
suppressMessages({ library(data.table); library(ggplot2); library(ragg); library(grid); library(gridExtra) })
source("/Users/mac/Documents/Codex/2026-08-31/xian/work/tools/fig_lang.R")
B    <- "/Users/mac/Documents/Codex/2026-08-31/xian/work"
RES  <- file.path(B, "外部复现-ELSA与HRS/res")
MRD  <- file.path(B, "mr-finn/res")

## ---------- 数据装配（观察性部分改由结果文件现算，不再硬编码）----------
O7 <- fread(file.path(RES, "07_多队列对照_跌倒.csv"), encoding = "UTF-8")
O8 <- fread(file.path(RES, "08_多队列合并.csv"), encoding = "UTF-8")
OC <- data.table::data.table(队列 = sub(" .*$", "", sub("^MR.*", "MR", O7$证据源)))
OBS <- data.table(
  组 = "观察性（4 国女性队列）",
  标签 = c(TT(c("CHARLS 中国", "ELSA 英国", "HRS 美国", "SHARE 欧洲"),
              c("CHARLS (China)", "ELSA (England)", "HRS (USA)", "SHARE (Europe)")),
           TT("四队列合并（DL）", "Pooled (4 cohorts, DL)")),
  OR = c(O7$OR[1:4], O8[结局 == "跌倒（任意）", 合并OR]),
  下限 = c(O7$下限[1:4], O8[结局 == "跌倒（任意）", 下限]),
  上限 = c(O7$上限[1:4], O8[结局 == "跌倒（任意）", 上限]),
  p = c(O7$p[1:4], O8[结局 == "跌倒（任意）", p]),
  主行 = c(FALSE, FALSE, FALSE, FALSE, TRUE))
OLD  <- fread(file.path(MRD, "04_对照_旧分析_UKB同源.csv"))
OVL  <- data.table(组 = "遗传学·样本重叠",
  标签 = TT(c("握力（左）→ UKB 跌倒", "握力（左）→ UKB Falling risk"),
            c("Grip strength (left) -> UKB falls", "Grip strength (left) -> UKB falling risk")),
  OR = OLD$or, 下限 = OLD$lci, 上限 = OLD$uci, p = OLD$p, 主行 = FALSE)
NEW  <- fread(file.path(MRD, "01_FinnGen_MR_主结果.csv"))[grepl("握力", 暴露) & 方法 == "Inverse variance weighted"]
IND  <- data.table(组 = "遗传学·独立结局",
  标签 = TT(c("握力（左）→ FinnGen 跌倒", "握力（右）→ FinnGen 跌倒", "低握力（EWGSOP）→ FinnGen"),
            c("Grip strength (left) -> FinnGen falls", "Grip strength (right) -> FinnGen falls",
              "Low grip strength (EWGSOP) -> FinnGen")),
  OR = NEW$or, 下限 = NEW$or_lci95, 上限 = NEW$or_uci95, p = NEW$pval, 主行 = FALSE)
## 阳/阴性对照数值取自 mr-finn 的对照结果文件（不手写，台账 C104）
FGALL <- fread(file.path(MRD, "01_FinnGen_MR_主结果.csv"))
BMI1  <- FGALL[grepl("BMI", 暴露) & 方法 == "Inverse variance weighted"][1]
if (nrow(BMI1) == 0) stop("01_FinnGen_MR_主结果.csv 中未找到 BMI 阳性对照行")
LDL1 <- fread(file.path(MRD, "08_阴性对照_LDL.csv"))[方法 == "Inverse variance weighted"][1]
CTL  <- data.table(
  组 = "质控对照",
  标签 = TT(c("BMI → FinnGen（阳性对照）", "LDL → FinnGen（阴性对照）"),
            c("BMI -> FinnGen (positive control)", "LDL -> FinnGen (negative control)")),
  OR   = c(BMI1$or, LDL1$OR), 下限 = c(BMI1$or_lci95, LDL1$下限),
  上限 = c(BMI1$or_uci95, LDL1$上限), p = c(BMI1$pval, LDL1$p), 主行 = FALSE)
GRP  <- c("观察性（4 国女性队列）", "遗传学·样本重叠", "遗传学·独立结局", "质控对照")
GRP_D <- TT(GRP, c("Observational (4 national cohorts of women)", "Genetic - overlapping samples",
                   "Genetic - independent outcome", "Positive and negative controls"))
D <- rbindlist(list(OBS, OVL, IND, CTL))
D[, 组 := factor(组, levels = GRP)]
D[, 标签 := factor(标签, levels = rev(c(OBS$标签, OVL$标签, IND$标签, CTL$标签)))]
COL <- setNames(c("#6A4C93", "#C0504D", "#2E74B5", "#4E7B4E"), GRP_D)
## 图中一切计数都由数据现算
n_row <- nrow(D); n_obs <- nrow(OBS); n_mr <- nrow(OVL) + nrow(IND); n_ctl <- nrow(CTL)

## ---------- 面板 A：森林图（单面板 + 分组标题行，避免分面条裁切）----------
SEQ <- list()
for (g in GRP) {
  SEQ[[length(SEQ) + 1]] <- data.table(表头 = TRUE, 组 = g, 标签 = paste0("| ", GRP_D[match(g, GRP)]),
    OR = NA_real_, 下限 = NA_real_, 上限 = NA_real_, p = NA_real_, 主行 = FALSE)
  SEQ[[length(SEQ) + 1]] <- copy(D[组 == g])[, 表头 := FALSE]
}
E <- rbindlist(SEQ, fill = TRUE)
E[, y := rev(seq_len(.N))]                      # 第一组显示在最上方
EH <- E[表头 == TRUE]; EB <- E[表头 == FALSE]
SEP <- EH$y - 0.5
gA <- ggplot() +
  geom_hline(yintercept = SEP, colour = "grey88", linewidth = 0.35) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35") +
  geom_errorbarh(data = EB, aes(x = OR, y = y, xmin = 下限, xmax = 上限, colour = 组, linewidth = 主行),
                 height = 0.26, show.legend = FALSE) +
  geom_point(data = EB, aes(x = OR, y = y, colour = 组, size = 主行, shape = 主行), show.legend = FALSE) +
  geom_text(data = EB, aes(x = 1.33, y = y, label = sprintf("%.3f (%.3f-%.3f)", OR, 下限, 上限)),
            hjust = 1, size = 2.6, colour = "#333333") +
  geom_text(data = EH, aes(x = 0.705, y = y, label = 标签), hjust = 0, size = 3.0, colour = "#1F3A5F") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3.0)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18)) +
  scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.1)) +
  scale_x_log10(breaks = c(0.70, 0.80, 0.90, 1.00, 1.10, 1.20), limits = c(0.70, 1.34)) +
  scale_y_continuous(breaks = EB$y, labels = EB$标签, expand = expansion(add = 0.9)) +
  labs(title = TT("A  三类证据同图对照：握力（每 +1 SD）与跌倒",
                  "A  Observational, overlapping and independent genetic estimates of the grip strength-falls association"),
       subtitle = TT(sprintf("共 %d 个估计：观察性 %d 个（发现队列 + 3 个外部复现 + 合并）、遗传学 %d 个（样本重叠 2 + 独立结局 3）、质控 %d 个。菱形 = 合并估计。\n遗传学暴露为 UK Biobank 握力（左/右 = 左手/右手）；重叠臂结局为 UKB 自报跌倒，独立臂为 FinnGen R12。\n均为两样本 MR 逆方差加权估计；观察性估计均经抽样设计校正（svydesign + svyglm）。", n_row, n_obs, n_mr, n_ctl),
                     sprintf("All %d estimates: %d observational (discovery cohort + 3 external replications + pooled), %d genetic (2 overlapping, 3 independent), and %d controls. Diamond = pooled estimate.\nGenetic exposure = UK Biobank grip strength (left/right = left/right hand); the overlapping arm uses self-reported falls in UK Biobank, the independent arm uses FinnGen R12.\nAll genetic estimates are two-sample inverse-variance weighted; all observational estimates use design-based standard errors (svydesign + svyglm).", n_row, n_obs, n_mr, n_ctl)),
       x = TT("跌倒 OR（对数刻度）", "Odds ratio for falls (log scale)"), y = NULL) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 12, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 7.4, colour = "#5A5A5A", lineheight = 1.15),
        axis.text.y = element_text(size = 7.6), axis.text.x = element_text(size = 8),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.margin = margin(t = 6, r = 26, b = 6, l = 26))

## ---------- 面板 B：精度对比（95%CI 宽度，对数刻度）----------
K <- D[标签 %in% c(OBS$标签[5], OVL$标签[1], IND$标签[1])]
K <- copy(K)
K <- K[match(c(OBS$标签[5], OVL$标签[1], IND$标签[1]), K$标签)]
KN <- TT(c("观察性：四队列合并", "遗传学：样本重叠（UKB→UKB）", "遗传学：独立结局（UKB→FinnGen）"),
         c("Observational: pooled 4 cohorts", "Genetic: overlapping (UKB to UKB)", "Genetic: independent (UKB to FinnGen)"))
K[, 名称 := KN]
K[, CI宽度 := log(上限 / 下限)]
K[, 名称 := factor(名称, levels = rev(KN))]
ratio <- K[2, CI宽度] / K[3, CI宽度]        # 重叠臂 ÷ 独立臂
gB <- ggplot(K, aes(x = CI宽度, y = 名称, fill = 名称)) +
  geom_col(width = 0.52) +
  geom_text(aes(label = TT(sprintf("95%%CI 宽度 = %.3f（对数刻度；OR 区间宽 %.2f 倍）", CI宽度, exp(CI宽度)),
                          sprintf("95%%CI width = %.3f (log scale; %.2f-fold range in OR)", CI宽度, exp(CI宽度)))),
            hjust = -0.06, size = 2.7, colour = "#333333") +
  scale_fill_manual(values = setNames(c("#6A4C93", "#C0504D", "#2E74B5"), KN), guide = "none") +
  scale_x_continuous(limits = c(0, max(K$CI宽度) * 2.9)) +
  labs(title = TT("B  精度差距：重叠分析窄得多，却不代表信息更多",
                  "B  Precision gap: the overlapping analysis is far narrower but not more informative"),
       subtitle = TT(sprintf("样本重叠臂的 95%%CI 宽度只有独立结局臂的 %.0f%%（即独立臂是它的 %.2f 倍）；两条估计在统计上相容（重叠臂点估计落在独立臂 95%%CI 内），\n提示前者是精度虚高而非额外信息。", 100 * ratio, 1 / ratio),
                     sprintf("The overlapping arm's 95%%CI is only %.0f%% as wide as the independent arm\n(which is %.2f times wider); the two estimates are statistically compatible,\nindicating inflated precision rather than extra information.",
                             100 * ratio, 1 / ratio)),
       x = TT("95%CI 宽度（对数 OR 尺度，越短越精确）", "95%CI width (log OR scale; shorter = more precise)"), y = NULL) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, colour = "#1F3A5F"),
        plot.subtitle = element_text(size = 7.2, colour = "#5A5A5A", lineheight = 1.15),
        axis.text = element_text(size = 8), panel.grid.minor = element_blank(),
        plot.margin = margin(t = 16, r = 16, b = 6, l = 6),
        panel.grid.major.y = element_blank())

## ---------- 面板 C：兼容性（独立结局 MR 的 95%CI 带）----------
ind <- IND[1]
obs <- OBS[主行 == TRUE, OR]
ovl <- OVL[1, OR]
KNM <- TT(c("观察性合并", "遗传学·样本重叠", "遗传学·独立结局"),
          c("Observational, pooled", "Genetic, overlapping", "Genetic, independent"))
gC <- ggplot() +
  annotate("rect", xmin = ind$下限, xmax = ind$上限, ymin = -Inf, ymax = Inf, fill = "#CFE0F0", alpha = 0.75) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey35") +
  geom_point(data = data.table(x = c(obs, ovl, ind$OR), y = c(3, 2, 1), 名 = KNM),
             aes(x = x, y = y, colour = 名), size = 3, show.legend = FALSE) +
  geom_text(data = data.table(x = c(obs, ovl, ind$OR), y = c(3, 2, 1),
                              名 = TT(c(sprintf("观察性合并 %.3f", obs), sprintf("样本重叠 MR %.3f", ovl), sprintf("独立结局 MR %.3f", ind$OR)),
                                      c(sprintf("Observational, pooled %.3f", obs), sprintf("Overlapping MR %.3f", ovl), sprintf("Independent MR %.3f", ind$OR)))),
            aes(x = x, y = y, label = 名), vjust = -0.9, size = 2.8, colour = "#333333") +
  scale_colour_manual(values = setNames(c("#6A4C93", "#C0504D", "#2E74B5"), KNM)) +
  scale_x_log10(breaks = c(0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10), limits = c(0.80, 1.11)) +
  scale_y_continuous(limits = c(0.5, 3.7), breaks = NULL) +
  labs(title = TT("C  兼容性：谁落在独立结局 MR 的置信区间里",
                  "C  Compatibility: which estimates fall inside the independent confidence interval"),
       subtitle = TT(sprintf("浅蓝带 = 独立结局 MR 的 95%%CI（%.3f–%.3f）。样本重叠 MR（%.3f）落在带内；观察性合并（%.3f）%s。",
                          ind$下限, ind$上限, ovl, obs,
                          ifelse(obs >= ind$下限 & obs <= ind$上限, "也落在带内",
                                 sprintf("落在带外，距下限仅 %.3f", ind$下限 - obs))),
                     sprintf("Pale blue band = 95%%CI of the independent-outcome MR estimate (%.3f-%.3f).\nThe overlapping MR estimate (%.3f) lies inside the band; the pooled observational estimate (%.3f) %s.",
                          ind$下限, ind$上限, ovl, obs,
                          ifelse(obs >= ind$下限 & obs <= ind$上限, "also lies inside the band",
                                 sprintf("lies %.3f below the lower limit", ind$下限 - obs)))),
       x = TT("跌倒 OR（对数刻度）", "Odds ratio for falls (log scale)"), y = NULL) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, colour = "#1F3A5F"),
        plot.margin = margin(t = 16, r = 16, b = 6, l = 6),
        plot.subtitle = element_text(size = 7.2, colour = "#5A5A5A", lineheight = 1.15),
        axis.text.x = element_text(size = 8), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

## ---------- 合成 ----------
right <- arrangeGrob(gB, gC, ncol = 1, heights = c(1, 1))
g <- arrangeGrob(gA, right, ncol = 2, widths = c(1.62, 1),
                 top = grid::textGrob(TT("图 4  握力与跌倒的三角验证：观察性四队列、基因预测握力（同源与独立结局）",
                                         "Figure 4  Triangulation of grip strength and falls: four observational cohorts and genetically predicted grip strength (overlapping and independent outcomes)"),
                                      gp = grid::gpar(fontsize = 13, col = "#1F3A5F"), x = 0.01, hjust = 0))
OUT <- file.path(B, "mr-finn", "res", paste0("50_Figure4_三角验证总图", FSFX, ".png"))
agg_png(OUT, width = 3400, height = 2050, res = 200, background = "white")
grid.draw(g); dev.off()

cat("图已生成:", OUT, "\n")
cat(sprintf("面板 A 行数=%d（观察性 %d、遗传学 %d、质控 %d）\n", n_row, n_obs, n_mr, n_ctl))
cat(sprintf("CI 宽度比（重叠 ÷ 独立）= %.3f / %.3f = %.3f\n", K[2, CI宽度], K[3, CI宽度], ratio))
cat(sprintf("观察性合并 %.4f；独立结局 MR 95%%CI %.4f–%.4f；落在带内? %s\n", obs, ind$下限, ind$上限,
            obs >= ind$下限 & obs <= ind$上限))
print(K[, .(名称, CI宽度 = round(CI宽度, 4))])
